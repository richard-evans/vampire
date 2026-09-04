//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <cmath>

// Vampire headers
#include "random.hpp"

// grains module headers
#include "internal.hpp"

namespace grains{
   namespace internal{

      namespace{

         //--------------------------------------------------------------------
         // Builds the free-standing polygon bounded by num_facets half-planes
         // with normals (nx[k],ny[k]) at perpendicular distance h[k] from
         // centre, by clipping a generous starting square. Used only to
         // bisect the Wulff scale s in apply_wulff_facets() below, to find
         // the facet distances that give the requested target area - the
         // faceted GRAIN itself is produced separately, by clipping the
         // original cell (see that function's own comment).
         //--------------------------------------------------------------------
         polygon_t facet_bound_polygon(point2_t centre, double half_width,
                                        const std::vector<double>& nx, const std::vector<double>& ny,
                                        const std::vector<double>& h){

            polygon_t result = rectangle(centre.x-half_width, centre.y-half_width,
                                          centre.x+half_width, centre.y+half_width);

            for(size_t k=0; k<h.size() && !result.empty(); k++){
               const double c = nx[k]*centre.x + ny[k]*centre.y + h[k];
               result = clip_halfplane(result, nx[k], ny[k], c);
            }

            return result;

         }

      } // end of anonymous namespace

      //--------------------------------------------------------------------
      // See internal.hpp for the full description of the blend. h_cell,k is
      // poly's own exact support distance from centre in direction k, so
      // clipping at exactly this distance removes nothing (strength=0 is an
      // exact no-op). h_wulff,k = s*gamma_k, with gamma_k the
      // alternating-family anisotropy factor and s found by bisection
      // (grain_geometry.cpp's round_polygon() uses the same approach) to
      // hit target_area.
      //--------------------------------------------------------------------
      polygon_t apply_wulff_facets(const polygon_t& poly, point2_t centre,
                                    int num_facets, double theta0,
                                    double anisotropy, double strength,
                                    double target_area){

         if(num_facets <= 0) return poly;
         if(strength <= 0.0) return poly;
         if(poly.size() < 3) return poly;
         if(!(target_area > 0.0)) return poly;

         std::vector<double> nx(num_facets), ny(num_facets), gamma(num_facets), h_cell(num_facets);

         for(int k=0; k<num_facets; k++){

            const double theta = theta0 + 2.0*M_PI*double(k)/double(num_facets);
            nx[k] = std::cos(theta);
            ny[k] = std::sin(theta);
            gamma[k] = (k % 2 == 0) ? 1.0 : anisotropy;

            double h = 0.0;
            for(size_t v=0; v<poly.size(); v++){
               const double proj = nx[k]*(poly[v].x-centre.x) + ny[k]*(poly[v].y-centre.y);
               if(proj > h) h = proj;
            }
            h_cell[k] = h;

         }

         // bracket and bisect the Wulff scale s so that the free-standing
         // facet polygon at h_wulff,k = s*gamma_k has area == target_area
         const double half_width = 4.0*max_vertex_radius(poly, centre) + 1.0;

         double s_lo = 0.0;
         double s_hi = half_width;
         while(polygon_area(facet_bound_polygon(centre, half_width, nx, ny,
                  [&](){ std::vector<double> h(num_facets); for(int k=0;k<num_facets;k++) h[k]=s_hi*gamma[k]; return h; }()))
               < target_area && s_hi < 1.0e6*half_width){
            s_hi *= 2.0;
         }

         for(int iter=0; iter<60; iter++){
            const double s_mid = 0.5*(s_lo + s_hi);
            std::vector<double> h_mid(num_facets);
            for(int k=0; k<num_facets; k++) h_mid[k] = s_mid*gamma[k];
            const double area_mid = polygon_area(facet_bound_polygon(centre, half_width, nx, ny, h_mid));
            if(area_mid < target_area) s_lo = s_mid; else s_hi = s_mid;
         }

         // blend and clip the ORIGINAL cell - this guarantees the faceted
         // result is always contained in poly and its area never exceeds
         // poly's, since every step below is a clip against poly, never a
         // free-standing construction of the Wulff shape.
         polygon_t result = poly;
         for(int k=0; k<num_facets && !result.empty(); k++){
            const double h_wulff = s_hi*gamma[k];
            const double h_blend = (1.0-strength)*h_cell[k] + strength*h_wulff;
            const double c = nx[k]*centre.x + ny[k]*centre.y + h_blend;
            result = clip_halfplane(result, nx[k], ny[k], c);
         }

         return result;

      }

      //--------------------------------------------------------------------
      // See internal.hpp. Must run before the boundary-spacing offset,
      // so grain_vertices_array is still in ABSOLUTE coordinates here.
      //--------------------------------------------------------------------
      void apply_grain_faceting(std::vector<std::vector<double> >& grain_coord_array,
                                 std::vector<std::vector<std::vector<double> > >& grain_vertices_array,
                                 const std::vector<double>& orientation,
                                 int num_facets, double anisotropy, double strength){

         if(num_facets <= 0 || strength <= 0.0) return;

         for(size_t grain=0; grain<grain_vertices_array.size(); grain++){

            const size_t nv = grain_vertices_array[grain].size();
            if(nv == 0) continue;

            polygon_t cell(nv);
            for(size_t v=0; v<nv; v++){
               cell[v] = point2_t(grain_vertices_array[grain][v][0], grain_vertices_array[grain][v][1]);
            }
            // apply_wulff_facets() assumes CCW vertex order; fix it up if needed
            cell = ensure_ccw(cell);

            const double theta0 = (grain < orientation.size()) ? orientation[grain] : 0.0;
            const point2_t centre(grain_coord_array[grain][0], grain_coord_array[grain][1]);

            // target the cell's own area, so faceting reshapes the grain
            // without shrinking the weight-fitted size distribution
            const double target_area = polygon_area(cell);

            const polygon_t faceted = apply_wulff_facets(cell, centre, num_facets, theta0, anisotropy, strength, target_area);

            // faceting can legitimately empty a cell at extreme strength/
            // anisotropy combinations (a thin sliver clipped away entirely);
            // zero vertices is the same "exclude this grain" signal used by
            // tessellation and spacing.
            grain_vertices_array[grain].assign(faceted.size(), std::vector<double>(2));
            for(size_t v=0; v<faceted.size(); v++){
               grain_vertices_array[grain][v][0] = faceted[v].x;
               grain_vertices_array[grain][v][1] = faceted[v].y;
            }

         }

      }

      //--------------------------------------------------------------------
      // See internal.hpp.
      //--------------------------------------------------------------------
      std::vector<double> generate_grain_facet_orientations(int num_grains, int num_facets,
                                                              grain_facet_orientation_t mode,
                                                              double angle_degrees, double spread_degrees,
                                                              MTRand& rng){

         std::vector<double> theta(std::max(0, num_grains), 0.0);
         if(num_facets <= 0) return theta;

         const double angle0 = angle_degrees * M_PI/180.0;
         const double spread = spread_degrees * M_PI/180.0;
         const double period = 2.0*M_PI/double(num_facets);

         for(size_t i=0; i<theta.size(); i++){
            switch(mode){
               case facet_orientation_fixed:
                  theta[i] = angle0;
                  break;
               case facet_orientation_textured:
                  theta[i] = angle0 + spread*mtrandom::gaussianc(rng);
                  break;
               case facet_orientation_random:
               default:
                  theta[i] = period*rng();
                  break;
            }
         }

         return theta;

      }

   } // end of internal namespace
} // end of grains namespace
