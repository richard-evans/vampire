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

// grains module headers
#include "internal.hpp"

namespace grains{
   namespace internal{

      namespace{

         //--------------------------------------------------------------------
         // Local helpers for building clipped polygons without leaving
         // duplicate vertices at plane boundaries (see clip_halfplane()).
         //--------------------------------------------------------------------
         inline bool nearly_equal(const point2_t& a, const point2_t& b){
            const double dx = a.x - b.x;
            const double dy = a.y - b.y;
            return (dx*dx + dy*dy) < (grain_geometry_epsilon*grain_geometry_epsilon);
         }

         inline void push_unique(polygon_t& poly, const point2_t& p){
            if(poly.empty() || !nearly_equal(poly.back(), p)) poly.push_back(p);
         }

         //--------------------------------------------------------------------
         // Clips poly by num_facets evenly-spaced half-planes, each at
         // perpendicular distance rho[k] from centre in direction
         // theta_k = 2*pi*k/num_facets. The retained area shrinks
         // monotonically as every rho[k] is scaled down toward 0 together,
         // which is what round_polygon() below bisects on.
         //--------------------------------------------------------------------
         polygon_t clip_by_facet_radii(const polygon_t& poly, point2_t centre, const std::vector<double>& rho, int num_facets){

            polygon_t result = poly;
            for(int k=0; k<num_facets && !result.empty(); k++){
               const double theta = 2.0*M_PI*double(k)/double(num_facets);
               const double nx = std::cos(theta);
               const double ny = std::sin(theta);
               const double c = nx*centre.x + ny*centre.y + rho[k];
               result = clip_halfplane(result, nx, ny, c);
            }
            return result;

         }

         // Common-radius special case of clip_by_facet_radii() above, used
         // by the plain (circular) round_polygon() overload.
         polygon_t clip_by_regular_facets(const polygon_t& poly, point2_t centre, double rho, int num_facets){
            return clip_by_facet_radii(poly, centre, std::vector<double>(num_facets, rho), num_facets);
         }

      } // end of anonymous namespace

      //--------------------------------------------------------------------
      // Shoelace formula. Signed: positive for CCW input, negative for CW.
      //--------------------------------------------------------------------
      double polygon_area(const polygon_t& poly){

         const size_t n = poly.size();
         if(n < 3) return 0.0;

         double area2 = 0.0;
         for(size_t i=0; i<n; i++){
            const point2_t& p0 = poly[i];
            const point2_t& p1 = poly[(i+1)%n];
            area2 += p0.x*p1.y - p1.x*p0.y;
         }

         return 0.5*area2;

      }

      //--------------------------------------------------------------------
      // Standard polygon centroid formula. Falls back to the vertex average
      // for degenerate (near-zero-area) polygons, where the area-weighted
      // formula would divide by ~zero.
      //--------------------------------------------------------------------
      point2_t polygon_centroid(const polygon_t& poly){

         const size_t n = poly.size();
         if(n == 0) return point2_t(0.0, 0.0);
         if(n == 1) return poly[0];
         if(n == 2) return point2_t(0.5*(poly[0].x+poly[1].x), 0.5*(poly[0].y+poly[1].y));

         double area2 = 0.0;
         double cx = 0.0;
         double cy = 0.0;
         for(size_t i=0; i<n; i++){
            const point2_t& p0 = poly[i];
            const point2_t& p1 = poly[(i+1)%n];
            const double cross = p0.x*p1.y - p1.x*p0.y;
            area2 += cross;
            cx += (p0.x+p1.x)*cross;
            cy += (p0.y+p1.y)*cross;
         }

         if(std::fabs(area2) < grain_geometry_epsilon){
            double sx = 0.0, sy = 0.0;
            for(size_t i=0; i<n; i++){ sx += poly[i].x; sy += poly[i].y; }
            return point2_t(sx/double(n), sy/double(n));
         }

         const double area6 = 3.0*area2; // 6*A, since area2 = 2*A
         return point2_t(cx/area6, cy/area6);

      }

      //--------------------------------------------------------------------
      // Maximum distance from origin to any vertex of poly (0 if empty).
      //--------------------------------------------------------------------
      double max_vertex_radius(const polygon_t& poly, point2_t origin){

         double max_r2 = 0.0;
         for(size_t i=0; i<poly.size(); i++){
            const double dx = poly[i].x - origin.x;
            const double dy = poly[i].y - origin.y;
            const double r2 = dx*dx + dy*dy;
            if(r2 > max_r2) max_r2 = r2;
         }

         return std::sqrt(max_r2);

      }

      //--------------------------------------------------------------------
      // Sutherland-Hodgman clip of a single convex polygon against a single
      // half-plane nx*x + ny*y <= c. Points within grain_geometry_epsilon of
      // the plane are treated as inside, so a polygon that lies entirely on
      // the inside of the plane is returned unchanged, and clipping twice by
      // the same plane is idempotent. push_unique() suppresses the duplicate
      // vertex that a naive implementation leaves behind when an existing
      // vertex sits exactly on the clip plane.
      //--------------------------------------------------------------------
      polygon_t clip_halfplane(const polygon_t& poly, double nx, double ny, double c){

         polygon_t result;
         const size_t n = poly.size();
         if(n == 0) return result;

         for(size_t i=0; i<n; i++){

            const point2_t& curr = poly[i];
            const point2_t& prev = poly[(i + n - 1) % n];

            const double d_curr = nx*curr.x + ny*curr.y - c;
            const double d_prev = nx*prev.x + ny*prev.y - c;

            const bool curr_inside = d_curr <= grain_geometry_epsilon;
            const bool prev_inside = d_prev <= grain_geometry_epsilon;

            if(curr_inside != prev_inside){
               const double t = d_prev / (d_prev - d_curr);
               push_unique(result, point2_t(prev.x + t*(curr.x - prev.x),
                                             prev.y + t*(curr.y - prev.y)));
            }

            if(curr_inside) push_unique(result, curr);

         }

         // guard the wraparound edge: the first and last vertices pushed can
         // coincide when the polygon closes exactly on the clip plane
         if(result.size() > 1 && nearly_equal(result.front(), result.back())) result.pop_back();

         return result;

      }

      //--------------------------------------------------------------------
      // See internal.hpp: normalises orientation to CCW so that every
      // orientation-sensitive primitive below can assume it.
      //--------------------------------------------------------------------
      polygon_t ensure_ccw(const polygon_t& poly){

         if(poly.size() < 3) return poly;
         if(polygon_area(poly) < 0.0) return polygon_t(poly.rbegin(), poly.rend());
         return poly;

      }

      //--------------------------------------------------------------------
      // Moves every edge of poly inward by the perpendicular distance delta,
      // by clipping against the inward-shifted half-plane of every original
      // edge in turn. Degenerate (near-zero-length) edges are skipped rather
      // than contributing a spurious half-plane. Collapses cleanly to an
      // empty polygon when delta removes more than the polygon contains.
      //--------------------------------------------------------------------
      polygon_t offset_inward(const polygon_t& poly, double delta){

         const size_t n = poly.size();
         if(n < 3) return polygon_t();

         polygon_t result = poly;

         for(size_t i=0; i<n && !result.empty(); i++){

            const point2_t& p0 = poly[i];
            const point2_t& p1 = poly[(i+1)%n];

            const double ex = p1.x - p0.x;
            const double ey = p1.y - p0.y;
            const double len = std::sqrt(ex*ex + ey*ey);
            if(len < grain_geometry_epsilon) continue;

            // outward unit normal of edge p0->p1 for a CCW polygon
            const double nx = ey/len;
            const double ny = -ex/len;
            const double c = nx*p0.x + ny*p0.y - delta;

            result = clip_halfplane(result, nx, ny, c);

         }

         return result;

      }

      //--------------------------------------------------------------------
      // Axis-aligned rectangle [xmin,xmax] x [ymin,ymax], CCW.
      //--------------------------------------------------------------------
      polygon_t rectangle(double xmin, double ymin, double xmax, double ymax){

         polygon_t poly;
         poly.push_back(point2_t(xmin, ymin));
         poly.push_back(point2_t(xmax, ymin));
         poly.push_back(point2_t(xmax, ymax));
         poly.push_back(point2_t(xmin, ymax));

         return poly;

      }

      //--------------------------------------------------------------------
      // Bisects on the facet radius rho (see clip_by_regular_facets() above)
      // to find the radius that retains exactly
      // area_fraction * polygon_area(poly), then returns the clip at that
      // radius.
      //--------------------------------------------------------------------
      polygon_t round_polygon(const polygon_t& poly, point2_t centre, double area_fraction, int num_facets){

         if(poly.size() < 3) return poly;
         if(area_fraction >= 1.0) return poly;
         if(area_fraction <= 0.0) return polygon_t();

         const double A = polygon_area(poly);
         if(!(A > 0.0)) return poly;

         const double target = area_fraction*A;

         double rho_lo = 0.0;
         double rho_hi = max_vertex_radius(poly, centre);

         for(int iter=0; iter<60; iter++){
            const double rho_mid = 0.5*(rho_lo + rho_hi);
            const double area_mid = polygon_area(clip_by_regular_facets(poly, centre, rho_mid, num_facets));
            if(area_mid < target) rho_lo = rho_mid; else rho_hi = rho_mid;
         }

         return clip_by_regular_facets(poly, centre, rho_hi, num_facets);

      }

      //--------------------------------------------------------------------
      // Organic-rounding overload: see internal.hpp. Draws a random
      // low-order cosine series shape(theta_k), normalised so
      // shape(theta) in [1-roughness_amplitude, 1+roughness_amplitude], then
      // bisects an overall scale s so that clipping at rho[k] = s*shape[k]
      // retains exactly area_fraction*polygon_area(poly).
      //--------------------------------------------------------------------
      polygon_t round_polygon(const polygon_t& poly, point2_t centre, double area_fraction,
                               double roughness_amplitude, int roughness_modes, MTRand& rng,
                               int num_facets){

         if(roughness_amplitude <= 0.0 || roughness_modes <= 0){
            return round_polygon(poly, centre, area_fraction, num_facets);
         }

         if(poly.size() < 3) return poly;
         if(area_fraction >= 1.0) return poly;
         if(area_fraction <= 0.0) return polygon_t();

         const double A = polygon_area(poly);
         if(!(A > 0.0)) return poly;

         const double target = area_fraction*A;

         // Per-mode amplitude tapered as 1/m (so low modes dominate, giving
         // a few broad organic lobes rather than high-frequency wobble) and
         // a random phase, then L1-normalised so shape(theta) stays within
         // [1-roughness_amplitude, 1+roughness_amplitude] for every theta.
         std::vector<double> amp(roughness_modes+1, 0.0), phase(roughness_modes+1, 0.0);
         double norm = 0.0;
         for(int m=1; m<=roughness_modes; m++){
            amp[m] = (2.0*rng() - 1.0)/double(m);
            phase[m] = 2.0*M_PI*rng();
            norm += std::fabs(amp[m]);
         }
         if(norm < 1.0e-12) norm = 1.0; // degenerate all-zero draw: shape stays flat (==1)

         std::vector<double> shape(num_facets);
         for(int k=0; k<num_facets; k++){
            const double theta = 2.0*M_PI*double(k)/double(num_facets);
            double h = 0.0;
            for(int m=1; m<=roughness_modes; m++) h += amp[m]*std::cos(double(m)*theta + phase[m]);
            shape[k] = 1.0 + roughness_amplitude*(h/norm);
         }

         const double min_shape = *std::min_element(shape.begin(), shape.end());

         double s_lo = 0.0;
         // Scaling the plain overload's upper bound by 1/min_shape
         // guarantees nothing is clipped at s_hi regardless of how shape
         // dips below 1 in some directions.
         double s_hi = max_vertex_radius(poly, centre) / std::max(min_shape, 1.0e-6);

         std::vector<double> rho(num_facets);
         for(int iter=0; iter<60; iter++){
            const double s_mid = 0.5*(s_lo + s_hi);
            for(int k=0; k<num_facets; k++) rho[k] = s_mid*shape[k];
            const double area_mid = polygon_area(clip_by_facet_radii(poly, centre, rho, num_facets));
            if(area_mid < target) s_lo = s_mid; else s_hi = s_mid;
         }

         for(int k=0; k<num_facets; k++) rho[k] = s_hi*shape[k];
         return clip_by_facet_radii(poly, centre, rho, num_facets);

      }

   } // end of internal namespace
} // end of grains namespace
