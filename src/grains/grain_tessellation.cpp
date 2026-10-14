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
#include <cstdlib>

// grains module headers
#include "internal.hpp"

namespace grains{
   namespace internal{

      namespace{

         // Uniform background grid, used to walk sites outward from a query
         // site ring by ring, nearest first.
         struct background_grid_t{

            double x0, y0;
            double cell_size;
            int nx, ny;
            std::vector<std::vector<size_t> > bins;

            background_grid_t(const std::vector<point2_t>& sites,
                               double xmin, double ymin, double xmax, double ymax){

               const double area = std::max(1.0e-12, (xmax-xmin)*(ymax-ymin));
               const double span = std::max(xmax-xmin, ymax-ymin);
               const size_t n = sites.size();

               // aim for roughly one site per grid cell on average
               cell_size = (n > 0) ? std::sqrt(area/double(n)) : span;
               if(!(cell_size > 1.0e-9)) cell_size = (span > 1.0e-9) ? span : 1.0;

               x0 = xmin;
               y0 = ymin;
               nx = std::max(1, int(std::ceil((xmax-xmin)/cell_size)));
               ny = std::max(1, int(std::ceil((ymax-ymin)/cell_size)));

               bins.resize(size_t(nx)*size_t(ny));
               for(size_t i=0; i<n; i++){
                  const int gx = cell_index_x(sites[i].x);
                  const int gy = cell_index_y(sites[i].y);
                  bins[size_t(gy)*size_t(nx)+size_t(gx)].push_back(i);
               }

            }

            int cell_index_x(double x) const{
               int gx = int(std::floor((x-x0)/cell_size));
               if(gx < 0) gx = 0;
               if(gx >= nx) gx = nx-1;
               return gx;
            }

            int cell_index_y(double y) const{
               int gy = int(std::floor((y-y0)/cell_size));
               if(gy < 0) gy = 0;
               if(gy >= ny) gy = ny-1;
               return gy;
            }

            const std::vector<size_t>& bin(int gx, int gy) const{
               return bins[size_t(gy)*size_t(nx)+size_t(gx)];
            }

         };

      } // end of anonymous namespace

      // Builds every site's power cell by starting from the domain rectangle
      // and clipping against the radical axis of each neighbour in turn,
      // nearest rings first, stopping early once a ring can no longer reach
      // inside the (already partially clipped) cell.
      std::vector<polygon_t> build_power_cells(
         const std::vector<point2_t>& sites,
         const std::vector<double>&   weights,
         double xmin, double ymin, double xmax, double ymax,
         bool periodic){

         const size_t n = sites.size();
         std::vector<polygon_t> cells(n);
         if(n == 0) return cells;

         if(!periodic){

            const background_grid_t grid(sites, xmin, ymin, xmax, ymax);

            for(size_t i=0; i<n; i++){

               polygon_t cell = rectangle(xmin, ymin, xmax, ymax);

               const int gx = grid.cell_index_x(sites[i].x);
               const int gy = grid.cell_index_y(sites[i].y);
               const int max_ring = std::max(grid.nx, grid.ny);

               const double pix2 = sites[i].x*sites[i].x + sites[i].y*sites[i].y - weights[i];

               for(int ring=0; ring<=max_ring && !cell.empty(); ring++){

                  if(ring > 0){
                     const double ring_min_dist = double(ring-1)*grid.cell_size;
                     if(ring_min_dist > 2.0*max_vertex_radius(cell, sites[i])) break;
                  }

                  for(int dy=-ring; dy<=ring && !cell.empty(); dy++){
                     const int cy = gy+dy;
                     if(cy < 0 || cy >= grid.ny) continue;

                     for(int dx=-ring; dx<=ring; dx++){
                        if(std::max(std::abs(dx), std::abs(dy)) != ring) continue; // ring boundary only
                        const int cx = gx+dx;
                        if(cx < 0 || cx >= grid.nx) continue;

                        const std::vector<size_t>& candidates = grid.bin(cx, cy);
                        for(size_t k=0; k<candidates.size(); k++){

                           const size_t j = candidates[k];
                           if(j == i) continue;

                           const double nx_ = 2.0*(sites[j].x - sites[i].x);
                           const double ny_ = 2.0*(sites[j].y - sites[i].y);
                           const double pjx2 = sites[j].x*sites[j].x + sites[j].y*sites[j].y - weights[j];
                           const double c = pjx2 - pix2;

                           cell = clip_halfplane(cell, nx_, ny_, c);
                           if(cell.empty()) break;
                        }
                        if(cell.empty()) break;
                     }
                  }
               }

               cells[i] = cell;

            }

            return cells;

         }

         // Periodic path: re-runs the same ring-walk against every site
         // replicated over the 3x3 tiling of domain translations, so a
         // periodic cell can extend beyond the domain box and wrap around
         // to meet neighbours (or their images) on the far side.
         const double Lx = xmax - xmin;
         const double Ly = ymax - ymin;
         const double tx_list[3] = { -Lx, 0.0, Lx };
         const double ty_list[3] = { -Ly, 0.0, Ly };
         const size_t self_offset = 4; // (a,b) = (1,1) i.e. (tx,ty) = (0,0) in the loop below

         struct image_t{ point2_t p; double w; };
         std::vector<image_t> images;
         images.reserve(n*9);
         for(size_t i=0; i<n; i++){
            for(int a=0; a<3; a++){
               for(int b=0; b<3; b++){
                  images.push_back(image_t{ point2_t(sites[i].x+tx_list[a], sites[i].y+ty_list[b]), weights[i] });
               }
            }
         }

         std::vector<point2_t> image_points(images.size());
         for(size_t k=0; k<images.size(); k++) image_points[k] = images[k].p;

         const background_grid_t grid(image_points, xmin-Lx, ymin-Ly, xmax+Lx, ymax+Ly);

         for(size_t i=0; i<n; i++){

            const size_t self_idx = i*9 + self_offset;

            // large enough that the true (minimum-image-bounded) cell can
            // never reach its edges
            polygon_t cell = rectangle(sites[i].x-1.5*Lx, sites[i].y-1.5*Ly, sites[i].x+1.5*Lx, sites[i].y+1.5*Ly);

            const int gx = grid.cell_index_x(sites[i].x);
            const int gy = grid.cell_index_y(sites[i].y);
            const int max_ring = std::max(grid.nx, grid.ny);

            const double pix2 = sites[i].x*sites[i].x + sites[i].y*sites[i].y - weights[i];

            for(int ring=0; ring<=max_ring && !cell.empty(); ring++){

               if(ring > 0){
                  const double ring_min_dist = double(ring-1)*grid.cell_size;
                  if(ring_min_dist > 2.0*max_vertex_radius(cell, sites[i])) break;
               }

               for(int dy=-ring; dy<=ring && !cell.empty(); dy++){
                  const int cy = gy+dy;
                  if(cy < 0 || cy >= grid.ny) continue;

                  for(int dx=-ring; dx<=ring; dx++){
                     if(std::max(std::abs(dx), std::abs(dy)) != ring) continue; // ring boundary only
                     const int cx = gx+dx;
                     if(cx < 0 || cx >= grid.nx) continue;

                     const std::vector<size_t>& candidates = grid.bin(cx, cy);
                     for(size_t k=0; k<candidates.size(); k++){

                        const size_t j = candidates[k];
                        if(j == self_idx) continue;

                        const double nx_ = 2.0*(images[j].p.x - sites[i].x);
                        const double ny_ = 2.0*(images[j].p.y - sites[i].y);
                        const double pjx2 = images[j].p.x*images[j].p.x + images[j].p.y*images[j].p.y - images[j].w;
                        const double c = pjx2 - pix2;

                        cell = clip_halfplane(cell, nx_, ny_, c);
                        if(cell.empty()) break;
                     }
                     if(cell.empty()) break;
                  }
               }
            }

            cells[i] = cell;

         }

         return cells;

      }

      // Adapts build_power_cells() to the grain_coord_array/
      // grain_vertices_array format used throughout this module.
      // grain_coord_array[i] is overwritten with the built cell's own
      // centroid on exit. use_laguerre_weights selects real per-seed
      // weights vs uniform (plain Voronoi) weights.
      void populate_vertex_points_power(std::vector <std::vector <double> > & grain_coord_array,
                                         std::vector <std::vector <std::vector <double> > > &  grain_vertices_array,
                                         const std::vector<double>& weights_in,
                                         double domain_x, double domain_y,
                                         bool use_laguerre_weights,
                                         bool periodic){

         const size_t n = grain_coord_array.size();

         std::vector<point2_t> sites(n);
         for(size_t i=0; i<n; i++) sites[i] = point2_t(grain_coord_array[i][0], grain_coord_array[i][1]);

         std::vector<double> weights(n, 1.0);
         if(use_laguerre_weights && weights_in.size() == n) weights = weights_in;

         const std::vector<polygon_t> cells = build_power_cells(sites, weights, 0.0, 0.0, domain_x, domain_y, periodic);

         grain_vertices_array.assign(n, std::vector<std::vector<double> >());

         for(size_t i=0; i<n; i++){

            const polygon_t& cell = cells[i];
            if(cell.size() < 3){
               // vanished cell (weight too small relative to neighbours): leave empty to exclude this grain
               grain_coord_array[i][0] = 0.0;
               grain_coord_array[i][1] = 0.0;
               continue;
            }

            for(size_t v=0; v<cell.size(); v++){
               std::vector<double> vertex(2);
               vertex[0] = cell[v].x;
               vertex[1] = cell[v].y;
               grain_vertices_array[i].push_back(vertex);
            }

            const point2_t centroid = polygon_centroid(cell);
            grain_coord_array[i][0] = centroid.x;
            grain_coord_array[i][1] = centroid.y;
         }

      }

      //--------------------------------------------------------------------
      int remove_boundary_grains(std::vector<std::vector<double> >& grain_coord_array,
                                  std::vector<std::vector<std::vector<double> > >& grain_vertices_array,
                                  double xmin, double ymin, double xmax, double ymax){

         const double eps = 1.0e-6 * std::max(xmax-xmin, ymax-ymin);
         int removed = 0;

         for(size_t g = 0; g < grain_vertices_array.size(); g++){

            const size_t nv = grain_vertices_array[g].size();
            if(nv < 3) continue; // already empty (e.g. weight-fitting vanished)

            bool touches_boundary = false;
            for(size_t v = 0; v < nv && !touches_boundary; v++){
               const double x = grain_vertices_array[g][v][0];
               const double y = grain_vertices_array[g][v][1];
               if(std::fabs(x-xmin) < eps || std::fabs(x-xmax) < eps ||
                  std::fabs(y-ymin) < eps || std::fabs(y-ymax) < eps){
                  touches_boundary = true;
               }
            }

            if(touches_boundary){
               grain_vertices_array[g].clear();
               grain_coord_array[g][0] = 0.0;
               grain_coord_array[g][1] = 0.0;
               removed++;
            }

         }

         return removed;

      }

      // Carves out a boundary gap of width delta between grains by
      // offsetting every grain's cell inward by delta, giving a constant
      // physical gap width regardless of grain size.
      void apply_grain_spacing(std::vector <std::vector <double> > & grain_coord_array,
                                std::vector <std::vector <std::vector <double> > > &  grain_vertices_array,
                                const std::vector<double>& delta){

         for(size_t grain=0; grain<grain_vertices_array.size(); grain++){

            const size_t nv = grain_vertices_array[grain].size();
            if(nv == 0) continue;

            polygon_t cell(nv);
            for(size_t v=0; v<nv; v++){
               cell[v] = point2_t(grain_vertices_array[grain][v][0], grain_vertices_array[grain][v][1]);
            }
            cell = ensure_ccw(cell);

            const polygon_t offset = offset_inward(cell, delta[grain]);

            const double x0 = grain_coord_array[grain][0];
            const double y0 = grain_coord_array[grain][1];

            grain_vertices_array[grain].assign(offset.size(), std::vector<double>(2));
            for(size_t v=0; v<offset.size(); v++){
               grain_vertices_array[grain][v][0] = offset[v].x - x0;
               grain_vertices_array[grain][v][1] = offset[v].y - y0;
            }

         }

      }

      // scalar-delta overload: every grain gets the same offset
      void apply_grain_spacing(std::vector <std::vector <double> > & grain_coord_array,
                                std::vector <std::vector <std::vector <double> > > &  grain_vertices_array,
                                double delta){

         apply_grain_spacing(grain_coord_array, grain_vertices_array,
                              std::vector<double>(grain_vertices_array.size(), delta));

      }

   } // end of internal namespace
} // end of grains namespace
