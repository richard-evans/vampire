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

         // Uniform background grid for overlap-force queries during
         // relaxation; rebuilt from scratch each step, queried per disc for
         // its 3x3 neighbourhood. periodic wraps at the domain edge.
         class pack_grid_t{

         public:

            pack_grid_t(const double xmin, const double ymin,
                        const double xmax, const double ymax,
                        const double cell_size, const bool periodic):
               xmin_(xmin), ymin_(ymin), cell_size_(cell_size), periodic_(periodic)
            {
               nx_ = std::max(1, int(std::ceil((xmax - xmin) / cell_size)));
               ny_ = std::max(1, int(std::ceil((ymax - ymin) / cell_size)));
               cells_.resize(size_t(nx_) * size_t(ny_));
            }

            void clear(){
               for(size_t c = 0; c < cells_.size(); c++) cells_[c].clear();
            }

            void insert(const double x, const double y, const int index){
               cells_[cell_index(cell_x(x), cell_y(y))].push_back(index);
            }

            // calls visitor(j) once for every disc index registered in the
            // 3x3 block of cells centred on (x,y)
            template<class F>
            void for_each_candidate(const double x, const double y, F visitor) const {

               const int cx = cell_x(x);
               const int cy = cell_y(y);

               for(int dj = -1; dj <= 1; dj++){
                  int j = cy + dj;
                  if(periodic_) j = ((j % ny_) + ny_) % ny_;
                  else if(j < 0 || j >= ny_) continue;

                  for(int di = -1; di <= 1; di++){
                     int i = cx + di;
                     if(periodic_) i = ((i % nx_) + nx_) % nx_;
                     else if(i < 0 || i >= nx_) continue;

                     const std::vector<int>& cell = cells_[cell_index(i, j)];
                     for(size_t k = 0; k < cell.size(); k++) visitor(cell[k]);
                  }
               }
            }

         private:

            // clamps to a safe finite range before converting to int
            int cell_x(const double x) const {
               double rel = (x - xmin_) / cell_size_;
               if(!std::isfinite(rel)) rel = 0.0;
               rel = std::max(-1.0e9, std::min(1.0e9, rel));
               int i = int(std::floor(rel));
               if(periodic_) i = ((i % nx_) + nx_) % nx_;
               else { if(i < 0) i = 0; if(i >= nx_) i = nx_ - 1; }
               return i;
            }

            int cell_y(const double y) const {
               double rel = (y - ymin_) / cell_size_;
               if(!std::isfinite(rel)) rel = 0.0;
               rel = std::max(-1.0e9, std::min(1.0e9, rel));
               int j = int(std::floor(rel));
               if(periodic_) j = ((j % ny_) + ny_) % ny_;
               else { if(j < 0) j = 0; if(j >= ny_) j = ny_ - 1; }
               return j;
            }

            int cell_index(const int i, const int j) const { return i + j * nx_; }

            double xmin_, ymin_, cell_size_;
            bool periodic_;
            int nx_, ny_;
            std::vector< std::vector<int> > cells_;

         };

         // shortest x-displacement from xi to xj under minimum-image
         // periodic wrapping of a domain of length L; a no-op when !periodic
         inline double min_image(double d, const double L, const bool periodic){
            if(periodic) d -= L * std::round(d / L);
            return d;
         }

         // Soft pairwise repulsion between overlapping discs, proportional to
         // overlap depth. Non-periodic mode also repels discs from the
         // current box walls, herding them inward; periodic mode uses
         // disc-disc forces alone, evaluated at minimum image.
         void calculate_forces(const std::vector<double>& x, const std::vector<double>& y,
                              const std::vector<double>& radius,
                              std::vector<double>& force_x, std::vector<double>& force_y,
                              const double xmin, const double ymin, const double xmax, const double ymax,
                              const bool periodic, const double force_const){

            const size_t n = x.size();
            const double Lx = xmax - xmin;
            const double Ly = ymax - ymin;

            double max_d = 0.0;
            for(size_t i = 0; i < n; i++) max_d = std::max(max_d, 2.0*radius[i]);
            const double cell_size = std::max(max_d, 1.0e-6);

            pack_grid_t grid(xmin, ymin, xmax, ymax, cell_size, periodic);
            for(size_t i = 0; i < n; i++) grid.insert(x[i], y[i], int(i));

            for(size_t i = 0; i < n; i++){

               double fx = 0.0, fy = 0.0;
               const double ri = radius[i];

               grid.for_each_candidate(x[i], y[i], [&](int j){
                  if(size_t(j) == i) return;
                  double dx = x[j] - x[i];
                  double dy = y[j] - y[i];
                  if(periodic){
                     dx = min_image(dx, Lx, true);
                     dy = min_image(dy, Ly, true);
                  }
                  const double min_dist = ri + radius[j];
                  const double dist2 = dx*dx + dy*dy;
                  if(dist2 < min_dist*min_dist){
                     const double dist = std::sqrt(dist2);
                     const double overlap = min_dist - dist;
                     double ux, uy; // unit vector from i toward j
                     if(dist > 1.0e-9*min_dist){
                        ux = dx/dist;
                        uy = dy/dist;
                     }
                     else{
                        // near-exact coincidence: deterministic push direction from the pair's indices
                        const double angle = std::fmod(12.9898*double(std::min(i, size_t(j))) +
                                                        78.233*double(std::max(i, size_t(j))), 2.0*M_PI);
                        const double sign = (size_t(j) > i) ? 1.0 : -1.0;
                        ux = sign * std::cos(angle);
                        uy = sign * std::sin(angle);
                     }
                     fx += -force_const * overlap * ux;
                     fy += -force_const * overlap * uy;
                  }
               });

               if(!periodic){
                  if(x[i] - xmin < ri) fx += force_const * (ri - (x[i] - xmin));
                  if(xmax - x[i] < ri) fx -= force_const * (ri - (xmax - x[i]));
                  if(y[i] - ymin < ri) fy += force_const * (ri - (y[i] - ymin));
                  if(ymax - y[i] < ri) fy -= force_const * (ri - (ymax - y[i]));
               }

               force_x[i] = fx;
               force_y[i] = fy;

            }

         }

         // overdamped integration step, capping displacement to avoid overshoot
         void move_particles(std::vector<double>& x, std::vector<double>& y,
                             const std::vector<double>& fx, const std::vector<double>& fy,
                             const double dt, const double max_step){
            const size_t n = x.size();
            for(size_t i = 0; i < n; i++){
               double dx = fx[i] * dt;
               double dy = fy[i] * dt;
               const double d = std::sqrt(dx*dx + dy*dy);
               if(d > max_step){
                  const double scale = max_step / d;
                  dx *= scale;
                  dy *= scale;
               }
               x[i] += dx;
               y[i] += dy;
            }
         }

         // wraps every disc centre back into [0,Lx) x [0,Ly)
         void wrap_positions(std::vector<double>& x, std::vector<double>& y,
                             const double Lx, const double Ly){
            const size_t n = x.size();
            for(size_t i = 0; i < n; i++){
               x[i] = std::fmod(x[i], Lx); if(x[i] < 0.0) x[i] += Lx;
               y[i] = std::fmod(y[i], Ly); if(y[i] < 0.0) y[i] += Ly;
            }
         }

         // compresses disc positions and box extent toward (target_x,target_y) by one step
         void shrink(std::vector<double>& x, std::vector<double>& y,
                    double& Lx, double& Ly,
                    const double target_x, const double target_y,
                    const int steps_remaining){
            const double dLx = (Lx - target_x) / double(steps_remaining);
            const double dLy = (Ly - target_y) / double(steps_remaining);
            const double fx = (Lx - dLx) / Lx;
            const double fy = (Ly - dLy) / Ly;
            const size_t n = x.size();
            for(size_t i = 0; i < n; i++){
               x[i] *= fx;
               y[i] *= fy;
            }
            Lx -= dLx;
            Ly -= dLy;
         }

         // largest fractional penetration depth 1 - dist/(r_i+r_j) among any
         // currently-overlapping pair, 0.0 if every disc is clear of every other
         double max_overlap_fraction(const std::vector<double>& x, const std::vector<double>& y,
                                     const std::vector<double>& radius,
                                     const double xmin, const double ymin, const double xmax, const double ymax,
                                     const bool periodic){

            const size_t n = x.size();
            if(n == 0) return 0.0;

            const double Lx = xmax - xmin;
            const double Ly = ymax - ymin;

            double max_d = 0.0;
            for(size_t i = 0; i < n; i++) max_d = std::max(max_d, 2.0*radius[i]);
            const double cell_size = std::max(max_d, 1.0e-6);

            pack_grid_t grid(xmin, ymin, xmax, ymax, cell_size, periodic);
            for(size_t i = 0; i < n; i++) grid.insert(x[i], y[i], int(i));

            double worst = 0.0;

            for(size_t i = 0; i < n; i++){
               grid.for_each_candidate(x[i], y[i], [&](int j){
                  if(size_t(j) <= i) return; // count each pair once
                  double dx = x[j] - x[i];
                  double dy = y[j] - y[i];
                  if(periodic){
                     dx = min_image(dx, Lx, true);
                     dy = min_image(dy, Ly, true);
                  }
                  const double min_dist = radius[i] + radius[j];
                  const double dist = std::sqrt(dx*dx + dy*dy);
                  if(dist < min_dist){
                     const double frac = 1.0 - dist/min_dist;
                     if(frac > worst) worst = frac;
                  }
               });
            }

            return worst;

         }

      } // end of anonymous namespace

      //--------------------------------------------------------------------
      grain_packing_result_t pack_grain_seeds(const std::vector<double>& diameters,
                                               double domain_x, double domain_y,
                                               bool periodic, MTRand& rng){

         grain_packing_result_t result;
         const size_t n = diameters.size();

         if(n == 0){
            result.achieved_max_overlap = 0.0;
            result.iterations = 0;
            result.converged = true;
            return result;
         }

         std::vector<double> radius(n);
         double mean_radius = 0.0;
         for(size_t i = 0; i < n; i++){
            radius[i] = 0.5 * diameters[i];
            mean_radius += radius[i];
         }
         mean_radius /= double(n);

         // internal tuning constants for the coarse-then-fine relaxation schedule
         const double inflate_factor = 2.2;
         const int coarse_steps = 400;
         const int fine_steps = 400;
         const double dt_coarse = 0.1 * mean_radius;
         const double dt_fine = 0.05 * mean_radius;
         const double force_const = 0.1;
         const double convergence_tol = 0.02;
         const double max_step = 0.5 * mean_radius;

         double Lx = inflate_factor * domain_x;
         double Ly = inflate_factor * domain_y;

         std::vector<double> x(n), y(n);
         for(size_t i = 0; i < n; i++){
            x[i] = rng() * Lx;
            y[i] = rng() * Ly;
         }

         std::vector<double> fx(n), fy(n);
         int iterations = 0;

         // coarse phase: relax while compressing the box from (Lx,Ly) down
         // to exactly (domain_x,domain_y)
         for(int step = 0; step < coarse_steps; step++){
            calculate_forces(x, y, radius, fx, fy, 0.0, 0.0, Lx, Ly, periodic, force_const);
            move_particles(x, y, fx, fy, dt_coarse, max_step);
            if(periodic) wrap_positions(x, y, Lx, Ly);
            shrink(x, y, Lx, Ly, domain_x, domain_y, coarse_steps - step);
            iterations++;
         }
         Lx = domain_x;
         Ly = domain_y;
         if(periodic) wrap_positions(x, y, Lx, Ly);

         // fine phase: relax at the fixed final box size until converged or budget exhausted
         bool converged = false;
         for(int step = 0; step < fine_steps; step++){
            calculate_forces(x, y, radius, fx, fy, 0.0, 0.0, Lx, Ly, periodic, force_const);
            move_particles(x, y, fx, fy, dt_fine, max_step);
            if(periodic) wrap_positions(x, y, Lx, Ly);
            iterations++;
            if(max_overlap_fraction(x, y, radius, 0.0, 0.0, Lx, Ly, periodic) <= convergence_tol){
               converged = true;
               break;
            }
         }

         if(!periodic){
            for(size_t i = 0; i < n; i++){
               x[i] = std::min(std::max(x[i], 0.0), domain_x);
               y[i] = std::min(std::max(y[i], 0.0), domain_y);
            }
         }

         result.sites.resize(n);
         for(size_t i = 0; i < n; i++) result.sites[i] = point2_t(x[i], y[i]);

         result.achieved_max_overlap = max_overlap_fraction(x, y, radius, 0.0, 0.0, Lx, Ly, periodic);
         result.iterations = iterations;
         result.converged = converged;

         return result;

      }

   } // end of internal namespace
} // end of grains namespace
