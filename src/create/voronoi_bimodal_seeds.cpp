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
#include <vector>

// Vampire headers
#include "random.hpp"

// create module headers
#include "internal.hpp"

namespace create{
namespace internal{

namespace{

   // background grid used for O(1) amortised neighbour queries during
   // random sequential addition (RSA) disc packing
   class neighbour_grid_t{

   public:

      neighbour_grid_t(const double xmin, const double ymin,
                        const double xmax, const double ymax,
                        const double cell_size):
         xmin_(xmin), ymin_(ymin), cell_size_(cell_size)
      {
         nx_ = std::max(1, int(std::ceil((xmax - xmin) / cell_size)));
         ny_ = std::max(1, int(std::ceil((ymax - ymin) / cell_size)));
         cells_.resize(nx_ * ny_);
      }

      // insert seed index i with coordinates (x,y) and radius r
      void insert(const double x, const double y, const int index){
         cells_[cell_index(x, y)].push_back(index);
      }

      // returns true if a disc at (x,y) with radius r overlaps any
      // already-inserted seed, i.e. dist(seed) < r + r_seed
      bool collides(const double x, const double y, const double r,
                    const std::vector<double>& xs, const std::vector<double>& ys,
                    const std::vector<double>& rs) const {

         const int cx = cell_x(x);
         const int cy = cell_y(y);

         for(int j = cy - 1; j <= cy + 1; j++){
            if(j < 0 || j >= ny_) continue;
            for(int i = cx - 1; i <= cx + 1; i++){
               if(i < 0 || i >= nx_) continue;
               const std::vector<int>& cell = cells_[i + j * nx_];
               for(size_t k = 0; k < cell.size(); k++){
                  const int idx = cell[k];
                  const double dx = xs[idx] - x;
                  const double dy = ys[idx] - y;
                  const double min_dist = r + rs[idx];
                  if(dx*dx + dy*dy < min_dist*min_dist) return true;
               }
            }
         }
         return false;
      }

   private:

      int cell_x(const double x) const {
         int i = int((x - xmin_) / cell_size_);
         if(i < 0) i = 0;
         if(i >= nx_) i = nx_ - 1;
         return i;
      }

      int cell_y(const double y) const {
         int j = int((y - ymin_) / cell_size_);
         if(j < 0) j = 0;
         if(j >= ny_) j = ny_ - 1;
         return j;
      }

      int cell_index(const double x, const double y) const {
         return cell_x(x) + cell_y(y) * nx_;
      }

      double xmin_, ymin_, cell_size_;
      int nx_, ny_;
      std::vector< std::vector<int> > cells_;

   };

} // end of anonymous namespace

//-----------------------------------------------------------------------------
// Generates seed points for a bimodal (two grain-size population) Voronoi
// film using random sequential addition (RSA) disc packing on a background
// grid for fast neighbour queries. Seeds are appended to grain_coord_array
// (which may already contain entries) and a parallel is_small_grain flag is
// appended to record which population each accepted seed belongs to. Returns
// the number of seeds placed.
//-----------------------------------------------------------------------------
int generate_bimodal_voronoi_seeds(std::vector <std::vector <double> >& grain_coord_array,
                                    std::vector<bool>& is_small_grain,
                                    const double domain_x, const double domain_y,
                                    const double std_diameter,   const double std_variance,
                                    const double small_diameter, const double small_variance,
                                    const double small_fraction){

   const double r_std   = 0.5 * std_diameter;
   const double r_small = 0.5 * small_diameter;
   const double r_max   = std::max(r_std, r_small);

   const double pad = 2.0 * r_max;

   const double xmin = -pad;
   const double xmax = domain_x + pad;
   const double ymin = -pad;
   const double ymax = domain_y + pad;

   const double cell_size = r_std + r_small;

   neighbour_grid_t grid(xmin, ymin, xmax, ymax, cell_size);

   std::vector<double> xs;
   std::vector<double> ys;
   std::vector<double> rs;

   const int FAILURE_LIMIT = 5000;
   int consecutive_failures = 0;
   int placed = 0;
   int placed_small = 0;

   while(consecutive_failures < FAILURE_LIMIT){

      // choose which population to attempt next based on which is currently
      // furthest below its target share, rather than an independent coin
      // flip per attempt. A per-attempt Bernoulli draw tracks the requested
      // fraction poorly here: small grains need much less clearance than
      // standard ones, so they keep succeeding long after standard grains
      // have effectively jammed out, and the realised fraction drifts far
      // above the request (observed ~0.3 requested -> ~0.85 realised on a
      // typical 5nm/2nm film). Always attempting the currently-deficient
      // population keeps the realised fraction close to the request, at the
      // cost of placing fewer total grains once the larger population jams.
      const bool is_small = (double(placed_small) < small_fraction * double(placed + 1));
      const double mean_d = is_small ? small_diameter : std_diameter;
      const double sigma  = is_small ? small_variance : std_variance;

      double r_i = 0.5 * std::exp(std::log(mean_d) + sigma * mtrandom::gaussian());
      // guard against pathological lognormal tail draws
      const double r_min_clamp = 0.2 * mean_d * 0.5;
      const double r_max_clamp = 2.0 * mean_d * 0.5;
      if(r_i < r_min_clamp) r_i = r_min_clamp;
      if(r_i > r_max_clamp) r_i = r_max_clamp;

      const double x = xmin + mtrandom::grnd() * (xmax - xmin);
      const double y = ymin + mtrandom::grnd() * (ymax - ymin);

      if(!grid.collides(x, y, r_i, xs, ys, rs)){
         xs.push_back(x);
         ys.push_back(y);
         rs.push_back(r_i);
         grid.insert(x, y, int(xs.size()) - 1);

         grain_coord_array.push_back(std::vector<double>());
         grain_coord_array.back().push_back(x);
         grain_coord_array.back().push_back(y);
         is_small_grain.push_back(is_small);

         placed++;
         if(is_small) placed_small++;
         consecutive_failures = 0;
      }
      else{
         consecutive_failures++;
      }
   }

   return placed;

}

} // end of internal namespace
} // end of create namespace
