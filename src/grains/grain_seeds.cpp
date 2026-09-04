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
#include <fstream>

// Vampire headers
#include "errors.hpp"
#include "random.hpp"
#include "vio.hpp"

// grains module headers
#include "internal.hpp"

namespace grains{
   namespace internal{

      namespace{

         // Mean disc area pi/4*E[d^2] for a diameter d drawn from a lognormal
         // distribution parameterised so its arithmetic mean is mean_d.
         inline double lognormal_mean_area(const double mean_d, const double sigma){
            return (M_PI/4.0) * mean_d*mean_d * std::exp(sigma*sigma);
         }

         //--------------------------------------------------------------------
         // helpers for grain-size statistics reporting
         //--------------------------------------------------------------------
         polygon_t to_polygon(const std::vector<std::vector<double> >& verts){
            polygon_t p;
            p.reserve(verts.size());
            for(size_t i = 0; i < verts.size(); i++) p.push_back(point2_t(verts[i][0], verts[i][1]));
            return p;
         }

         // true if an edge's two endpoints both lie on the same domain box side
         bool on_same_box_side(const point2_t& a, const point2_t& b,
                               double xmin, double ymin, double xmax, double ymax, double eps){
            if(std::fabs(a.x-xmin) < eps && std::fabs(b.x-xmin) < eps) return true;
            if(std::fabs(a.x-xmax) < eps && std::fabs(b.x-xmax) < eps) return true;
            if(std::fabs(a.y-ymin) < eps && std::fabs(b.y-ymin) < eps) return true;
            if(std::fabs(a.y-ymax) < eps && std::fabs(b.y-ymax) < eps) return true;
            return false;
         }

      } // end of anonymous namespace

      //--------------------------------------------------------------------
      std::vector<double> read_grain_size_distribution_file(const std::string& filename){

         std::vector<double> diameters;

         std::ifstream ifile(filename.c_str());
         if(!ifile.is_open()){
            terminaltextcolor(RED);
            std::cerr << "Error: could not open create:grain-size-distribution-file \"" << filename << "\"" << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error: could not open create:grain-size-distribution-file \"" << filename << "\"" << std::endl;
            err::vexit();
         }

         double d;
         while(ifile >> d){
            if(d > 0.0) diameters.push_back(d);
         }

         if(diameters.empty()){
            terminaltextcolor(RED);
            std::cerr << "Error: create:grain-size-distribution-file \"" << filename
                      << "\" contains no valid (positive) diameters" << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error: create:grain-size-distribution-file \"" << filename
                 << "\" contains no valid (positive) diameters" << std::endl;
            err::vexit();
         }

         return diameters;

      }

      //--------------------------------------------------------------------
      double mean_grain_area(const grain_size_distribution_params_t& params){

         switch(params.distribution){

            case grain_size_delta:
               return (M_PI/4.0) * params.mean_diameter * params.mean_diameter;

            case grain_size_normal: {
               const double sigma = params.sd * params.mean_diameter;
               return (M_PI/4.0) * (params.mean_diameter*params.mean_diameter + sigma*sigma);
            }

            case grain_size_lognormal:
               return lognormal_mean_area(params.mean_diameter, params.sd);

            case grain_size_bimodal: {
               const double a1 = lognormal_mean_area(params.mean_diameter, params.sd);
               const double a2 = lognormal_mean_area(params.second_mean_diameter, params.second_sd);
               return (1.0 - params.second_fraction)*a1 + params.second_fraction*a2;
            }

            case grain_size_file:
            default: {
               if(params.file_diameters.empty()) return 1.0;
               double sum = 0.0;
               for(size_t i = 0; i < params.file_diameters.size(); i++){
                  const double d = params.file_diameters[i];
                  sum += (M_PI/4.0) * d * d;
               }
               return sum / double(params.file_diameters.size());
            }

         }

      }

      //--------------------------------------------------------------------
      double bimodal_area_fraction(const grain_size_distribution_params_t& params){

         const double a1 = lognormal_mean_area(params.mean_diameter, params.sd);
         const double a2 = lognormal_mean_area(params.second_mean_diameter, params.second_sd);
         const double f2 = params.second_fraction;
         const double total = (1.0 - f2)*a1 + f2*a2;

         return (total > 0.0) ? (f2*a2/total) : 0.0;

      }

      //--------------------------------------------------------------------
      int choose_grain_count(double domain_x, double domain_y, const grain_size_distribution_params_t& params){

         if(params.distribution == grain_size_file) return int(params.file_diameters.size());

         const double domain_area = domain_x * domain_y;
         const double area = mean_grain_area(params);

         int n = int(std::lround(domain_area / std::max(area, 1.0e-12)));
         if(n < 1) n = 1;

         return n;

      }

      //--------------------------------------------------------------------
      std::vector<grain_size_sample_t> sample_grain_diameters(int n, const grain_size_distribution_params_t& params, MTRand& rng){

         std::vector<grain_size_sample_t> out;

         // file mode: the file IS the target list, exactly as given
         if(params.distribution == grain_size_file){
            out.reserve(params.file_diameters.size());
            for(size_t i = 0; i < params.file_diameters.size(); i++){
               grain_size_sample_t s;
               s.diameter = params.file_diameters[i];
               s.population = 0;
               out.push_back(s);
            }
            return out;
         }

         out.reserve(std::max(n, 0));

         for(int i = 0; i < n; i++){

            int population = 0;
            double mean_d = params.mean_diameter;
            double sigma  = params.sd;

            if(params.distribution == grain_size_bimodal && rng() < params.second_fraction){
               population = 1;
               mean_d = params.second_mean_diameter;
               sigma  = params.second_sd;
            }

            double d;

            if(params.distribution == grain_size_delta){
               d = params.mean_diameter;
            }
            else if(params.distribution == grain_size_normal){
               const double abs_sigma = params.sd * params.mean_diameter;
               d = params.mean_diameter + abs_sigma * mtrandom::gaussianc(rng);
               // guard against a pathological negative/near-zero draw
               if(d < 0.05*params.mean_diameter) d = 0.05*params.mean_diameter;
            }
            else{ // lognormal (also used per-population inside bimodal)
               // mu_log chosen so the sampled mean diameter equals mean_d
               const double mu_log = std::log(mean_d) - 0.5*sigma*sigma;
               d = std::exp(mu_log + sigma * mtrandom::gaussianc(rng));
            }

            grain_size_sample_t s;
            s.diameter = d;
            s.population = population;
            out.push_back(s);

         }

         return out;

      }

      //--------------------------------------------------------------------
      std::vector<point2_t> lloyd_polish_seeds(const std::vector<point2_t>& sites_in,
                                                double xmin, double ymin, double xmax, double ymax,
                                                double regularity, bool periodic){

         std::vector<point2_t> sites = sites_in;
         const size_t n = sites.size();
         if(n == 0 || regularity <= 0.0) return sites;

         // relax toward a centroidal Voronoi tessellation (Lloyd, 1982)
         const int MAX_SWEEPS = 4;
         const int sweeps = std::max(1, int(std::lround(regularity * MAX_SWEEPS)));
         const std::vector<double> uniform_weights(n, 1.0);
         const double Lx = xmax - xmin;
         const double Ly = ymax - ymin;

         for(int s = 0; s < sweeps; s++){

            const std::vector<polygon_t> cells = build_power_cells(sites, uniform_weights, xmin, ymin, xmax, ymax, periodic);

            for(size_t i = 0; i < n; i++){
               if(cells[i].size() < 3) continue;
               const point2_t c = polygon_centroid(cells[i]);
               sites[i].x += regularity * (c.x - sites[i].x);
               sites[i].y += regularity * (c.y - sites[i].y);
            }

            // keep positions in the canonical [xmin,xmax) x [ymin,ymax) range
            if(periodic){
               for(size_t i = 0; i < n; i++){
                  sites[i].x = std::fmod(sites[i].x - xmin, Lx); if(sites[i].x < 0.0) sites[i].x += Lx; sites[i].x += xmin;
                  sites[i].y = std::fmod(sites[i].y - ymin, Ly); if(sites[i].y < 0.0) sites[i].y += Ly; sites[i].y += ymin;
               }
            }

         }

         return sites;

      }

      //--------------------------------------------------------------------
      std::vector<double> generate_grain_spacing_jitter(int num_grains, double base_delta,
                                                          double relative_sd, MTRand& rng){

         std::vector<double> delta(std::max(num_grains, 0), base_delta);
         if(relative_sd <= 0.0) return delta;

         for(size_t i = 0; i < delta.size(); i++){
            const double d = base_delta * (1.0 + relative_sd * mtrandom::gaussianc(rng));
            delta[i] = std::max(0.0, d);
         }

         return delta;

      }

      //--------------------------------------------------------------------
      weight_fit_result_t fit_grain_weights(const std::vector<point2_t>& sites,
                                             const std::vector<double>& target_diameters,
                                             double xmin, double ymin, double xmax, double ymax,
                                             double tol, int max_iter, double kappa, bool periodic){

         const size_t n = sites.size();

         weight_fit_result_t result;
         result.iterations = 0;
         result.achieved_tol = 0.0;
         result.converged = (n == 0);
         result.num_vanished = 0;

         std::vector<double> target_area(n);
         std::vector<double> weights(n);
         double min_weight_scale = 1.0e300;
         for(size_t i = 0; i < n; i++){
            const double r = 0.5 * target_diameters[i];
            weights[i] = r*r;
            target_area[i] = M_PI * r * r;
            min_weight_scale = std::min(min_weight_scale, target_area[i]/M_PI);
         }
         if(n == 0) min_weight_scale = 1.0;

         for(int iter = 0; iter < max_iter && n > 0; iter++){

            result.iterations = iter + 1;

            const std::vector<polygon_t> cells = build_power_cells(sites, weights, xmin, ymin, xmax, ymax, periodic);

            std::vector<double> area(n, 0.0);
            double max_rel_err = 0.0;
            for(size_t i = 0; i < n; i++){
               area[i] = (cells[i].size() >= 3) ? polygon_area(cells[i]) : 0.0;
               const double rel_err = std::fabs(area[i] - target_area[i]) / target_area[i];
               if(rel_err > max_rel_err) max_rel_err = rel_err;
            }

            result.achieved_tol = max_rel_err;
            if(max_rel_err < tol){ result.converged = true; break; }

            for(size_t i = 0; i < n; i++){
               const double delta_w = kappa * (target_area[i] - area[i]) / M_PI;
               double new_w = weights[i] + delta_w;
               // step capped by the smallest target area in the seed set
               const double max_step = min_weight_scale;
               if(new_w > weights[i] + max_step) new_w = weights[i] + max_step;
               if(new_w < weights[i] - max_step) new_w = weights[i] - max_step;
               weights[i] = new_w;
            }

         }

         result.weights = weights;

         if(n > 0){
            const std::vector<polygon_t> final_cells = build_power_cells(sites, weights, xmin, ymin, xmax, ymax, periodic);
            int vanished = 0;
            for(size_t i = 0; i < n; i++) if(final_cells[i].size() < 3) vanished++;
            result.num_vanished = vanished;
         }

         return result;

      }

      //--------------------------------------------------------------------
      grain_seed_set_t generate_grain_seeds(double domain_x, double domain_y,
                                             const grain_size_distribution_params_t& params,
                                             MTRand& rng){

         grain_seed_set_t out;

         const int n = choose_grain_count(domain_x, domain_y, params);
         const std::vector<grain_size_sample_t> samples = sample_grain_diameters(n, params, rng);

         const size_t num_seeds = samples.size();
         std::vector<double> diameters(num_seeds);
         out.population.resize(num_seeds);
         out.target_diameter.resize(num_seeds);
         for(size_t i = 0; i < num_seeds; i++){
            diameters[i] = samples[i].diameter;
            out.population[i] = samples[i].population;
            out.target_diameter[i] = samples[i].diameter;
         }

         // seed placement always runs periodic internally, regardless of the
         // requested final output periodicity (params.periodic below)
         const grain_packing_result_t packing = pack_grain_seeds(diameters, domain_x, domain_y, true, rng);
         out.sites = lloyd_polish_seeds(packing.sites, 0.0, 0.0, domain_x, domain_y, params.regularity, true);
         out.pack = packing;
         out.fit = fit_grain_weights(out.sites, diameters, 0.0, 0.0, domain_x, domain_y, 0.03, 500, 0.5, params.periodic);
         out.weights = out.fit.weights;

         return out;

      }

      //--------------------------------------------------------------------
      grain_statistics_summary_t compute_grain_size_statistics(
         const std::vector<std::vector<double> >& grain_coord_array,
         const std::vector<std::vector<std::vector<double> > >& grain_vertices_array,
         const std::vector<int>& population){

         grain_statistics_summary_t s;
         s.placed_primary = 0;
         s.placed_second = 0;
         s.num_vanished = 0;

         double sum_d1 = 0.0, sum_d1_sq = 0.0;
         double sum_d2 = 0.0, sum_d2_sq = 0.0;
         double area_primary = 0.0, area_second = 0.0;

         for(size_t g = 0; g < grain_coord_array.size(); g++){

            const size_t nv = grain_vertices_array[g].size();
            if(nv < 3){ s.num_vanished++; continue; }

            const polygon_t poly = to_polygon(grain_vertices_array[g]);
            const double area = std::fabs(polygon_area(poly));
            const double diameter = 2.0 * std::sqrt(area/M_PI);

            const bool is_second = (g < population.size()) ? (population[g] == 1) : false;

            if(is_second){
               s.placed_second++;
               sum_d2 += diameter;
               sum_d2_sq += diameter*diameter;
               area_second += area;
            }
            else{
               s.placed_primary++;
               sum_d1 += diameter;
               sum_d1_sq += diameter*diameter;
               area_primary += area;
            }
         }

         s.mean_diameter_primary = (s.placed_primary > 0) ? sum_d1/s.placed_primary : 0.0;
         s.mean_diameter_second  = (s.placed_second  > 0) ? sum_d2/s.placed_second  : 0.0;

         s.sd_diameter_primary = (s.placed_primary > 1) ?
            std::sqrt(std::max(0.0, sum_d1_sq/s.placed_primary - s.mean_diameter_primary*s.mean_diameter_primary)) : 0.0;
         s.sd_diameter_second = (s.placed_second > 1) ?
            std::sqrt(std::max(0.0, sum_d2_sq/s.placed_second - s.mean_diameter_second*s.mean_diameter_second)) : 0.0;

         const int placed = s.placed_primary + s.placed_second;
         s.realised_number_fraction = (placed > 0) ? double(s.placed_second)/double(placed) : 0.0;

         const double total_area = area_primary + area_second;
         s.realised_area_fraction = (total_area > 0.0) ? area_second/total_area : 0.0;

         return s;

      }

      //--------------------------------------------------------------------
      void write_grain_statistics_file(const std::string& filename,
         const std::vector<std::vector<double> >& grain_coord_array,
         const std::vector<std::vector<std::vector<double> > >& grain_vertices_array,
         const std::vector<int>& population,
         const std::vector<double>& target_diameter,
         double xmin, double ymin, double xmax, double ymax){

         std::ofstream ofile(filename.c_str());
         if(!ofile.is_open()) return;

         ofile << "# id\tx\ty\tarea\tdiameter\tvertices\tneighbours\tpopulation\ttarget_diameter" << std::endl;

         const double eps = 1.0e-6 * std::max(xmax-xmin, ymax-ymin);

         for(size_t g = 0; g < grain_coord_array.size(); g++){

            const size_t nv = grain_vertices_array[g].size();
            if(nv < 3) continue;

            const polygon_t poly = to_polygon(grain_vertices_array[g]);
            const double area = std::fabs(polygon_area(poly));
            const double diameter = 2.0 * std::sqrt(area/M_PI);

            int boundary_edges = 0;
            for(size_t v = 0; v < nv; v++){
               const point2_t& a = poly[v];
               const point2_t& b = poly[(v+1)%nv];
               if(on_same_box_side(a, b, xmin, ymin, xmax, ymax, eps)) boundary_edges++;
            }
            const int neighbours = int(nv) - boundary_edges;

            ofile << g << "\t" << grain_coord_array[g][0] << "\t" << grain_coord_array[g][1] << "\t"
                  << area << "\t" << diameter << "\t" << nv << "\t" << neighbours << "\t"
                  << (g < population.size() ? population[g] : 0) << "\t"
                  << (g < target_diameter.size() ? target_diameter[g] : 0.0) << std::endl;

         }

      }

   } // end of internal namespace
} // end of grains namespace
