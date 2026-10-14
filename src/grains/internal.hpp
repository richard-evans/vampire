#ifndef GRAINS_INTERNAL_H_
#define GRAINS_INTERNAL_H_
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

//---------------------------------------------------------------------
// Internal data structures and functions for the grains module. Not
// to be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <fstream>
#include <functional>
#include <list>
#include <string>
#include <vector>

// Vampire headers
#include "grains.hpp"
#include "create.hpp"
#include "mtrand.hpp"

namespace grains{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Internal data type definitions
      //-----------------------------------------------------------------------------

      // create:grain-tessellation: plain Voronoi or weighted Laguerre (power) diagram
      enum grain_tessellation_t { tessellation_voronoi, tessellation_laguerre };

      // create:grain-size-distribution
      enum grain_size_distribution_t { grain_size_delta, grain_size_normal, grain_size_lognormal, grain_size_bimodal, grain_size_file };

      // create:grain-facet-orientation
      enum grain_facet_orientation_t { facet_orientation_random, facet_orientation_fixed, facet_orientation_textured };

      //-----------------------------------------------------------------------------
      // Per-grain magnetic-property bookkeeping, populated by set_properties()
      // (grains::) once atom creation has assigned atoms::grain_array.
      //-----------------------------------------------------------------------------

      extern std::vector<int> grain_size_array;   // atoms per grain

      extern std::vector<double> x_coord_array;   // mean grain coordinates
      extern std::vector<double> y_coord_array;
      extern std::vector<double> z_coord_array;

      extern std::vector<double> sat_mag_array;   // summed saturation magnetic moment per grain

      extern std::vector<double> orientation_array; // per-grain Wulff facet orientation, radians

      // Renumbers grains so no empty (atom-less) grain id remains.
      void remix_grain_numbers();

      // Registers the raw (pre-remix) per-grain orientation array.
      void set_raw_orientation(const std::vector<double>& raw_orientation);

      //-----------------------------------------------------------------------------
      // Convex polygon primitives (CCW vertex lists) used to build and clip
      // grain cells. clip_halfplane() is the fundamental primitive used
      // throughout tessellation, spacing, rounding and faceting below.
      //-----------------------------------------------------------------------------

      // tolerance for on/inside/outside plane classification
      constexpr double grain_geometry_epsilon = 1.0e-10;

      // simple 2D point
      struct point2_t{
         double x;
         double y;
         point2_t(double xi = 0.0, double yi = 0.0) : x(xi), y(yi) {}
      };

      // CCW-ordered convex polygon
      typedef std::vector<point2_t> polygon_t;

      // shoelace formula; signed, positive for CCW input
      double polygon_area(const polygon_t& poly);

      // area-weighted centroid
      point2_t polygon_centroid(const polygon_t& poly);

      // maximum distance from origin to any vertex of poly
      double max_vertex_radius(const polygon_t& poly, point2_t origin);

      // Sutherland-Hodgman clip of poly against the half-plane nx*x+ny*y <= c
      polygon_t clip_halfplane(const polygon_t& poly, double nx, double ny, double c);

      // fixes vertex winding to CCW if needed
      polygon_t ensure_ccw(const polygon_t& poly);

      // shrinks poly by moving every edge inward by delta
      polygon_t offset_inward(const polygon_t& poly, double delta);

      // axis-aligned rectangle [xmin,xmax] x [ymin,ymax], CCW
      polygon_t rectangle(double xmin, double ymin, double xmax, double ymax);

      // facet count used to approximate a circle when rounding a grain outline
      constexpr int grain_rounding_num_facets = 64;

      // Rounds poly's corners by clipping to a bisected radius so the
      // retained area is exactly area_fraction * polygon_area(poly).
      polygon_t round_polygon(const polygon_t& poly, point2_t centre, double area_fraction,
                               int num_facets = grain_rounding_num_facets);

      // Organic-rounding overload: as above but clips toward a random
      // low-order cosine outline instead of a common circular radius.
      polygon_t round_polygon(const polygon_t& poly, point2_t centre, double area_fraction,
                               double roughness_amplitude, int roughness_modes, MTRand& rng,
                               int num_facets = grain_rounding_num_facets);

      //-----------------------------------------------------------------------------
      // Physically-relaxed disc-packing seed placement for a granular film.
      //-----------------------------------------------------------------------------

      struct grain_packing_result_t{
         std::vector<point2_t> sites;
         double achieved_max_overlap; // largest remaining overlap fraction between any pair
         int iterations;
         bool converged;
      };

      // Places one seed per entry of diameters by relaxing a random disc
      // packing down to the [0,domain_x] x [0,domain_y] rectangle.
      grain_packing_result_t pack_grain_seeds(const std::vector<double>& diameters,
                                               double domain_x, double domain_y,
                                               bool periodic, MTRand& rng);

      //-----------------------------------------------------------------------------
      // Laguerre (power) diagram construction. Each site i carries a weight
      // w_i; a point belongs to cell i if |x-p_i|^2 - w_i is smallest. Equal
      // weights reduce this to an ordinary Voronoi diagram.
      //-----------------------------------------------------------------------------

      // Builds the power cell of every site, clipped to [xmin,xmax]x[ymin,ymax].
      // periodic builds against the minimum-image tiling instead.
      std::vector<polygon_t> build_power_cells(
         const std::vector<point2_t>& sites,
         const std::vector<double>&   weights,
         double xmin, double ymin, double xmax, double ymax,
         bool periodic);

      // Removes every grain whose cell touches the domain boundary (non-
      // periodic only), leaving a void margin at the film edges. Returns
      // the number of grains removed.
      int remove_boundary_grains(std::vector<std::vector<double> >& grain_coord_array,
                                  std::vector<std::vector<std::vector<double> > >& grain_vertices_array,
                                  double xmin, double ymin, double xmax, double ymax);

      //-----------------------------------------------------------------------------
      // Wulff faceting: clips a tessellated grain cell toward the
      // equilibrium crystallite shape (Wulff construction), blended by a
      // strength parameter between the unclipped cell and the pure facet shape.
      //-----------------------------------------------------------------------------

      // Clips poly by num_facets half-planes at angles theta0 + 2*pi*k/num_facets,
      // blended between poly's own outline (strength=0) and a pure Wulff
      // polygon of area target_area (strength=1). anisotropy scales
      // alternating facet families relative to each other.
      polygon_t apply_wulff_facets(const polygon_t& poly, point2_t centre,
                                    int num_facets, double theta0,
                                    double anisotropy, double strength,
                                    double target_area);

      // Array-format adapter over apply_wulff_facets(); must run before the
      // boundary-spacing offset since grain_vertices_array is still absolute here.
      void apply_grain_faceting(std::vector<std::vector<double> >& grain_coord_array,
                                 std::vector<std::vector<std::vector<double> > >& grain_vertices_array,
                                 const std::vector<double>& orientation,
                                 int num_facets, double anisotropy, double strength);

      // Per-grain facet orientation theta_i, selected by
      // create:grain-facet-orientation (random/fixed/textured).
      std::vector<double> generate_grain_facet_orientations(int num_grains, int num_facets,
                                                              grain_facet_orientation_t mode,
                                                              double angle_degrees, double spread_degrees,
                                                              MTRand& rng);

      //-----------------------------------------------------------------------------
      // Grain-size distribution sampling, seed placement and weight fitting.
      //-----------------------------------------------------------------------------

      // create:grain-size-* input parameters needed to sample a target diameter.
      struct grain_size_distribution_params_t{
         grain_size_distribution_t distribution;
         double mean_diameter;               // create:grain-size
         double sd;                          // create:grain-size-sd (lognormal log-sd; normal sd as a fraction of the mean)
         double second_mean_diameter;        // create:grain-size-second
         double second_sd;                   // create:grain-size-second-sd
         double second_fraction;             // create:grain-size-second-fraction, BY NUMBER
         std::vector<double> file_diameters; // create:grain-size-distribution-file, pre-read
         double regularity;                  // create:grain-regularity
         bool periodic;                      // create:grain-periodic-boundaries

         grain_size_distribution_params_t():
            distribution(grain_size_normal), mean_diameter(100.0), sd(0.3),
            second_mean_diameter(20.0), second_sd(0.15), second_fraction(0.3),
            regularity(0.0), periodic(false)
         {}
      };

      // one sampled diameter, tagged with population (0 = primary, 1 = second/bimodal)
      struct grain_size_sample_t{
         double diameter;
         int population;
      };

      // Reads a whitespace/newline-separated list of positive diameters (Angstroms).
      std::vector<double> read_grain_size_distribution_file(const std::string& filename);

      // mean disc area of the requested distribution, in closed form
      double mean_grain_area(const grain_size_distribution_params_t& params);

      // bimodal number-fraction -> area-fraction conversion
      double bimodal_area_fraction(const grain_size_distribution_params_t& params);

      // number of seeds needed so the tessellation realises the requested size distribution
      int choose_grain_count(double domain_x, double domain_y, const grain_size_distribution_params_t& params);

      // draws n target diameters from the requested distribution
      std::vector<grain_size_sample_t> sample_grain_diameters(int n, const grain_size_distribution_params_t& params, MTRand& rng);

      // Relaxes sites toward their own Voronoi cell centroids (Lloyd, 1982),
      // polishing pack_grain_seeds()'s output toward a more uniform pattern.
      std::vector<point2_t> lloyd_polish_seeds(const std::vector<point2_t>& sites,
                                                double xmin, double ymin, double xmax, double ymax,
                                                double regularity, bool periodic);

      // Per-grain jittered inward-offset distances for create:grain-spacing-sd.
      std::vector<double> generate_grain_spacing_jitter(int num_grains, double base_delta,
                                                          double relative_sd, MTRand& rng);

      struct weight_fit_result_t{
         std::vector<double> weights;
         double achieved_tol;
         int iterations;
         bool converged;
         int num_vanished; // sites whose final power cell has < 3 vertices
      };

      // Fits each site's weight so its power-cell area matches its target
      // grain area, by damped fixed-point iteration.
      weight_fit_result_t fit_grain_weights(const std::vector<point2_t>& sites,
                                             const std::vector<double>& target_diameters,
                                             double xmin, double ymin, double xmax, double ymax,
                                             double tol = 0.03, int max_iter = 500, double kappa = 0.5,
                                             bool periodic = false);

      struct grain_seed_set_t{
         std::vector<point2_t> sites;
         std::vector<double> weights;
         std::vector<double> target_diameter;
         std::vector<int> population;
         grain_packing_result_t pack;
         weight_fit_result_t fit;
      };

      // Top-level driver: sample diameters -> choose N -> pack seeds -> Lloyd polish -> fit weights.
      grain_seed_set_t generate_grain_seeds(double domain_x, double domain_y,
                                             const grain_size_distribution_params_t& params,
                                             MTRand& rng);

      //-----------------------------------------------------------------------
      // Grain-size statistics reporting, run on the freshly tessellated
      // cells before spacing/rounding/faceting shrink them.
      //-----------------------------------------------------------------------
      struct grain_statistics_summary_t{
         int placed_primary, placed_second;
         double mean_diameter_primary, sd_diameter_primary;
         double mean_diameter_second, sd_diameter_second;
         double realised_number_fraction; // second population, by number
         double realised_area_fraction;   // second population, by area
         int num_vanished;
      };

      grain_statistics_summary_t compute_grain_size_statistics(
         const std::vector<std::vector<double> >& grain_coord_array,
         const std::vector<std::vector<std::vector<double> > >& grain_vertices_array,
         const std::vector<int>& population);

      // One line per grain: id, centre, area, diameter, vertex/neighbour count, population, target diameter.
      void write_grain_statistics_file(const std::string& filename,
         const std::vector<std::vector<double> >& grain_coord_array,
         const std::vector<std::vector<std::vector<double> > >& grain_vertices_array,
         const std::vector<int>& population,
         const std::vector<double>& target_diameter,
         double xmin, double ymin, double xmax, double ymax);

      //-----------------------------------------------------------------------------
      // Scaffolding shared by every driver that assigns atoms to grains:
      // bins atoms into a supercell grid so each grain only scans its own
      // bounding box rather than every atom in the system.
      //-----------------------------------------------------------------------------

      // 2D binning of atom indices by (unit-cell-x, unit-cell-y) supercell
      typedef std::vector<std::vector<std::vector<int> > > supercell_array_t;

      supercell_array_t build_supercell_array(const std::vector<cs::catom_t>& catom_array,
                                               double unit_cell_dim_x, double unit_cell_dim_y,
                                               int num_unit_cells_x, int num_unit_cells_y);

      // Bounding supercell range of one grain's polygon (relative to its
      // centre x0,y0), plus vertices flattened for vmath::point_in_polygon*().
      struct grain_footprint_t{
         int minx, maxx, miny, maxy;
         std::vector<double> px, py;
         bool periodic = false;
         int num_unit_cells_x = 0;
         int num_unit_cells_y = 0;
      };

      grain_footprint_t compute_grain_footprint(const std::vector<std::vector<double> >& vertices,
                                                 double x0, double y0,
                                                 double unit_cell_dim_x, double unit_cell_dim_y,
                                                 bool periodic = false,
                                                 int num_unit_cells_x = 0,
                                                 int num_unit_cells_y = 0);

      // wraps d to its minimum-image value under a periodic domain of width `domain`
      double wrap_minimum_image(double d, double domain);

      // Visits every atom binned into a supercell within fp's bounding box.
      void assign_atoms_in_footprint(const supercell_array_t& supercell_array,
                                      const grain_footprint_t& fp,
                                      const std::function<void(int atom)>& visit);

      //-----------------------------------------------------------------------------
      // Internal shared variables for the grains module
      //-----------------------------------------------------------------------------

      // own random number generator, reseeded from grain_structure_seed each call
      extern MTRand grnd;

      extern bool generate_voronoi_substructure;
      extern double voronoi_grain_substructure_size;    // mean grain size of the substructure within a particle
      extern double voronoi_grain_substructure_spacing; // spacing between substructure grains
      extern double voronoi_grain_substructure_crystallization_radius;
      extern double voronoi_grain_substructure_overlap_factor;

      // hexagonal-lattice seeding parameters for the grain SUBSTRUCTURE path only
      extern bool parity;
      extern double voronoi_sd;

      extern double voronoi_elliptical_rounding;        // degree of grain cap rounding (0 = vertical/"shear" sides, 1 = fully rounded/ellipsoidal cap)
      extern double voronoi_elliptical_rounding_height; // fraction of grain_film_height, measured from the grain base, at which the flat column ends and the rounded cap begins
      extern double grain_film_height;                  // physical grain/film thickness (Angstroms) that voronoi_elliptical_rounding_height is measured against; < 0 (unset) falls back to cs::system_dimensions[2] for backwards compatibility

      extern bool output_gv_file; // toggle output of grain positions to file

      extern grain_tessellation_t grain_tessellation; // create:grain-tessellation engine selector

      //-----------------------------------------------------------------------------
      // Grain size distribution and weight fitting (create:grain-size-*).
      // Sentinel < 0 means "not set", falling back to dimensions:particle-size/spacing.
      //-----------------------------------------------------------------------------
      extern double grain_size;                    // create:grain-size, mean diameter (Angstroms)
      extern double grain_spacing;                  // create:grain-spacing, boundary gap between grains (Angstroms)
      extern grain_size_distribution_t grain_size_distribution; // create:grain-size-distribution
      extern double grain_size_sd;                  // create:grain-size-sd
      extern double grain_size_second;               // create:grain-size-second
      extern double grain_size_second_sd;            // create:grain-size-second-sd
      extern double grain_size_second_fraction;      // create:grain-size-second-fraction, BY NUMBER

      // true once a non-delta-only create:grain-size-* keyword has been parsed
      extern bool grain_size_shape_param_explicitly_set;
      extern std::string grain_size_distribution_file; // create:grain-size-distribution-file
      extern double grain_regularity;                // create:grain-regularity
      extern int grain_structure_seed;                // create:grain-structure-random-seed
      extern bool output_grain_statistics_file;       // create:grain-statistics-output

      extern double grain_rounding; // create:grain-rounding, area fraction retained (0 = off)
      extern double grain_boundary_roughness; // create:grain-boundary-roughness, organic-rounding amplitude (0 = off)
      extern int grain_boundary_roughness_modes; // create:grain-boundary-roughness-modes
      extern double grain_spacing_sd; // create:grain-spacing-sd, fractional jitter on the boundary gap (0 = off)
      extern bool grain_periodic_boundaries; // create:grain-periodic-boundaries

      //-----------------------------------------------------------------------------
      // Wulff faceting (create:grain-facet-*). grain_facet_symmetry <= 0 is the master off-switch.
      //-----------------------------------------------------------------------------
      extern int grain_facet_symmetry;              // create:grain-facet-symmetry: 0 (off), 4, 6 or 8
      extern double grain_facet_strength;            // create:grain-facet-strength: 0 (untouched) .. 1 (pure Wulff habit)
      extern double grain_facet_anisotropy;          // create:grain-facet-anisotropy: odd/even facet-family distance ratio
      extern grain_facet_orientation_t grain_facet_orientation; // create:grain-facet-orientation
      extern double grain_facet_orientation_angle;   // create:grain-facet-orientation-angle, degrees
      extern double grain_facet_orientation_spread;  // create:grain-facet-orientation-spread, degrees (textured mode only)

      //-----------------------------------------------------------------------------
      // Internal functions for grains module
      //-----------------------------------------------------------------------------

      extern void voronoi_substructure(std::vector<cs::catom_t> & catom_array);

      // rounds every grain's corners to area_fraction of its own area (create:grain-rounding)
      void voronoi_grain_rounding(std::vector <std::vector <double> > & grain_coord_array,
                                  std::vector <std::vector <std::vector <double> > > &  grain_vertices_array,
                                  double area_fraction);

      // organic-rounding overload: each grain gets an independent random blob outline
      void voronoi_grain_rounding(std::vector <std::vector <double> > & grain_coord_array,
                                  std::vector <std::vector <std::vector <double> > > &  grain_vertices_array,
                                  double area_fraction,
                                  double roughness_amplitude, int roughness_modes, MTRand& rng);

      double elliptical_rounding_factor(const double z);

      // Builds the tessellation via build_power_cells() and writes it into
      // the grain_coord_array/grain_vertices_array format used throughout
      // this module. use_laguerre_weights selects weighted vs plain Voronoi.
      void populate_vertex_points_power(std::vector <std::vector <double> > & grain_coord_array,
                                         std::vector <std::vector <std::vector <double> > > &  grain_vertices_array,
                                         const std::vector<double>& weights,
                                         double domain_x, double domain_y,
                                         bool use_laguerre_weights,
                                         bool periodic);

      // inward-offsets every grain cell by delta, giving neighbours a constant boundary gap
      void apply_grain_spacing(std::vector <std::vector <double> > & grain_coord_array,
                                std::vector <std::vector <std::vector <double> > > &  grain_vertices_array,
                                double delta);

      // per-grain overload: delta[i] is grain i's own inward offset (create:grain-spacing-sd)
      void apply_grain_spacing(std::vector <std::vector <double> > & grain_coord_array,
                                std::vector <std::vector <std::vector <double> > > &  grain_vertices_array,
                                const std::vector<double>& delta);

      extern void write_grain_vertices(int id, double dx, double dy, std::ofstream& ofile, std::vector< std::vector <double> >& vertices);

      extern void print_grain_progress(unsigned int grain, unsigned int total_grains);

   } // end of internal namespace
} // end of grains namespace

#endif //GRAINS_INTERNAL_H_
