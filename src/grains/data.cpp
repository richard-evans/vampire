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

// Vampire headers
#include "grains.hpp"

// grains module headers
#include "internal.hpp"

namespace grains{

   namespace internal{

      //----------------------------------------------------------------------------
      // Shared variables used within grains module
      //----------------------------------------------------------------------------

      // Own random number generator, reseeded from grain_structure_seed at
      // the start of every voronoi_film()/voronoi_substructure() call.
      MTRand grnd;

      bool generate_voronoi_substructure = false;
      double voronoi_grain_substructure_size = 50.0;    // mean grain size of the substructure within a particle
      double voronoi_grain_substructure_spacing = 10.0; // spacing between substructure grains
      double voronoi_grain_substructure_crystallization_radius = 1.2;
      double voronoi_grain_substructure_overlap_factor = 1.0;

      // Default hexagonal-lattice jitter for the grain substructure path
      // (grain_substructure.cpp); 0.1 is the standard deviation of seed
      // displacement as a fraction of the lattice spacing.
      bool parity = 0;
      double voronoi_sd = 0.1;

      double voronoi_elliptical_rounding = 0.0;        // 0.0 = vertical grain walls (default), 1.0 = fully rounded cap
      double voronoi_elliptical_rounding_height = 0.5; // fraction of grain_film_height where the flat column ends and the cap begins
      double grain_film_height = -1.0;                 // < 0 (unset) falls back to cs::system_dimensions[2]

      bool output_gv_file = false; // toggle output of grain positions to file (create:grain-shape-output)

      grain_tessellation_t grain_tessellation = tessellation_laguerre; // create:grain-tessellation, defaults to the weighted Laguerre engine

      //-------------------------------------------------------------------
      // Grain size distribution and weight fitting (create:grain-size-*)
      //-------------------------------------------------------------------
      double grain_size = -1.0; // sentinel: fall back to dimensions:particle-size
      double grain_spacing = -1.0; // sentinel: fall back to dimensions:particle-spacing
      // Default to a normal distribution with a modest 0.3 sd rather than
      // delta (monodisperse): this gives a naturally varied-looking film
      // out of the box, without requiring a user who never touches
      // create:grain-size-distribution/grain-size-sd to get perfectly
      // uniform grains.
      grain_size_distribution_t grain_size_distribution = grain_size_normal;
      double grain_size_sd = 0.3;
      double grain_size_second = 20.0; // 2 nm in Angstroms
      double grain_size_second_sd = 0.15;
      double grain_size_second_fraction = 0.3;
      std::string grain_size_distribution_file = "";
      bool grain_size_shape_param_explicitly_set = false;
      double grain_regularity = 0.0;
      int grain_structure_seed = 1951218893; // matches random.hpp's legacy voronoi_seed default
      bool output_grain_statistics_file = false; // create:grain-statistics-output
      double grain_rounding = 0.0; // create:grain-rounding, 0 = off (default)
      double grain_boundary_roughness = 0.0; // create:grain-boundary-roughness, 0 = off (default, plain circular rounding)
      int grain_boundary_roughness_modes = 4; // create:grain-boundary-roughness-modes
      double grain_spacing_sd = 0.0; // create:grain-spacing-sd, 0 = off (default, constant-width boundary)
      bool grain_periodic_boundaries = false; // create:grain-periodic-boundaries

      //-------------------------------------------------------------------
      // Wulff faceting (create:grain-facet-*)
      //-------------------------------------------------------------------
      int grain_facet_symmetry = 0; // create:grain-facet-symmetry, 0 = off (default)
      double grain_facet_strength = 1.0; // create:grain-facet-strength; only reached once symmetry > 0
      double grain_facet_anisotropy = 1.0; // isotropic (regular polygon) by default
      grain_facet_orientation_t grain_facet_orientation = facet_orientation_random;
      double grain_facet_orientation_angle = 0.0;  // degrees
      double grain_facet_orientation_spread = 15.0; // degrees, textured mode only

   } // end of internal namespace

} // end of grains namespace
