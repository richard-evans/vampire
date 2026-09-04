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
#include <string>

// Vampire headers
#include "errors.hpp"
#include "grains.hpp"
#include "spininitialize.hpp"
#include "vio.hpp"
#include "random.hpp"

// grains module headers
#include "internal.hpp"

namespace grains{

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for grains module. All
   // keywords below still use the "create:" prefix, matching every other
   // create:* keyword handled by the create module.
   //-----------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="create";
      if(key!=prefix) return false;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      std::string test="voronoi-grain-substructure";
      if(word==test){
         internal::generate_voronoi_substructure = true;
         return true;
      }
      //--------------------------------------------------------------------
      // create:voronoi-size-variance is deprecated: despite its name it is a
      // standard deviation, not a variance, and maps onto
      // create:grain-size-sd. Does NOT set create:grain-regularity (the
      // unrelated seed-placement jitter), so warn the user explicitly.
      //--------------------------------------------------------------------
      test="voronoi-size-variance";
      if(word==test || word == "grain-size-variance"){
         double vsd=atof(value.c_str());
         vin::check_for_valid_value(vsd, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::grain_size_sd=vsd;
         internal::grain_size_shape_param_explicitly_set=true;
         zlog << zTs() << "Warning: create:" << word << " is deprecated and does not map 1:1 onto its replacement. "
              << "Setting create:grain-size-sd = " << vsd << "; if you relied on the old hexagonal-lattice jitter "
              << "behaviour, set create:grain-regularity explicitly as well." << std::endl;
         return true;
      }
      // create:voronoi-row-offset is deprecated and has no equivalent: grain
      // seeds are no longer placed on a hexagonal lattice. Warn and ignore.
      test="voronoi-row-offset";
      if(word==test){
         zlog << zTs() << "Warning: create:voronoi-row-offset is deprecated and has no effect - "
              << "the hexagonal seed lattice it offset has been retired." << std::endl;
         return true;
      }
      //--------------------------------------------------------------------
      test="voronoi-random-seed";
      if(word==test || word == "grain-structure-random-seed"){
         int vs=atoi(value.c_str());
         vin::check_for_valid_int(vs, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
         internal::grain_structure_seed=vs;
         // also set the seed for the grain substructure path (grain_substructure.cpp)
         mtrandom::voronoi_seed=vs;
         if(word == "voronoi-random-seed"){
            zlog << zTs() << "Warning: create:voronoi-random-seed is deprecated, use create:grain-structure-random-seed instead." << std::endl;
         }
         return true;
      }
      // create:voronoi-rounded-grains + create:voronoi-rounded-grains-area
      // are deprecated and map onto the single create:grain-rounding fraction.
      test="voronoi-rounded-grains";
      if(word==test || word == "rounded-grains"){
         internal::grain_rounding=0.9;
         zlog << zTs() << "Warning: create:" << word << " is deprecated, use create:grain-rounding = <area fraction> instead. "
              << "Setting create:grain-rounding = 0.9." << std::endl;
         return true;
      }
      // create:voronoi-include-boundary-grains is deprecated and has no
      // equivalent: every grain cell touching the domain boundary of a
      // non-periodic film is now unconditionally removed (see
      // remove_boundary_grains()). Warn and ignore.
      test="voronoi-include-boundary-grains";
      if(word==test || word == "include-boundary-grains"){
         zlog << zTs() << "Warning: create:" << word << " is deprecated and has no effect - "
              << "grain cells touching the film boundary are always removed for a non-periodic film "
              << "(leaving a void margin), rather than being clipped and kept." << std::endl;
         return true;
      }
      //-------------------------------------------------------------------
      test="voronoi-rounded-grains-area";
      if(word==test || word == "rounded-grains-area"){
         double vsd=atof(value.c_str());
         vin::check_for_valid_value(vsd, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::grain_rounding=vsd;
         zlog << zTs() << "Warning: create:" << word << " is deprecated, use create:grain-rounding instead." << std::endl;
         return true;
      }
      //-------------------------------------------------------------------
      test="voronoi-elliptical-rounding";
      if(word==test || word == "elliptical-rounding"){
         double er=atof(value.c_str());
         vin::check_for_valid_value(er, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::voronoi_elliptical_rounding=er;
         return true;
      }
      //-------------------------------------------------------------------
      test="voronoi-elliptical-rounding-height";
      if(word==test || word == "elliptical-rounding-height"){
         double erh=atof(value.c_str());
         vin::check_for_valid_value(erh, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::voronoi_elliptical_rounding_height=erh;
         return true;
      }
      // create:voronoi-bimodal-grains and its small-grain-* companions are
      // deprecated and map onto the general grain-size-distribution
      // mechanism below.
      test="voronoi-bimodal-grains";
      if(word==test){
         internal::grain_size_distribution=internal::grain_size_bimodal;
         zlog << zTs() << "Warning: create:voronoi-bimodal-grains is deprecated, use create:grain-size-distribution = bimodal instead." << std::endl;
         return true;
      }
      //-------------------------------------------------------------------
      test="voronoi-small-grain-diameter";
      if(word==test){
         double sgd=atof(value.c_str());
         vin::check_for_valid_value(sgd, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         internal::grain_size_second=sgd;
         internal::grain_size_shape_param_explicitly_set=true;
         zlog << zTs() << "Warning: create:voronoi-small-grain-diameter is deprecated, use create:grain-size-second instead." << std::endl;
         return true;
      }
      //-------------------------------------------------------------------
      test="voronoi-small-grain-fraction";
      if(word==test){
         double sgf=atof(value.c_str());
         vin::check_for_valid_value(sgf, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::grain_size_second_fraction=sgf;
         internal::grain_size_shape_param_explicitly_set=true;
         zlog << zTs() << "Warning: create:voronoi-small-grain-fraction is deprecated, use create:grain-size-second-fraction instead." << std::endl;
         return true;
      }
      //-------------------------------------------------------------------
      test="voronoi-small-grain-size-variance";
      if(word==test){
         double sgv=atof(value.c_str());
         vin::check_for_valid_value(sgv, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::grain_size_second_sd=sgv;
         internal::grain_size_shape_param_explicitly_set=true;
         zlog << zTs() << "Warning: create:voronoi-small-grain-size-variance is deprecated, use create:grain-size-second-sd instead." << std::endl;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-size";
      if(word==test){
         double gs=atof(value.c_str());
         vin::check_for_valid_value(gs, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         internal::grain_size=gs;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-size-sd";
      if(word==test){
         double gsd=atof(value.c_str());
         vin::check_for_valid_value(gsd, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::grain_size_sd=gsd;
         internal::grain_size_shape_param_explicitly_set=true;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-size-second";
      if(word==test){
         double gs2=atof(value.c_str());
         vin::check_for_valid_value(gs2, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         internal::grain_size_second=gs2;
         internal::grain_size_shape_param_explicitly_set=true;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-size-second-sd";
      if(word==test){
         double gsd2=atof(value.c_str());
         vin::check_for_valid_value(gsd2, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::grain_size_second_sd=gsd2;
         internal::grain_size_shape_param_explicitly_set=true;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-size-second-fraction";
      if(word==test){
         double gsf=atof(value.c_str());
         vin::check_for_valid_value(gsf, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::grain_size_second_fraction=gsf;
         internal::grain_size_shape_param_explicitly_set=true;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-size-distribution-file";
      if(word==test){
         internal::grain_size_distribution_file=value;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-periodic-boundaries";
      if(word==test){
         internal::grain_periodic_boundaries=true;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-spacing";
      if(word==test){
         double gsp=atof(value.c_str());
         vin::check_for_valid_value(gsp, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0 - 1 millimetre");
         internal::grain_spacing=gsp;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-spacing-sd";
      if(word==test){
         double gssd=atof(value.c_str());
         vin::check_for_valid_value(gssd, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::grain_spacing_sd=gssd;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-rounding";
      if(word==test){
         double gr=atof(value.c_str());
         vin::check_for_valid_value(gr, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::grain_rounding=gr;
         return true;
      }
      // create:grain-boundary-roughness clamps to 0.6 rather than 1.0: an
      // amplitude approaching 1.0 would let the perturbed radius approach
      // zero in some directions, producing extreme, near-degenerate concavity.
      test="grain-boundary-roughness";
      if(word==test){
         double gbr=atof(value.c_str());
         vin::check_for_valid_value(gbr, word, line, prefix, unit, "none", 0.0, 0.6,"input","0.0 - 0.6");
         internal::grain_boundary_roughness=gbr;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-boundary-roughness-modes";
      if(word==test){
         int gbrm=atoi(value.c_str());
         vin::check_for_valid_int(gbrm, word, line, prefix, 1, 20,"input","1 - 20");
         internal::grain_boundary_roughness_modes=gbrm;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-regularity";
      if(word==test){
         double gr=atof(value.c_str());
         vin::check_for_valid_value(gr, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::grain_regularity=gr;
         return true;
      }
      // Wulff faceting (create:grain-facet-*). grain-facet-symmetry is the
      // master off-switch (0, default); the valid non-zero values (4, 6, 8)
      // are even facet counts, since the 2D Wulff construction alternates
      // an anisotropic facet energy between two interleaved facet families.
      test="grain-facet-symmetry";
      if(word==test){
         int gfs=atoi(value.c_str());
         vin::check_for_valid_int(gfs, word, line, prefix, 0, 8,"input","0, 4, 6 or 8");
         if(gfs!=0 && gfs!=4 && gfs!=6 && gfs!=8){
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'" << prefix << ":" << word << "\' must be one of 0 (off), 4, 6 or 8." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error - value for \'" << prefix << ":" << word << "\' must be one of 0 (off), 4, 6 or 8." << std::endl;
            err::vexit();
         }
         internal::grain_facet_symmetry=gfs;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-facet-strength";
      if(word==test){
         double gfst=atof(value.c_str());
         vin::check_for_valid_value(gfst, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
         internal::grain_facet_strength=gfst;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-facet-anisotropy";
      if(word==test){
         double gfa=atof(value.c_str());
         vin::check_for_valid_value(gfa, word, line, prefix, unit, "none", 0.0, 10.0,"input","0.0 - 10.0");
         internal::grain_facet_anisotropy=gfa;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-facet-orientation-angle";
      if(word==test){
         double gfoa=atof(value.c_str());
         vin::check_for_valid_value(gfoa, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
         internal::grain_facet_orientation_angle=gfoa;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-facet-orientation-spread";
      if(word==test){
         double gfos=atof(value.c_str());
         vin::check_for_valid_value(gfos, word, line, prefix, unit, "none", 0.0, 180.0,"input","0.0 - 180.0 degrees");
         internal::grain_facet_orientation_spread=gfos;
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-facet-orientation";
      if(word==test){
         test="random";
         if(value==test){
            internal::grain_facet_orientation = internal::facet_orientation_random;
            return true;
         }
         test="fixed";
         if(value==test){
            internal::grain_facet_orientation = internal::facet_orientation_fixed;
            return true;
         }
         test="textured";
         if(value==test){
            internal::grain_facet_orientation = internal::facet_orientation_textured;
            return true;
         }
         // otherwise throw an error
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'" << prefix << ":" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"random\"" << std::endl;
         std::cerr << "\t\"fixed\"" << std::endl;
         std::cerr << "\t\"textured\"" << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error - value for \'" << prefix << ":" << word << "\' must be one of:" << std::endl;
         zlog << zTs() << "\t\"random\"" << std::endl;
         zlog << zTs() << "\t\"fixed\"" << std::endl;
         zlog << zTs() << "\t\"textured\"" << std::endl;
         err::vexit();
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-statistics-output";
      if(word==test){
         internal::output_grain_statistics_file=true; // default
         std::string VFalse="false";
         if(value==VFalse){
            internal::output_grain_statistics_file=false;
         }
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-size-distribution";
      if(word==test){
         test="lognormal";
         if(value==test){
            internal::grain_size_distribution = internal::grain_size_lognormal;
            return true;
         }
         test="normal";
         if(value==test){
            internal::grain_size_distribution = internal::grain_size_normal;
            return true;
         }
         test="bimodal";
         if(value==test){
            internal::grain_size_distribution = internal::grain_size_bimodal;
            return true;
         }
         test="delta";
         if(value==test){
            internal::grain_size_distribution = internal::grain_size_delta;
            return true;
         }
         test="file";
         if(value==test){
            internal::grain_size_distribution = internal::grain_size_file;
            return true;
         }
         // otherwise throw an error
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'" << prefix << ":" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"lognormal\"" << std::endl;
         std::cerr << "\t\"normal\"" << std::endl;
         std::cerr << "\t\"bimodal\"" << std::endl;
         std::cerr << "\t\"delta\"" << std::endl;
         std::cerr << "\t\"file\"" << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error - value for \'" << prefix << ":" << word << "\' must be one of:" << std::endl;
         zlog << zTs() << "\t\"lognormal\"" << std::endl;
         zlog << zTs() << "\t\"normal\"" << std::endl;
         zlog << zTs() << "\t\"bimodal\"" << std::endl;
         zlog << zTs() << "\t\"delta\"" << std::endl;
         zlog << zTs() << "\t\"file\"" << std::endl;
         err::vexit();
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-shape-output";
      if(word==test){
         internal::output_gv_file=true; // default
         // also check for value
         std::string VFalse="false";
         if(value==VFalse){
            internal::output_gv_file=false;
         }
         return true;
      }
      //-------------------------------------------------------------------
      test="grain-tessellation";
      if(word==test){
         test="voronoi"; // half-plane clipping engine with uniform weights - an ordinary Voronoi diagram
         if(value==test){
            internal::grain_tessellation = internal::tessellation_voronoi;
            return true;
         }
         test="laguerre"; // half-plane clipping engine with real per-seed weights (w_i = r_i^2) - a Laguerre/power diagram
         if(value==test){
            internal::grain_tessellation = internal::tessellation_laguerre;
            return true;
         }
         // otherwise throw an error
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'" << prefix << ":" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"voronoi\"" << std::endl;
         std::cerr << "\t\"laguerre\"" << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error - value for \'" << prefix << ":" << word << "\' must be one of:" << std::endl;
         zlog << zTs() << "\t\"voronoi\"" << std::endl;
         zlog << zTs() << "\t\"laguerre\"" << std::endl;
         err::vexit();
         return true;
      }
      //--------------------------------------------------------------------
      // create:grain-magnetisation-direction = material | alternating
      //
      // Controls how the initial spin direction (material[#]:initial-spin-
      // direction) is applied across grains: "material" (default) leaves
      // every atom as specified by its material's texture; "alternating"
      // reverses the spin direction for every odd-numbered grain, giving
      // neighbouring grains opposite magnetisation directions.
      //--------------------------------------------------------------------
      test="grain-magnetisation-direction";
      if(word==test){
         std::string loctest="alternating";
         if(value==loctest){
            spininitialize::set_grain_magnetisation_mode(1); // grain_mode_alternating
            return true;
         }
         else{
            // default: "material" (or any unrecognised value falls back to
            // the default of no grain-level post-processing)
            spininitialize::set_grain_magnetisation_mode(0); // grain_mode_material
            return true;
         }
      }
      //--------------------------------------------------------------------
      test="voronoi-grain-substructure-crystallization-radius";
      if(word==test){
         double rsize=atof(value.c_str());
         vin::check_for_valid_value(rsize, word, line, prefix, unit, "none", 0.01, 2.0,"input","0.01 - 2");
         internal::voronoi_grain_substructure_crystallization_radius=rsize;
         return true;
      }
      //--------------------------------------------------------------------
      test="voronoi-grain-substructure-overlap-factor";
      if(word==test){
         double ol=atof(value.c_str());
         vin::check_for_valid_value(ol, word, line, prefix, unit, "none", 0.1, 3.0,"input","0.1 - 3");
         internal::voronoi_grain_substructure_overlap_factor = ol;
         return true;
      }
      //--------------------------------------------------------------------
      test="voronoi-grain-substructure-size";
      if(word==test){
         double psize=atof(value.c_str());
         vin::check_for_valid_value(psize, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         internal::voronoi_grain_substructure_size=psize;
         return true;
      }
      else
      //--------------------------------------------------------------------
      test="voronoi-grain-substructure-spacing";
      if(word==test){
         double pspacing=atof(value.c_str());
         vin::check_for_valid_value(pspacing, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.0 Angstroms - 1 millimetre");
         internal::voronoi_grain_substructure_spacing=pspacing;
         return true;
      }
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

   //---------------------------------------------------------------------------
   // Function to process material parameters for grains module
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      // No material-level keywords belong to this module: the one
      // grain-substructure material parameter that might seem to
      // (create:material[#]:voronoi-grain-substructure-nucleation-height)
      // sets a field on create::internal::mp_t, which lives in the create
      // module (shared with every particle-shape file), so it is parsed by
      // create::match_material_parameter instead.
      (void)word; (void)value; (void)unit; (void)line; (void)super_index; (void)sub_index;

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of grains namespace
