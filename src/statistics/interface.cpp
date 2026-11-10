//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "stats.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "vio.hpp"

// statistics module headers
#include "internal.hpp"

namespace stats{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for statistics module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="output";
      if(key!=prefix) return false;

      //--------------------------------------------------------------------
      std::string test="geofencing-cubic-threshold"; // magnetisation is in a cubic state if m . e_a >= threshold
      if(word==test){
         double threshold = atof(value.c_str());
         vin::check_for_valid_value(threshold, word, line, prefix, unit, "none", 0.0, 1.0, "input", "0.0 - 1.0");
         internal::geofencing[cubic_geofencing].threshold = threshold;
         return true;
      }
      //--------------------------------------------------------------------
      test="geofencing-cubic-mode"; // <111> easy axes or <100> easy axes when <111> directions are hard
      if(word==test){
         if(value == "111-easy") internal::cubic_geofencing_111_hard = false;
         else if(value == "111-hard") internal::cubic_geofencing_111_hard = true;
         else err::fatal("Error - value for \'" + prefix + ":" + word + "\' must be \"111-easy\" or \"111-hard\" on line " + std::to_string(line) + " of input file");
         return true;
      }
      //--------------------------------------------------------------------
      test="geofencing-cubic-material-id"; // geofence material magnetisation instead of system magnetisation
      if(word==test){
         int mat = atoi(value.c_str());
         vin::check_for_valid_int(mat, word, line, prefix, 1, mp::max_materials, "input", "1 - 100");
         internal::geofencing[cubic_geofencing].material = mat - 1; // save material ID as m-1
         stats::calculate_material_magnetization = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="geofencing-uniaxial-threshold"; // magnetisation is in a uniaxial state if m . e_a >= threshold
      if(word==test){
         double threshold = atof(value.c_str());
         vin::check_for_valid_value(threshold, word, line, prefix, unit, "none", 0.0, 1.0, "input", "0.0 - 1.0");
         internal::geofencing[uniaxial_geofencing].threshold = threshold;
         return true;
      }
      //--------------------------------------------------------------------
      test="geofencing-uniaxial-material-id"; // geofence material magnetisation instead of system magnetisation
      if(word==test){
         int mat = atoi(value.c_str());
         vin::check_for_valid_int(mat, word, line, prefix, 1, mp::max_materials, "input", "1 - 100");
         internal::geofencing[uniaxial_geofencing].material = mat - 1; // save material ID as m-1
         stats::calculate_material_magnetization = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="geofencing-uniaxial-axis"; // easy axis for uniaxial geofencing
      if(word==test){
         std::vector<double> u = vin::doubles_from_string(value);
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input"); // also normalises u
         internal::uniaxial_geofencing_axis = u;
         return true;
      }
      //--------------------------------------------------------------------
      test="geofencing-per-spin-thresholds"; // list of m . e_a thresholds for per-spin geofencing
      if(word==test){
         std::vector<double> thresholds = vin::doubles_from_string(value);
         if(thresholds.empty()) err::fatal("Error - value for \'" + prefix + ":" + word + "\' must include at least one threshold on line " + std::to_string(line) + " of input file");
         for(double& threshold : thresholds) vin::check_for_valid_value(threshold, word, line, prefix, unit, "none", 0.0, 1.0, "input", "0.0 - 1.0");
         internal::per_spin_geofencing_thresholds = thresholds;
         return true;
      }

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      // add prefix string
      std::string prefix="material:";

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of stats namespace
