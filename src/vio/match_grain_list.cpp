//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <sstream>
#include <iostream>

// vampire headers
#include "stats.hpp"
#include "vio.hpp"

// vio module headers
#include "internal.hpp"

namespace vout{

   // funtion to process input parameters for grain output file
   int match_vout_grain_list(std::string const word, std::string const value, int const line, std::vector<grain::output_t> & output_list){

      std::string prefix="grain:";

      std::string test="time-steps";
      if(word==test){
         grain::output_list.push_back(grain::time_steps);
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      test="real-time";
      if(word==test){
         grain::output_list.push_back(grain::real_time);
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      test="temperature";
      if(word==test){
         grain::output_list.push_back(grain::temperature);
         return EXIT_SUCCESS;
      }
      //-------------------------------------------------------------------
      test="electron-temperature"; // identical to temperature
      if(word==test){
         grain::output_list.push_back(grain::temperature);
         return EXIT_SUCCESS;
      }
      //-------------------------------------------------------------------
      test="phonon-temperature";
      if(word==test){
         grain::output_list.push_back(grain::phonon_temperature);
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      test="applied-field-strength";
      if(word==test){
         grain::output_list.push_back(grain::applied_field);
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      test="applied-field-unit-vector";
      if(word==test){
         grain::output_list.push_back(grain::applied_field_unit_vector);
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      test="magnetisation";
      if(word==test){
         // enable statistics collection for grains
         stats::calculate_grain_magnetization = true;
         grain::output_list.push_back(grain::magnetisation);
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      test="mean-magnetisation-length";
      if(word==test){
         // enable statistics collection for grains
         stats::calculate_grain_magnetization = true;
         grain::output_list.push_back(grain::mean_magnetisation_length);
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      test="mean-susceptibility";
      if(word==test){
         // enable statistics collection for grains
         stats::calculate_grain_magnetization  = true;
         stats::calculate_grain_susceptibility = true;
         grain::output_list.push_back(grain::mean_susceptibility);
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      test="mean-specific-heat";
      if(word==test){
         // enable statistics collection for grains
         stats::calculate_grain_energy        = true;
         stats::calculate_grain_specific_heat = true;
         grain::output_list.push_back(grain::mean_specific_heat);
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      test="mean-torque";
      if(word==test){
         // enable statistics collection for grains
         stats::calculate_grain_torque = true;
         grain::output_list.push_back(grain::mean_torque);
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      test="material-magnetisation";
      if(word==test){
         // enable statistics collection for grains
         stats::calculate_material_grain_magnetization  = true;
         grain::output_list.push_back(grain::material_magnetisation);
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      test="material-height-magnetisation";
      if(word==test){
         // enable statistics collection for grains
         stats::calculate_material_grain_height_magnetization  = true;
         grain::output_list.push_back(grain::material_height_magnetisation);
         return EXIT_SUCCESS;
      }
      //-------------------------------------------------------------------
      test="output-rate";
      if(word==test){
         int r=atoi(value.c_str());
         vin::check_for_valid_int(r, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
         vout::grain::output_rate = r;
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      else{
         terminaltextcolor(RED);
         std::cerr << "Error - Unknown control statement \'grain:" << word << "\' on line " << line << " of input file" << std::endl;
         terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }

      return EXIT_SUCCESS;

   }

}
