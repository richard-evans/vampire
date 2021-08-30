//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "micromagnetic.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "units.hpp"
#include "errors.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for micromagnetic module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="micromagnetic";
      if(key!=prefix) return false;

      //--------------------------------------------------------------------
      std::string test="integrator";
      if(word==test){
         test="llg";
         if(value==test){
            micromagnetic::integrator=0;
            return true;
         }
         test="llb";
         if(value==test){
            micromagnetic::integrator=1;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"llg\"" << std::endl;
            std::cerr << "\t\"llb\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //--------------------------------------------------------------------
      test="atomistic-steps-per-micromagnetic-step";
      if(word==test){
         double dt=atof(value.c_str());
         vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1, 100000,"input","1 step - 100000 steps");
         micromagnetic::num_atomic_steps_mm =dt;
         return true;
      }
      //--------------------------------------------------------------------
      test="discretisation"; // whether the material is micromagnetic or atomistic
      if(word==test){

         // check for type of discretisation
         test="micromagnetic";
         if(value==test){
            discretisation_type = 1;
            return true;
         }
         test="atomistic"; // runs simualtion as normal (atomistics)
         if(value==test){
            discretisation_type = 0;
            return true;
         }
         test="multiscale"; // at the moment just runs a normal atomsitic simulation
         if(value==test){
            discretisation_type = 2;
            return true;
         }
         // otherwise throw an error
         terminaltextcolor(RED);
         std::cerr << "Error - value for" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"micromagnetic\"" << std::endl;
         std::cerr << "\t\"atomistic\"" << std::endl;
         std::cerr << "\t\"multiscale\"" << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error - value for" << word << "\' must be one of:" << std::endl;
         zlog << zTs() << "\t\"micromagnetic\"" << std::endl;
         zlog << zTs() <<"\t\"atomistic\"" << std::endl;
         zlog << zTs() << "\t\"multiscale\"" << std::endl;
         err::vexit();

         return true;
      }
      //--------------------------------------------------------------------
      test="pinning-field-correction";
      if(word==test){
         micromagnetic::internal::mm_correction = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="resistance-GMR";
      if(word==test){
         double h = atof(value.c_str());
         vin::check_for_valid_value(h, word, line, prefix, unit, "none", 0, 100000,"input","0 - 100000");
         micromagnetic::internal::res_GMR = h;
         return true;
      }
      //--------------------------------------------------------------------
      test="resistance-RA";
      if(word==test){
         double h = atof(value.c_str());
         vin::check_for_valid_value(h, word, line, prefix, unit, "none", 0, 100000,"input","0 - 100000");
         micromagnetic::internal::res_RA = h;
         return true;
      }
      //--------------------------------------------------------------------
      test="bias-magnets";
      if(word==test){
         micromagnetic::internal::bias_magnets = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="bias-magnets-max-height";
      if(word==test){
         double h = atof(value.c_str());
         micromagnetic::internal::bias_magnets_max_height = h;
         vin::check_for_valid_value(h, word, line, prefix, unit, "none", 0, 1,"input","0 - 1");
         return true;
      }
      //--------------------------------------------------------------------
      test="bias-magnets-min-height";
      if(word==test){
         double h = atof(value.c_str());
         micromagnetic::internal::bias_magnets_min_height = h;
         vin::check_for_valid_value(h, word, line, prefix, unit, "none", 0, 1,"input","0 - 1");
         return true;
      }
      //--------------------------------------------------------------------
      test="bias-magnets-max-width";
      if(word==test){
         double h = atof(value.c_str());
         micromagnetic::internal::bias_magnets_max_width = h;
         vin::check_for_valid_value(h, word, line, prefix, unit, "none", 0, 1,"input","0 - 1");
         return true;
      }
      //--------------------------------------------------------------------
      test="bias-magnets-min-width";
      if(word==test){
         double h = atof(value.c_str());
         micromagnetic::internal::bias_magnets_min_width = h;
         vin::check_for_valid_value(h, word, line, prefix, unit, "none", 0, 1,"input","0 - 1");
         return true;
      }
      //--------------------------------------------------------------------
      test="bias-magnets-gap";
      if(word==test){
         double h = atof(value.c_str());
         vin::check_for_valid_value(h, word, line, prefix, unit, "length", 0, 100000000,"input","1 A - 100000 A");
         micromagnetic::internal::bias_magnets_gap = h;
         return true;
      }
      //--------------------------------------------------------------------
      test="pinning-field-height";
      if(word==test){
         double h=atof(value.c_str());
         vin::check_for_valid_value(h, word, line, prefix, unit, "length", 1, 100000000,"input","1 A - 100000 A");
         micromagnetic::internal::pinning_field_height =h;
         return true;
      }
      //--------------------------------------------------------------------
      test="bias-magnet-Ms";
      if(word==test){
         double h=atof(value.c_str());
         vin::check_for_valid_value(h, word, line, prefix, unit, "mm_energy", -100000000, 100000000,"input","-100000000 - 100000");
         micromagnetic::internal::bias_magnet_ms_input =h;
         return true;
      }
      //--------------------------------------------------------------------
      test="temperature-dependent-parameters";
      if(word==test){
         bool tf = vin::check_for_valid_bool(value, word, line, prefix, "input");
         micromagnetic::internal::temperature_dependent_parameters = tf;
         return true;
      }
      test="output-magnetisation";
      if(word==test){
         micromagnetic::internal::output_m = true;
         micromagnetic::internal::output_list.push_back(0);
         return true;
      }
      test="output-applied-field";
      if(word==test){
         micromagnetic::internal::output_applied_field = true;
         micromagnetic::internal::output_list.push_back(1);

         return true;
      }

            test="output-resistance";
      if(word==test){
         micromagnetic::internal::output_resistance = true;
         micromagnetic::internal::output_list.push_back(4);

         return true;
      }

      test="output-temperature";
      if(word==test){
         micromagnetic::internal::output_temperature = true;
         micromagnetic::internal::output_list.push_back(2);

         return true;
      }
      test="output-time";
      if(word==test){
         micromagnetic::internal::output_time = true;
         micromagnetic::internal::output_list.push_back(3);

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
      std::string test="micromagnetic-discretisation-enabled";
      if(word==test){
         vin::read_material[super_index].micromagnetic_enabled=true;
         //std::cout << super_index << '\t' << vin::read_material[super_index].micromagnetic_enabled << std::endl;
         return true;
      }


      test="pinning-field-strength";
      if(word==test){
         double K=atof(value.c_str());
         vin::check_for_valid_value(K, word, line, prefix, unit, "mm_energy", 1e-28, 10*1.0e-10,"material"," 0.000000 - 1");
         vin::read_material[super_index].pinning_field_strength = K;
         return true;
      }
      test="pinning-field-unit-vector";
      if(word==test){
         // temporary storage container
         std::vector<double> u(3);

         // read values from string
         u=vin::doubles_from_string(value);
         vin::check_for_valid_unit_vector(u, word, line, prefix,"length");
         // Copy sanitised unit vector to material
         vin::read_material[super_index].pinning_field_unit_vector[0] =u.at(0);
         vin::read_material[super_index].pinning_field_unit_vector[1] =u.at(1);
         vin::read_material[super_index].pinning_field_unit_vector[2] =u.at(2);
         return true;
      }

      test="SAF-exchange-coupling"; // new and preferred form
      if( (word == test)){
         double J = atof(value.c_str());
         vin::check_for_valid_value(J, word, line, prefix, unit, "mm_energy", 0.0000000000000001*1.0e-3/1.0e16, 10*1.0e-3/1.0e16,"material"," 0.000000 - 1");
         vin::read_material[super_index].SAF[sub_index] = J;
         vin::read_material[super_index].enable_SAF = true;
         vin::read_material[sub_index].enable_SAF = true;
         std::cout << super_index << '\t' << sub_index << std::endl;
         return true;
      }

      test="micromagnetic-exchange"; // new and preferred form
      if( (word == test)){
         double J = atof(value.c_str());
         vin::check_for_valid_value(J, word, line, prefix, unit, "exchange", 1e-56, 10e16,"material"," 0.000000 - 1");
         vin::read_material[super_index].EF_MM[sub_index] = J;
         vin::read_material[super_index].override_atomsitic[sub_index] = true;
         return true;
      }


      test="resistance"; // new and preferred form
      if( (word == test)){
        internal::resistance_layer_1 = sub_index;
        internal::resistance_layer_2 = super_index;
        enable_resistance = true;
         return true;
      }


      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of micromagnetic namespace
