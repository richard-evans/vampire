//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "exchange.hpp"
#include "errors.hpp"
#include "vio.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for exchange module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="exchange";
      if(key!=prefix) return false;

      //-------------------------------------------------------------------
      std::string test="dmi-cutoff-range";
      if(word==test){
          double cr = atof(value.c_str());
          // Test for valid range
          vin::check_for_valid_value(cr, word, line, prefix, unit, "length", 0.0, 1.0e9,"input","0.0 - 1e9");
          internal::dmi_cutoff_range = cr;
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
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index, const int max_materials){

      // add prefix string
      std::string prefix="material:";

      // Check for empty material parameter array and resize to avoid segmentation fault
      if(internal::mp.size() == 0){
         internal::mp.resize(max_materials);
      }

      //------------------------------------------------------------
      // Check for material properties
      //------------------------------------------------------------
      std::string test = "dmi-constant"; // short form
      std::string test2 = "dzyaloshinskii-moriya-interaction-constant"; // long form
      if( (word == test) || (word == test2) ){
         double dmi = atof(value.c_str());
         vin::check_for_valid_value(dmi, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e17");
         internal::mp[super_index].dmi[sub_index] = dmi;
         internal::enable_dmi = true; // Switch on dmi calculation and fully unrolled tensorial anisotropy
         return true;
      }
      test = "exchange-matrix";
      if(word==test){
         // extract comma separated values from string
         std::vector<double> Jij = vin::doubles_from_string(value);
         if(Jij.size() == 1){
            vin::check_for_valid_value(Jij[0], word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e18");
            // set all components in case vectorial form is needed later
            vin::read_material[super_index].Jij_matrix_SI[sub_index][0] = Jij[0]; // Import exchange as field
            vin::read_material[super_index].Jij_matrix_SI[sub_index][1] = Jij[0];
            vin::read_material[super_index].Jij_matrix_SI[sub_index][2] = Jij[0];
            return true;
         }
         else if(Jij.size() == 3){
            vin::check_for_valid_vector(Jij, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e18");
            vin::read_material[super_index].Jij_matrix_SI[sub_index][0] = Jij[0]; // Import exchange as field
            vin::read_material[super_index].Jij_matrix_SI[sub_index][1] = Jij[1];
            vin::read_material[super_index].Jij_matrix_SI[sub_index][2] = Jij[2];
            // set vectorial anisotropy
            internal::exchange_type = internal::vectorial;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error in input file - material[" << super_index << "]:exchange_matrix[" << sub_index << "] must have one or three values." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error in input file - material[" << super_index << "]:exchange_matrix[" << sub_index << "] must have one or three values." << std::endl;
            err::vexit();
            return false;
         }
      }
      test = "biquadratic-exchange";
      if( word == test ){
         double bqe = atof(value.c_str());
         vin::check_for_valid_value(bqe, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e17");
         internal::mp[super_index].bqe[sub_index] = bqe;
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of exchange namespace
