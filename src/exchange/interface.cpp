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

   // internal namespacve for exchange module
   namespace internal{

      //-----------------------------------------------------------------------------------------------------------------------
      // general function for reading exchange value from input file and storing in 4D format
      //-----------------------------------------------------------------------------------------------------------------------
      void read_exchange_values(int material_i, int material_j, int neighbour, std::string const word, std::string const prefix, std::string const value, std::string const unit, int const line, exchange_matrix_4D_t& exchange_matrix){

         // extract comma separated values from string
         std::vector<double> Jij = vin::doubles_from_string(value);

         // optional factor 2 correction for ab-initio
         const double ef = internal::exchange_factor;
         for(auto& J : Jij) J *= ef;

         if(Jij.size() == 1){
            vin::check_for_valid_value(Jij[0], word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e18");
            // set exchange constants
            exchange_matrix.set_exchange_values(material_i, material_j, neighbour, Jij);
         }
         else if(Jij.size() == 3){
            vin::check_for_valid_vector(Jij, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e18");
            // set exchange constants
            exchange_matrix.set_exchange_values(material_i, material_j, neighbour, Jij);
            internal::minimum_needed_exchange_type = exchange::vectorial;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error in input file - material[" << material_i << "]:exchange_matrix[" << material_j << "] must have one or three values." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error in input file - material[" << material_i << "]:exchange_matrix[" << material_j << "] must have one or three values." << std::endl;
            err::vexit();
         }

         return;

      }

   } // end of internal namespace

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
          vin::check_for_valid_value(cr, word, line, prefix, unit, "length", 0.0, 1.0e3,"input","0.0 - 1e3");
          internal::dmi_cutoff_range = cr;
          return true;
      }
      //-------------------------------------------------------------------
      test="kitaev-cutoff-range";
      if(word==test){
          double cr = atof(value.c_str());
          // Test for valid range
          vin::check_for_valid_value(cr, word, line, prefix, unit, "length", 0.0, 1e3,"input","0.0 - 1e3");
          internal::kitaev_cutoff_range = cr;
          return true;
      }
      //--------------------------------------------------------------------
      if( word == "ab-initio" ){
         bool set_exchange = vin::check_for_valid_bool(value, word, line, prefix,"input");
         if(set_exchange) internal::exchange_factor = 2.0;
         else internal::exchange_factor = 1.0;
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
      //------------------------------------------------------------------------
      test = "exchange-matrix";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 0, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-1st-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 0, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-2nd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 1, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-3rd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 2, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-4th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 3, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-5th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 4, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-6th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 5, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-7th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 6, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-8th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 6, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-9th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 8, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-10th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 9, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 0, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-1st-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 0, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-2nd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 1, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-3rd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 2, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-4th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 3, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-5th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 4, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-6th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 5, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-7th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 6, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-8th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 6, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-9th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 8, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-10th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 9, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //--------------------------------------------------------------------
      test = "kitaev-constant"; // short form
      if( word == test ){
         double k = atof(value.c_str());
         vin::check_for_valid_value(k, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e17");
         internal::mp[super_index].kitaev[sub_index] = k;
         internal::enable_kitaev = true; // Switch on kitaev calculation and fully unrolled tensorial anisotropy
         return true;
      }
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of exchange namespace
