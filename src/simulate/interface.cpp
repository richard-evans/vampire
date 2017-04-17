//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "errors.hpp"
#include "sim.hpp"
#include "vio.hpp"

// Internal sim header
#include "internal.hpp"

namespace sim{

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for sim module
   //-----------------------------------------------------------------------------
   bool match_input_parameter(string const key, string const word, string const value, string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="sim";
      if(key!=prefix) return false;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------

      std::string test="slonczewski-spin-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::DoublesFromString(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::slonczewski_spin_polarization_unit_vector = u;
         return true;
      }
      test="monte-carlo-preconditioning-steps";
      if(word==test){
         int n = atoi(value.c_str());
         // Test for valid range
         vin::check_for_valid_int(n, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
         sim::internal::num_monte_carlo_preconditioning_steps = n;
         return true;
      }
      //--------------------------------------------------------------------
      // input parameter not found here
      return false;
   }

   //----------------------------------------------------------------------------------
   // material parameter match function
   //----------------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index){

      // add prefix string
      std::string prefix="material:";

      // Check for material id > current array size and if so dynamically expand mp array
      if((unsigned int) super_index + 1 > sim::internal::mp.size() && super_index + 1 < 101) sim::internal::mp.resize(super_index + 1);

      //------------------------------------------------------------
      std::string test="slonczewski-adiabatic-spin-torque";
      /*
         aj parameter for material in slonczewski torque calculation
         */
      if(word==test){
         double aj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(aj, word, line, prefix, unit, "field", 0.0, 1.0e2,"input","0 - 100T");
         sim::internal::mp[super_index].slonczewski_aj.set(aj);
         return true;
      }
      //------------------------------------------------------------
      test="slonczewski-non-adiabatic-spin-torque";
      /*
         bj parameter for material in slonczewski torque calculation
         */
      if(word==test){
         double bj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(bj, word, line, prefix, unit, "field", 0.0, 1.0e2,"input","0 - 100T");
         sim::internal::mp[super_index].slonczewski_bj.set(bj);
         return true;
      }
      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of namespace sim
