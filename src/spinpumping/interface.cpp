//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard Evans 2022. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "spinpumping.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "vio.hpp"

// spinpumping module headers
#include "internal.hpp"

namespace spin_pumping{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for spinpumping module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="spin-pumping";
      if(key!=prefix) return false;

      //------------------------------------------------------------------------
      // If spin-transport parameter is requested, then enable module
      //------------------------------------------------------------------------
      spin_pumping::internal::enabled = true;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      //------------------------------------------------------------------------
      std::string test = "cell-size-x";
      if( word == test ){
          // Set spin transport cell size along x-direction
          double csx = vin::str_to_double(value);
          vin::check_for_valid_value(csx, word, line, prefix, unit, "length", 1.0, 10000.0,"input","0.1 - 1000 nm");
          internal::cell_size_x = csx;
          return true;
      }
      //------------------------------------------------------------------------
      test = "cell-size-y";
      if( word == test ){
          // Set spin transport cell size along y-direction
          double csy = vin::str_to_double(value);
          vin::check_for_valid_value(csy, word, line, prefix, unit, "length", 1.0, 10000.0,"input","0.1 - 1000 nm");
          internal::cell_size_y = csy;
          return true;
      }
      //------------------------------------------------------------------------
      test = "cell-size-z";
      if( word == test ){
          // Set spin transport cell size along z-direction
          double csz = vin::str_to_double(value);
          vin::check_for_valid_value(csz, word, line, prefix, unit, "length", 1.0, 10000.0,"input","0.1 - 1000 nm");
          internal::cell_size_z = csz;
          return true;
      }
      //------------------------------------------------------------------------
      test = "cell-size";
      if( word == test ){
          // Set spin transport cell sizes along x,y,z-directions
          double cs = vin::str_to_double(value);
          vin::check_for_valid_value(cs, word, line, prefix, unit, "length", 1.0, 10000.0,"input","0.1 - 1000 nm");
          internal::cell_size_x = cs;
          internal::cell_size_y = cs;
          internal::cell_size_z = cs;
          return true;
      }
      //------------------------------------------------------------------------
      test = "update-rate";
      if( word == test ){
          // Set resistivity for environment (cells with no atoms)
          uint64_t ur = vin::str_to_uint64(value);
          vin::check_for_valid_int(ur, word, line, prefix, 1,10000000,"input","1 - 1E7 time steps");
          spin_pumping::internal::update_rate = ur;
          // set time counter to ensure initial calculation at start
          spin_pumping::internal::time_counter = ur;
          return true;
      }
      //------------------------------------------------------------------------
      test = "atomistic-spin-pumping";
      if( word == test ){
          // Set resistivity for environment (cells with no atoms)
          sim::compute_time_derivative = true;
          spin_pumping::internal::output_atomistic_spin_pumping_flag = true;
          return true;
      }
      //------------------------------------------------------------------------
      test = "cell-spin-pumping";
      if( word == test ){
          // Set resistivity for environment (cells with no atoms)
          sim::compute_time_derivative = true;
          spin_pumping::internal::output_cells_spin_pumping_flag = true;
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

      // Check for material id > current array size and if so dynamically expand mp array
      if((unsigned int) super_index + 1 > internal::mp.size() && super_index + 1 < 101) internal::mp.resize(super_index + 1);

      //---------------------------------------------------------------------------
      // Effective spin mixing conductance [1/m^2] = hbar/(e^2) G
      // G := spin-mixing conductance [ 1/Ohm 1/m^2 ], where h/e^2 = Ohm
      std::string test = "spin-mixing-conductance";
      if( word == test ){
         // Set resistivity for atom type
         double geff = vin::str_to_double(value);
         vin::check_for_valid_value(geff, word, line, prefix, unit, "spin_mixing_conductance", 1.0e-20, 1.0e20,"input","1E-20 - 1E20 per Ohm per metre squared");
         spin_pumping::internal::mp[super_index].spin_mix_conductance.set(geff);
         return true;
      }

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of spin_pumping namespace
