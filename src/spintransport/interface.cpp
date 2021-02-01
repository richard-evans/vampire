//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "spintransport.hpp"
#include "errors.hpp"
#include "vio.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for spintransport module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="spin-transport";
      if(key!=prefix) return false;

      //------------------------------------------------------------------------
      // If spin-transport parameter is requested, then enable module
      //------------------------------------------------------------------------
      st::internal::enabled = true;

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
      test = "current-direction";
      if( word == test ){
         //-------------------------------------
         // Check for valid current directions
         //-------------------------------------
         test = "+x";
         if( value == test ){
          // set current along +x direction
          st::internal::current_direction = st::internal::px;
          return true;
         }
         //-------------------------------------
         test = "-x";
         if( value == test ){
          // set current along -x direction
          st::internal::current_direction = st::internal::mx;
          return true;
         }
         //-------------------------------------
         test = "+y";
         if( value == test ){
          // set current along +y direction
          st::internal::current_direction = st::internal::py;
          return true;
         }
         //-------------------------------------
         test = "-y";
         if( value == test ){
          // set current along -y direction
          st::internal::current_direction = st::internal::my;
          return true;
         }
         //-------------------------------------
         test = "+z";
         if( value == test ){
          // set current along +z direction
          st::internal::current_direction = st::internal::pz;
          return true;
         }
         //-------------------------------------
         test = "-z";
         if( value == test ){
          // set current along -z direction
          st::internal::current_direction = st::internal::mz;
          return true;
         }
         //--------------------------------------------
         // If here then no valid direction specified - error
         terminaltextcolor(RED);
         std::cerr << "Error - value for " << word << " on line " << line << "of input file must be be one of:" << std::endl;
         std::cerr << "\t+x" << std::endl;
         std::cerr << "\t-x" << std::endl;
         std::cerr << "\t+y" << std::endl;
         std::cerr << "\t-y" << std::endl;
         std::cerr << "\t+z" << std::endl;
         std::cerr << "\t-z" << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error - value for " << word << " on line " << line << "of input file must be be one of:" << std::endl;
         zlog << zTs() << "\t+x" << std::endl;
         zlog << zTs() << "\t-x" << std::endl;
         zlog << zTs() << "\t+y" << std::endl;
         zlog << zTs() << "\t-y" << std::endl;
         zlog << zTs() << "\t+z" << std::endl;
         zlog << zTs() << "\t-z" << std::endl;
         err::vexit();
      }
      //------------------------------------------------------------------------
      test = "applied-voltage";
      if( word == test ){
          // Set applied voltage in spin transport model
          double V = vin::str_to_double(value);
          vin::check_for_valid_value(V, word, line, prefix, unit, "potential", 0.0, 1000.0,"input","+0 - +1000 V");
          internal::voltage = V;
          return true;
      }
      //------------------------------------------------------------------------
      test = "environment-resistivity";
      if( word == test ){
          // Set resistivity for environment (cells with no atoms)
          double rho = vin::str_to_double(value);
          vin::check_for_valid_value(rho, word, line, prefix, unit, "resistivity", 1.0e-10, 1.0e12,"input","1E-10 - 1E12 Ohm metres");
          internal::environment_resistivity = rho;
          return true;
      }
      //------------------------------------------------------------------------
      test = "update-rate";
      if( word == test ){
          // Set resistivity for environment (cells with no atoms)
          uint64_t ur = vin::str_to_uint64(value);
          vin::check_for_valid_int(ur, word, line, prefix, 1,10000000,"input","1 - 1E7 time steps");
          st::internal::update_rate = ur;
          // set time counter to ensure initial calculation at start
          st::internal::time_counter = ur;
          return true;
      }
      // channel length
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

      // tunnel-barrier = true/false (if sp, then propagate across)
      // specific-heat-capacity? for joule heating

      // Check for material id > current array size and if so dynamically expand mp array
      if((unsigned int) super_index + 1 > internal::mp.size() && super_index + 1 < 101) internal::mp.resize(super_index + 1);

      //---------------------------------------------------------------------------
      // resistivity of continuous material made of this atom type [Parallel state, (Ohm m)]
      // rho := resistance per unit length and per unit of cross-sectional area
      std::string test = "resistivity";
      if( word == test ){
         // Set resistivity for atom type
         double rho = vin::str_to_double(value);
         vin::check_for_valid_value(rho, word, line, prefix, unit, "resistivity", 1.0e-20, 1.0e20,"input","1E-20 - 1E20 Ohm metres");
         st::internal::mp[super_index].resistivity.set(rho);
         return true;
      }
      //---------------------------------------------------------------------------
      // spin-resistivity of continuous material made of this atom type [Anti-parallel state, (Ohm m)]
      // tunnel barrier has a high R and high SR
      // normal materials SR is low ~ a few %
      test = "spin-resistivity";
      if( word == test ){
         // Set resistivity for atom type
         double rho = vin::str_to_double(value);
         vin::check_for_valid_value(rho, word, line, prefix, unit, "resistivity", 0.0, 1.0e20,"input","0.0 - 1E20 Ohm metres");
         st::internal::mp[super_index].spin_resistivity.set(rho);
         return true;
      }
      //---------------------------------------------------------------------------
      /*test = "tunnel-barrier"; // defines a tunnel barrier material where spin information is propogated through
      if( word == test ){
         // Set resistivity for atom type
         double rho = vin::str_to_double(value);
         vin::check_for_valid_bool(rho, word, line, prefix, unit, "resistivity", 1.0e-10, 1.0e12,"input","1E-10 - 1E12 Ohm metres");
         st::internal::mp[super_index].tunnel_barrier = true;
         return true;
      }*/
      //------------------------------------------------------------
      test  = "spin-transport-relaxation-torque";
      // aj parameter for material in slonczewski torque calculation
      if( word==test ){
         double aj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(aj, word, line, prefix, unit, "field", -1.0e-2, 1.0e2,"input","-100 - 100T");
         st::internal::mp[super_index].stt_rj.set(aj);
         return true;
      }
      //------------------------------------------------------------
      test = "spin-transport-precession-torque";
      // bj parameter for material in slonczewski torque calculation
      if( word==test ){
         double bj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(bj, word, line, prefix, unit, "field", -1.0e-2, 1.0e2,"input","-100 - 100T");
         st::internal::mp[super_index].stt_pj.set(bj);
         return true;
      }
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of spin_transport namespace
