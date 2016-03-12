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
#include "create.hpp"
#include "vio.hpp"

// Internal sim header
#include "internal.hpp"

namespace create{

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for create module
   //-----------------------------------------------------------------------------
   bool match_input_parameter(string const key, string const word, string const value, string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="create";
      if(key!=prefix) return false;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------

      /*std::string test="slonczewski-spin-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::DoublesFromString(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::slonczewski_spin_polarization_unit_vector = u;
         return true;
      }*/
      //--------------------------------------------------------------------
      // input parameter not found here
      return false;
   }

   //----------------------------------------------------------------------------------
   // material parameter match function
   //----------------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      // add prefix string
      std::string prefix="material:";

      // Check for material id > current array size and if so dynamically expand mp array
      if((unsigned int) super_index + 1 > create::internal::mp.size() && super_index + 1 < 101) create::internal::mp.resize(super_index + 1);

      //------------------------------------------------------------
      std::string test="alloy-host"; // determines host material
      if(word==test){
         // if this keyword is set, then atoms of this type will be scanned for alloy materials
         create::internal::mp[super_index].alloy_master=true;
         terminaltextcolor(YELLOW);
         std::cout << "Warning: Keyword \'alloy-host\' is deprecated and may be removed in a future release. Please use \'host-alloy\' instead." << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Warning: Keyword \'alloy-host\' is deprecated and may be removed in a future release. Please use \'host-alloy\' instead." << std::endl;
         return true;
      }
      //--------------------------------------------------------------------
      else
      test="host-alloy"; // determines host material
      if(word==test){
         // if this keyword is set, then atoms of this type will be scanned for alloy materials
         create::internal::mp[super_index].alloy_master=true;
         // check for type of host alloy
         test=""; // blank (assume homogeneous)
         if(value==test) create::internal::mp[super_index].host_alloy_type = internal::homogeneous;
         else
         test="homogeneous"; // default
         if(value==test) create::internal::mp[super_index].host_alloy_type = internal::homogeneous;
         else
         test="random"; // localised distribution
         if(value==test) create::internal::mp[super_index].host_alloy_type = internal::random;
         else
         test="granular"; // create distribution from intrinsic granular structure
         if(value==test) create::internal::mp[super_index].host_alloy_type = internal::granular;
         else
         test="checker-board"; // create distribution from intrinsic granular structure
         if(value==test) create::internal::mp[super_index].host_alloy_type = internal::checkerboard;
         // otherwise throw an error
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'material[" << super_index << "]:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"homogeneous\"" << std::endl;
            std::cerr << "\t\"random\"" << std::endl;
            std::cerr << "\t\"granular\"" << std::endl;
            std::cerr << "\t\"checker-board\"" << std::endl;
            zlog << zTs() << "Error - value for \'material[" << super_index << "]:" << word << "\' must be one of:" << std::endl;
            zlog << zTs() << "\t\"homogeneous\"" << std::endl;
            zlog << zTs() << "\t\"random\"" << std::endl;
            zlog << zTs() << "\t\"granular\"" << std::endl;
            zlog << zTs() << "\t\"checker-board\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }

         return true;
      }
      //--------------------------------------------------------------------
      else
      test="alloy-fraction"; // determines %mixing for disordered alloys
      if(word==test){
         double af=atof(value.c_str());
         vin::check_for_valid_value(af, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
         create::internal::mp[super_index].slave_material[sub_index].fraction=af;
         return true;
      }

/*material[1]:alloy-type[2] = native, reciprocal, homogeneous
material[1]:alloy-variance[2] = 0,1 pm xx%
material[1]:host-alloy = random, homogeneous, granular, checker-board
material[1]:host-alloy-smoothness = sharp, standard, smooth, 0-1
material[1]:host-alloy-scale = xx !nm
material[1]:save-alloy-profile (= file.dat)*/


      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of namespace create
