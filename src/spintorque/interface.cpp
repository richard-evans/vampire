//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <fstream>
#include <iostream>

// Vampire headers
#include "spintorque.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{

int match_material(string const word, string const value, string const unit, int const line, int const super_index){

   // add prefix string
   std::string prefix="material:";

   // Check for material id > current array size and if so dynamically expand mp array
   if((unsigned int) super_index + 1 > st::internal::mp.size() && super_index + 1 < 101) st::internal::mp.resize(super_index + 1);

   //------------------------------------------------------------
   std::string test="spin-diffusion-length";
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double lsdl=atof(value.c_str());
      //std::cout << "lsdl" << super_index << "\t" << lsdl << "\t" << value.c_str() << std::endl;
      vin::check_for_valid_value(lsdl, word, line, prefix, unit, "length", 0.01, 1.0e10,"material"," 0.01 - 1e10 Angstroms");
      st::internal::mp[super_index].lambda_sdl=lsdl*1.e-10; // defined in metres
      std::cout << "lsdl" << super_index << "\t" << lsdl << std::endl;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="spin-polarisation-conductivity";
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double betac=atof(value.c_str());
      vin::check_for_valid_value(betac, word, line, prefix, unit, "none", 1.0e-3, 1.0e3,"material"," 0.001 - 1000");
      st::internal::mp[super_index].beta_cond=betac;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="spin-polarisation-diffusion";
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double betad=atof(value.c_str());
      vin::check_for_valid_value(betad, word, line, prefix, unit, "none", 1.0e-3, 1.0e3,"material"," 0.001 - 1000");
      st::internal::mp[super_index].beta_diff=betad;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="spin-accumulation";
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double sa=atof(value.c_str());
      vin::check_for_valid_value(sa, word, line, prefix, unit, "none", 1.0e-6, 1.0e10,"material"," 1.0e-6 - 1.0e10");
      st::internal::mp[super_index].sa_infinity=sa;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="diffusion-constant";
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double dc=atof(value.c_str());
      vin::check_for_valid_value(dc, word, line, prefix, unit, "none", 1.0e-9, 100.0,"material"," 1.0e-9 - 100"); //Angstroms^2/s
      st::internal::mp[super_index].diffusion=dc;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="sd-exchange-constant";
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double sd=atof(value.c_str());
      vin::check_for_valid_value(sd, word, line, prefix, unit, "energy", 1.0e-30, 1.e-17,"material"," 1.0e-30 - 1.e-21 J");
      st::internal::mp[super_index].sd_exchange=sd;
      return EXIT_SUCCESS;
   }

   //--------------------------------------------------------------------
   // keyword not found
   //--------------------------------------------------------------------
   return EXIT_FAILURE;

}


//spin-torque:current-direction = x/y/z
//spin-torque:micro-cell-size
//spin-torque:micro-cell-thickness
//spin-torque:current-density*/


} // end of namespace st
