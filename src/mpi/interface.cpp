//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

namespace vmpi{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for mpi module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="sim";
      if(key!=prefix) return false;

      std::string test="";

      //--------------------------------------------------------------------
      test="mpi-mode";
      if(word==test){
         test="geometric-decomposition";
         if(value==test){
            vmpi::mpi_mode=0;
            return true;
         }
         test="replicated-data";
         if(value==test){
            vmpi::mpi_mode=1;
            vmpi::replicated_data_staged=false;
            return true;
         }
         test="replicated-data-staged";
         if(value==test){
            vmpi::mpi_mode=1;
            vmpi::replicated_data_staged=true;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"geometric-decomposition\"" << std::endl;
            std::cerr << "\t\"replicated-data\"" << std::endl;
            std::cerr << "\t\"replicated-data-staged\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //--------------------------------------------------------------------
      test="mpi-ppn";
      if(word==test){
         int ppn=atoi(value.c_str());
         vin::check_for_valid_int(ppn, word, line, prefix, 1, 1024,"input","1 - 1024");
         vmpi::ppn=ppn;
         return true;
      }
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of vmpi namespace
