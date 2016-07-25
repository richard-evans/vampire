//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "errors.hpp"
#include "gpu.hpp"
#include "vio.hpp"

namespace gpu{

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for ltmp settings
   //-----------------------------------------------------------------------------
   bool match_input_parameter(string const key, string const word, string const value, string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="gpu";
      if(key!=prefix) return false;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      std::string test="acceleration";
      if(word==test){
			gpu::acceleration = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="calculate-statistics-on-cpu";
      if(word==test){
			gpu::cpu_stats = true;
         return true;
      }
      //--------------------------------------------------------------------
      else{
         terminaltextcolor(RED);
         std::cerr << "Error - Unknown control statement \'"<< prefix << ":" << word << "\' on line " << line << " of input file" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      return false;
   }
} // end of namespace gpu
