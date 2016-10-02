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

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic{

   //------------------------------------ ---------------------------------------
   // Function to process input file parameters for micromagnetic module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="micromagnetic";
      if(key!=prefix) return false;



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
      std::string test="discretisation"; // whether the material is micromagnetic or atomistic
      if(word==test){

         // check for type of discretisation
         test="micromagnetic";
         if(value==test){
      //      discretisation_micromagnetic = true;
            // call micromagnetic function?? somehow? create a bool - micromagnetic ==true? but where is the bool stpred?
            return true;
         }
         test="atomistic"; // runs simualtion as normal - do i need to type anything to make it so that?
         if(value==test){
         //   discretisation_micromagnetic= false;
            return true;
         }
         test="multiscale"; // at the moment just runs a normal atomsitic simulation
         if(value==test){
      //      discretisation_micromagnetic = false;
            std::cerr << "Error - multiscale not yet implemented" << std::endl;
            return true;
         }
         // otherwise throw an error
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'material[" << super_index << "]:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"micromagnetic\"" << std::endl;
         std::cerr << "\t\"atomistic\"" << std::endl;
         std::cerr << "\t\"multiscale\"" << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error - value for \'material[" << super_index << "]:" << word << "\' must be one of:" << std::endl;
         zlog << zTs() << "\t\"micromagnetic\"" << std::endl;
         zlog << zTs() <<"\t\"atomistic\"" << std::endl;
         zlog << zTs() << "\t\"multiscale\"" << std::endl;
         err::vexit();

         return true;
      }

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of micromagnetic namespace
