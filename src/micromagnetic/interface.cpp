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

   //---------------------------------------------------------------------------
   // Function to process input file parameters for micromagnetic module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="micromagnetic";
      if(key!=prefix) return false;


      std::string test="integrator";
      if(word==test){
          test="llg";
          if(value==test){
              micromagnetic::integrator=0;
              return true;
          }
          test="llb";
          if(value==test){
              micromagnetic::integrator=1;
              return true;
          }
          else{
          terminaltextcolor(RED);
              std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
              std::cerr << "\t\"llg\"" << std::endl;
              std::cerr << "\t\"llb\"" << std::endl;
          terminaltextcolor(WHITE);
              err::vexit();
          }
      }
      test="atomistic-steps-per-micromagnetic-step";
            std::cout << word << "\t" << test << std::endl;

      if(word==test){
        std::cout << "a" <<std::endl;
        double dt=atof(value.c_str());
        //vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
        micromagnetic::num_atomic_steps_mm =dt;
        return true;
      }


      test="discretisation"; // whether the material is micromagnetic or atomistic
      if(word==test){

         // check for type of discretisation
         test="micromagnetic";
         if(value==test){
            discretisation_type = 1;
            return true;
         }
         test="atomistic"; // runs simualtion as normal (atomistics)
         if(value==test){
            discretisation_type = 0;
            return true;
         }
         test="multiscale"; // at the moment just runs a normal atomsitic simulation
         if(value==test){
            discretisation_type = 2;
            return true;
         }
         // otherwise throw an error
         terminaltextcolor(RED);
         std::cerr << "Error - value for" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"micromagnetic\"" << std::endl;
         std::cerr << "\t\"atomistic\"" << std::endl;
         std::cerr << "\t\"multiscale\"" << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error - value for" << word << "\' must be one of:" << std::endl;
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

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      // add prefix string
      std::string prefix="material:";


      //--------------------------------------------------------------------
      std::string test="micromagnetic-discretisation-enabled";
            std::cout << word << "\t" << test << std::endl;
      if(word==test){
        std::cout << "a" << std::endl;
         double K=atof(value.c_str());
         vin::read_material[super_index].micromagnetic_enabled=true;
          std::cout << super_index << '\t' << vin::read_material[super_index].micromagnetic_enabled << std::endl;
            return true;
      }


      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of micromagnetic namespace
