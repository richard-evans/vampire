//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers
#include "environment.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for environment module
//--------------------------------------------------------------------------------
namespace environment{


     extern bool enabled;
     extern int demag_update_rate;
     extern std::vector < double > environment_field_x;
     extern std::vector < double > environment_field_y;
     extern std::vector < double > environment_field_z;


     extern std::vector < double > atomic_field_x;
     extern std::vector < double > atomic_field_y;
     extern std::vector < double > atomic_field_z;


     int LLB(double temperature,
             double Hx,
             double Hy,
             double Hz,
             double H,
             double dt);
   //-----------------------------------------------------------------------------
   // Function to initialise environment module
   //-----------------------------------------------------------------------------
   void initialize(   double system_dimensions_x,
                      double system_dimensions_y,
                      double system_dimensions_z);

   //---------------------------------------------------------------------------
   // Function to process input file parameters for environment module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of environment namespace

#endif //ENVIRONMENT_H_
