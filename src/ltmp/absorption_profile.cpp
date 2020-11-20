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
#include "ltmp.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vmath.hpp"

// Local temperature pulse headers
#include "internal.hpp"

//--------------------------------------------------
// class functions for tabulated absorption profile
//

//--------------------------------------------------
// Constructor
//
ltmp::abs_t::abs_t():
  profile_loaded(false)
  {}

//--------------------------------------------------
// Function to return if profile is set
//
bool ltmp::abs_t::is_set(){

   if(profile_loaded) return true;
   else return false;

}

//--------------------------------------------------
// Add point to list of input points
//
void ltmp::abs_t::add_point(double height, double absorption){

   // check for valid height value
   if( height < 0.0 || height > 10000.0 ){
      std::cerr << "Error: height value " << height << " in absorption profile file is invalid. height values must be in the range 0.0 - 10000.0 A. Exiting." << std::endl;
      zlog << zTs() << "Error: height value " << height << " in absorption profile file is invalid. height values must be in the range 0.0 - 10000.0 A. Exiting." << std::endl;
      err::vexit();
   }

   // check for valid absorption value
   if( absorption < 0.0 || absorption > 1.0 ){
      std::cerr << "Error: absorption value " << height << " in absorption profile file is invalid. absorption values must be in the range 0.0 - 1.0. Exiting." << std::endl;
      zlog << zTs() << "Error: absorption value " << height << " in absorption profile file is invalid. absorption values must be in the range 0.0 - 1.0. Exiting." << std::endl;
      err::vexit();
   }

   // add values to list
   z.push_back(height);
   A.push_back(absorption);

   profile_loaded = true;

   return;
}

//--------------------------------------------------
// Creates a lookup table of interpolated functions
// to calculate absorption constant
//
void ltmp::abs_t::set_interpolation_table(){

   // Output informative message to log
   zlog << zTs() << "Determining interpolation variables for tabulated absorption profile." << std::endl;

   // Check for undefined absorption profile
   if(z.size()==0){
      z_max=10000.0;
      A_max=1.0;
      return;
   }

   // check z(i+1) > z(i)
   for(unsigned int i=1; i<z.size(); i++){
      if(z[i]<z[i-1]){
         std::cerr << "Error: height value "<< z[i] <<" on line " << i+1 << " is less than the previous value " << z[i-1] << " and must be given in ascending order. Exiting" << std::endl;
         zlog << zTs() << "Error: height value "<< z[i] <<" on line " << i+1 << " is less than the previous value " << z[i-1] << " and must be given in ascending order. Exiting" << std::endl;
         err::vexit();
      }
   }

   // loop over all value pairs and determine m and c values starting from i+1
   for(unsigned int i=1; i<z.size(); i++){
      m.push_back(vmath::interpolate_m(z[i-1],A[i-1],z[i],A[i]));
      c.push_back(vmath::interpolate_c(z[i-1],A[i-1],z[i],A[i]));
   }

   // determine maximum height specified in interpolation table
   z_max=z[z.size()-1];

   return;

}

//--------------------------------------------------
// Returns interpolated absorption constant for
// a given z-height
//
double ltmp::abs_t::get_absorption_constant(double height){

   // check for value larger than Tmax
   if(height>=z_max) return A_max;

   // loop over all values
   for(unsigned int i=1; i<z.size(); i++){
      double zmin = z[i-1];
      double zmax = z[i];
      // check if within range
      if(height >= zmin && height <= zmax){
         // calculate interpolated value
         return height*m[i]+c[i];
      }
   }

   return 0.0;

}
