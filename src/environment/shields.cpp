//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
// Vampire headers
#include "environment.hpp"

// micromagnetic module headers
#include "internal.hpp"

// C++ headers
#include <math.h>

namespace environment{

namespace internal{

bool  in_x(double x, double z){

   const double xr = x - 500.0; // shift to zero
   const double zr = z - 300.0; // shift to base of side shield

   const double zmin = 200.0 * exp(-((fabs(xr)-190)*0.01));
   std::cout << x << '\t' << z << '\t'  << xr << '\t' << zr << '\t' << zmin << std::endl;
   if(zr > zmin && z < 500.0) return true;
   return false;

}

//------------------------------------------------------------------------------
// Function to calculate basic shield geometry for reader
//------------------------------------------------------------------------------
bool in_shield(double x, double y, double z){

   // height of inner sensor region
   const double stack_height = 200; // Angstroms

   const double xr = x;
   const double yr = y;
   const double zr = z; // reduced height

   // Bottom shield
   if(zr < 300.0 && zr > 0.0) return true;

   // Top shield (++z) 31-51 nm
   if(zr > 520.0 && zr < 720.0) return true;

   // Top shield (+z) 52-72 nm
   if(z > 740.0 && z < 940.0) return true;

   // side shields
   if(in_x(xr, zr)) return true;

   return false;

}

}

}
