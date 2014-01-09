//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2013 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify 
//  it under the terms of the GNU General Public License as published by 
//  the Free Software Foundation; either version 2 of the License, or 
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but 
//  WITHOUT ANY WARRANTY; without even the implied warranty of 
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License 
//  along with this program; if not, write to the Free Software Foundation, 
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
//-------------------------------------------------------------------
//
//    Implementation of various algorithms for Monte Carlo moves
//
//    Takes initial spin as argument and returns new spin positions
//
//    (c) R F L Evans 2013 University of York
//
//-------------------------------------------------------------------
// standard library header files
#include <valarray>

// vampire header files
#include "random.hpp"
#include "sim.hpp"

namespace sim{

// Function declarations
std::valarray<double> mc_gaussian(std::valarray<double>&);
std::valarray<double> mc_spin_flip(std::valarray<double>&);
std::valarray<double> mc_uniform(std::valarray<double>&);
std::valarray<double> mc_angle(std::valarray<double>&);
std::valarray<double> mc_hinzke_nowak(std::valarray<double>&);

///--------------------------------------------------------
///
///  Master function to call desired Monte Carlo move
///
///--------------------------------------------------------
std::valarray<double> mc_move(std::valarray<double>& old_spin){
   
   // Reference enum list for readability
   using namespace sim;
   //enum mc_algorithms { spin_flip, uniform, angle, hinzke_nowak};

   // Select algorithm using case statement
   switch(sim::mc_algorithm){
      
      case spin_flip:
         return mc_spin_flip(old_spin);
         break;
      case uniform:
         return mc_uniform(old_spin);
         break;
      case angle:
         return mc_angle(old_spin);
         break;
      case hinzke_nowak:
         return mc_hinzke_nowak(old_spin);
         break;
      default:
         return mc_hinzke_nowak(old_spin);
         break;
   }
}

/// Angle move
/// Move spin within cone near old position
std::valarray<double> mc_angle(std::valarray<double>& old_spin){
   
   // Declare new spin
   std::valarray<double> new_spin(3);

   new_spin[0]=old_spin[0]+mtrandom::gaussian()*sim::mc_delta_angle;
   new_spin[1]=old_spin[1]+mtrandom::gaussian()*sim::mc_delta_angle;
   new_spin[2]=old_spin[2]+mtrandom::gaussian()*sim::mc_delta_angle;

   // Calculate new spin length
   const double r = 1.0/sqrt (new_spin[0]*new_spin[0]+new_spin[1]*new_spin[1]+new_spin[2]*new_spin[2]); 
   
   // Apply normalisation
   new_spin*=r;

   return new_spin;

}

/// Spin flip move
/// Reverse spin direction
std::valarray<double> mc_spin_flip(std::valarray<double>& old_spin){
   
   // Declare new spin
   std::valarray<double> new_spin(3);

   new_spin[0]=-old_spin[0];
   new_spin[1]=-old_spin[1];
   new_spin[2]=-old_spin[2];

   return new_spin;  

}

/// Random move
/// Place spin randomly on unit sphere
std::valarray<double> mc_uniform(std::valarray<double>& old_spin){
   
  // Declare new spin
   std::valarray<double> new_spin(3);

   new_spin[0]=mtrandom::gaussian();
   new_spin[1]=mtrandom::gaussian();
   new_spin[2]=mtrandom::gaussian();

   // Calculate new spin length
   const double r = 1.0/sqrt (new_spin[0]*new_spin[0]+new_spin[1]*new_spin[1]+new_spin[2]*new_spin[2]); 

   // Apply normalisation
   new_spin*=r;

   return new_spin;

}

/// Combination move selecting random move from spin_flip, angle and random
///
/// D. Hinzke, U. Nowak, Computer Physics Communications 121–122 (1999) 334–337
/// "Monte Carlo simulation of magnetization switching in a Heisenberg model for small ferromagnetic particles"
/// 
std::valarray<double> mc_hinzke_nowak(std::valarray<double>& old_spin){
   
   // Select random move type
   const int pick_move=int(3.0*mtrandom::grnd());
   
      switch(pick_move){
      
         case 0:
            return mc_spin_flip(old_spin);
            break;
         case 1:
            return mc_uniform(old_spin);
            break;
         case 2:
            return mc_angle(old_spin);
            break;
         default:
            return mc_angle(old_spin);
            break;
      }
      
}

}