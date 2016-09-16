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
void mc_gaussian(const std::valarray<double>&, std::valarray<double>&);
void mc_spin_flip(const std::valarray<double>&, std::valarray<double>&);
void mc_uniform(std::valarray<double>&);
void mc_angle(const std::valarray<double>&, std::valarray<double>&);
void mc_hinzke_nowak(const std::valarray<double>&, std::valarray<double>&);

///--------------------------------------------------------
///
///  Master function to call desired Monte Carlo move
///
///--------------------------------------------------------
void mc_move(const std::valarray<double>& old_spin, std::valarray<double>& new_spin){

   // Reference enum list for readability
   using namespace sim;
   //enum mc_algorithms { spin_flip, uniform, angle, hinzke_nowak};

   // Select algorithm using case statement
   switch(sim::mc_algorithm){
      
      case spin_flip:
         mc_spin_flip(old_spin, new_spin);
         break;
      case uniform:
         mc_uniform(new_spin);
         break;
      case angle:
         mc_angle(old_spin, new_spin);
         break;
      case hinzke_nowak:
         mc_hinzke_nowak(old_spin, new_spin);
         break;
      default:
         mc_hinzke_nowak(old_spin, new_spin);
         break;
   }
   return;
}

/// Angle move
/// Move spin within cone near old position
void mc_angle(const std::valarray<double>& old_spin, std::valarray<double>& new_spin){

   new_spin[0]=old_spin[0]+mtrandom::gaussian()*sim::mc_delta_angle;
   new_spin[1]=old_spin[1]+mtrandom::gaussian()*sim::mc_delta_angle;
   new_spin[2]=old_spin[2]+mtrandom::gaussian()*sim::mc_delta_angle;

   // Calculate new spin length
   const double r = 1.0/sqrt (new_spin[0]*new_spin[0]+new_spin[1]*new_spin[1]+new_spin[2]*new_spin[2]); 

   // Apply normalisation
   new_spin*=r;

   return;

}

/// Spin flip move
/// Reverse spin direction
void mc_spin_flip(const std::valarray<double>& old_spin, std::valarray<double>& new_spin){

   new_spin[0]=-old_spin[0];
   new_spin[1]=-old_spin[1];
   new_spin[2]=-old_spin[2];

   return;

}

/// Random move
/// Place spin randomly on unit sphere
void mc_uniform(std::valarray<double>& new_spin){

   new_spin[0]=mtrandom::gaussian();
   new_spin[1]=mtrandom::gaussian();
   new_spin[2]=mtrandom::gaussian();

   // Calculate new spin length
   const double r = 1.0/sqrt (new_spin[0]*new_spin[0]+new_spin[1]*new_spin[1]+new_spin[2]*new_spin[2]); 

   // Apply normalisation
   new_spin*=r;

   return;

}

/// Combination move selecting random move from spin_flip, angle and random
///
/// D. Hinzke, U. Nowak, Computer Physics Communications 121–122 (1999) 334–337
/// "Monte Carlo simulation of magnetization switching in a Heisenberg model for small ferromagnetic particles"
/// 
void mc_hinzke_nowak(const std::valarray<double>& old_spin, std::valarray<double>& new_spin){

   // Select random move type
   const int pick_move=int(3.0*mtrandom::grnd());

      switch(pick_move){
         case 0:
            mc_spin_flip(old_spin, new_spin);
            break;
         case 1:
            mc_uniform(new_spin);
            break;
         case 2:
            mc_angle(old_spin, new_spin);
            break;
         default:
            mc_angle(old_spin, new_spin);
            break;
      }
      return;
}

}

