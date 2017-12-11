//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
// standard library header files
#include <vector>

// vampire header files
#include "random.hpp"
#include "sim.hpp"

// Internal header file
#include "internal.hpp"

namespace montecarlo{

namespace internal{

// Function declarations
void mc_gaussian(const std::vector<double>&, std::vector<double>&);
void mc_spin_flip(const std::vector<double>&, std::vector<double>&);
void mc_uniform(std::vector<double>&);
void mc_angle(const std::vector<double>&, std::vector<double>&);
void mc_hinzke_nowak(const std::vector<double>&, std::vector<double>&);

///--------------------------------------------------------
///
///  Master function to call desired Monte Carlo move
///
///--------------------------------------------------------
void mc_move(const std::vector<double>& old_spin, std::vector<double>& new_spin){

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
void mc_angle(const std::vector<double>& old_spin, std::vector<double>& new_spin){

   new_spin[0]=old_spin[0]+mtrandom::gaussian()*internal::mc_delta_angle;
   new_spin[1]=old_spin[1]+mtrandom::gaussian()*internal::mc_delta_angle;
   new_spin[2]=old_spin[2]+mtrandom::gaussian()*internal::mc_delta_angle;

   // Calculate new spin length
   const double r = 1.0/sqrt (new_spin[0]*new_spin[0]+new_spin[1]*new_spin[1]+new_spin[2]*new_spin[2]);

   // Apply normalisation
   for (int i=0; i < new_spin.size(); i++) {
      new_spin[i]*=r;
   }

   return;

}

/// Spin flip move
/// Reverse spin direction
void mc_spin_flip(const std::vector<double>& old_spin, std::vector<double>& new_spin){

   new_spin[0]=-old_spin[0];
   new_spin[1]=-old_spin[1];
   new_spin[2]=-old_spin[2];

   return;

}

/// Random move
/// Place spin randomly on unit sphere
void mc_uniform(std::vector<double>& new_spin){

   new_spin[0]=mtrandom::gaussian();
   new_spin[1]=mtrandom::gaussian();
   new_spin[2]=mtrandom::gaussian();

   // Calculate new spin length
   const double r = 1.0/sqrt (new_spin[0]*new_spin[0]+new_spin[1]*new_spin[1]+new_spin[2]*new_spin[2]);

   // Apply normalisation
   for (int i=0; i < new_spin.size(); i++) {
      new_spin[i]*=r;
   }

   return;

}

/// Combination move selecting random move from spin_flip, angle and random
///
/// D. Hinzke, U. Nowak, Computer Physics Communications 121–122 (1999) 334–337
/// "Monte Carlo simulation of magnetization switching in a Heisenberg model for small ferromagnetic particles"
///
void mc_hinzke_nowak(const std::vector<double>& old_spin, std::vector<double>& new_spin){

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

} //end of namespace internal

} //end of namespace montecarlo
