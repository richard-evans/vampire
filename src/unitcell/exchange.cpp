//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <iostream>

// Vampire headers
#include "unitcell.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{
namespace internal{

//------------------------------------------------------------------------
// Function to determine exchange energy based on interaction range and
// exchange interaction function
//------------------------------------------------------------------------

double exponential_exchange(double range_sq){
   return exp(-sqrt(range_sq)/uc::internal::exchange_decay);
}

double exponential_exchange(double range_sq, double epA, double epB, double epC){
   return epA*exp(-epB*sqrt(range_sq)/uc::internal::exchange_decay) + epC;
}

double exchange(double range_sq, double nn_cutoff_sq, int i_mat, int j_mat){

   // Select program to run
   switch(uc::internal::exchange_function){

      case nearest_neighbour:{
         if(range_sq <= nn_cutoff_sq) return 1.0;
         else return 0.0;
         break;
      }

      case shell:{
         return 1.0;
      }

      case exponential:{
         return exponential_exchange(range_sq, exchange_parameter_A, exchange_parameter_B, exchange_parameter_C);
         break;
      }

      case NdFeB_exponential:{
         double A=36.9434;
         double B=1.25094;
         double C=-0.229572;
         const double Fe_ratio=0.69618016759*1.07692307692; // 560/520 = 1.07692307692
         const double J0Nd=Fe_ratio*4.06835e-20/16.0;
         // NdFe interaction
         if ((i_mat == 0 && j_mat == 1) || (i_mat == 1 && j_mat == 0) ){
            return 0.45*J0Nd;
         }
         // FeFe interaction
         if ((i_mat == j_mat) && i_mat == 1){
            return 2.0*2.179872e-21*(A*exp(-B*sqrt(range_sq))+C)*Fe_ratio;
         }
      }

      default:{
         if(range_sq <= nn_cutoff_sq) return 1.0;
         else return 0.0;
      }
   }

   return 0.0;

}

} // end of internal namespace
} // end of unitcell namespace
