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

// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include <math.h>

namespace micromagnetic{
namespace internal{

// Fit Parameters (local scoped constants to minimize function overhead)
const double a0 = 0.8;
const double a1 =-2.2e-07;
const double a2 = 1.95e-13;
const double a3 =-1.3e-17;
const double a4 =-4e-23;
const double a5 =-6.5076312364e-32;
const double i4PI = 1.0/(4.0*3.14159265359);
const double i660 = 1.0/660.0;

// -----------------------------------------------------------------------------------
//        chi_para = a0/4pi . Tc/(Tc-T) + sum_i ai (Tc -T)^i          if T < TC
//        chi_para = b0 +1/4pi . Tc/(T-TC)                            if T > TC
// -----------------------------------------------------------------------------------
inline double chi_para_fn(const double T, const double Tc){

   const double TcmT = Tc - T;
   const double TcmT3 = TcmT*TcmT*TcmT;
   const double TcmT6 = TcmT3*TcmT3;
   const double TcmT9 = TcmT3*TcmT6;

   if(T<Tc) return 1.0 / ( 9.54393845712027*( (a0*Tc*i660*i4PI)/(TcmT) +
                                              a1*TcmT + a2*TcmT3 +
                                              a3*TcmT*TcmT*TcmT*TcmT +
                                              a4*TcmT6 + a5*TcmT9 ) );
   else return 1.0 / ( 9.54393845712027*( 1.1*1.4*Tc*i660*i4PI/(T-Tc) ) );

}

// -----------------------------------------------------------------------------------
// calculates the temperature dependant parallel component of the susceptability from the analytical equation
// -----------------------------------------------------------------------------------
void calculate_chi_para(int num_mm_cells,
                        std::vector<int>& list_of_mm_cells,
                        std::vector<double>& chi_para,
                        std::vector<double>& T,
                        std::vector<double>& Tc){

   // compute chi value for cells (on this processor)
   for (int lc = 0; lc < num_mm_cells; lc++){
      int cell = list_of_mm_cells[lc];
      chi_para[cell] = chi_para_fn(T[cell], Tc[cell]);
   }

   #ifdef MPICF
      //MPI_Allreduce(MPI_IN_PLACE, &chi[0],     num_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
   #endif

   return;            //returns the 1D vector for the susceptability,
}

} //closes the internal namspace
}  //closes the micromagnetic namespace
