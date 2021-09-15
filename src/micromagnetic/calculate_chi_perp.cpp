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
#include "material.hpp"

namespace mm = micromagnetic::internal;

namespace micromagnetic{
namespace internal{

// Fit Parameters (local scoped constants to minimize function overhead)
const double a0 = 0.00211549427182711;
const double a1 = 0.110224660661792;
const double a2 = -0.855153260915204;
const double a3 = 3.42088365387997;
const double a4 = -7.84821585896818;
const double a5 = 10.3247035469514;
const double a6 = -6.85608273303224;
const double a7 = 0.797453198330591;
const double a8 = 1.53787854178089;
const double a9 = -0.627128148404525;
const double i4PI = 1.0/(4.0*3.14159265359);
const double i660 = 1.0/660.0;
const double chi0 = (a0+ a1+ a2+ a3 + a4 + a5 + a6 + a7+ a8+ a9); //normalisation when T=0

inline double chi_perp_fn(const double T, const double Tc){

   const double k = (1.068*Tc - T) / 1.068*Tc;
   const double k2 = k*k;
   const double k4 = k2*k2;
   const double k6 = k2*k4;
   const double k8 = k4*k4;

   if(T < 1.068*Tc) return chi0 / ( a0+ a1*k + a2*k2 + a3*k2*k + a4*k4 + a5*k4*k + a6*k6 + a7*k6*k + a8*k8 + a9*k8*k );
   else return chi0 / ( 0.8*1.4*i660*Tc *i4PI / (T-Tc) );

   //calculates the perpendicular susceptability from the analytical expression
   /*if(T[cell]<(1.068*Tc[cell])) chi_perp[cell] = a0+ a1*     ((1.068*Tc[cell]-T[cell]) /(1.068*Tc[cell]))+
                                                a2*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),2.)+
                                                a3*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),3.)+
                                                a4*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),4.) +
                                                a5*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),5.) +
                                                a6*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),6.) +
                                                a7*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),7.)+
                                                a8*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),8.) +
                                                a9*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),9.);
   else chi_perp[cell] = (0.8*1.4/660.*Tc[cell])/(4*PI)/(T[cell]-Tc[cell]);*/

//chi_0 / (chi_perp[cell]);

}

// -----------------------------------------------------------------------------------
// calculates the temperature dependant perpendicular component of the susceptability from the analytical equation
// -----------------------------------------------------------------------------------
void calculate_chi_perp(int num_mm_cells,
                       std::vector<int>& list_of_mm_cells,
                       std::vector<double>& chi_perp,
                       std::vector<double>& T,
                       std::vector<double>& Tc){

   //Fit Parameters
   const double a0 = 0.00211549427182711;
   const double a1 = 0.110224660661792;
   const double a2 = -0.855153260915204;
   const double a3 = 3.42088365387997;
   const double a4 = -7.84821585896818;
   const double a5 = 10.3247035469514;
   const double a6 = -6.85608273303224;
   const double a7 = 0.797453198330591;
   const double a8 = 1.53787854178089;
   const double a9 = -0.627128148404525;
   const double PI= 3.14159;

   const double chi_0 = (a0+ a1+ a2+ a3 + a4 + a5 + a6 + a7+ a8+ a9); //normalisation when T=0


   // compute chi value for cells (on this processor)
   for (int lc = 0; lc < num_mm_cells; lc++){

      int cell = list_of_mm_cells[lc];

      //chi_perp[cell] = chi_perp_fn(T[cell], Tc[cell]);

      //calculates the perpendicular susceptability from the analytical expression
      if(T[cell]<(1.068*Tc[cell])) chi_perp[cell] = a0+ a1*     ((1.068*Tc[cell]-T[cell]) /(1.068*Tc[cell]))+
                                                   a2*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),2.)+
                                                   a3*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),3.)+
                                                   a4*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),4.) +
                                                   a5*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),5.) +
                                                   a6*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),6.) +
                                                   a7*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),7.)+
                                                   a8*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),8.) +
                                                   a9*pow((((1.068*Tc[cell])-T[cell])/(1.068*Tc[cell])),9.);
      else chi_perp[cell] = (0.8*1.4/660.*Tc[cell])/(4*PI)/(T[cell]-Tc[cell]);

      //normalises
      //returns one over perpendicularsusceptability
      //because its always used as 1/chi so we save it as 1/chi to save computaiton.
      //divided by ms because in the function we use unit vectors
      chi_perp[cell] = -chi_0 / (ms[cell] * chi_perp[cell]);
      //chi_perp[cell] = chi_0 / (chi_perp[cell]); - should actually be this one!

   }

   return;

}

} //closes the internal namspace
}  //closes the micromagnetic namespace
