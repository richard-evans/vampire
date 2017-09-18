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
#include "environment.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include <math.h>
namespace environment{

   namespace internal{
      //calcuaktes the temperature dependant perpendicular component of the susceptability from the analytical equation

      double calculate_chi_perp(double T){

         // -----------------------------------------------------------------------------------
         //               chi_perp = a0 + sum_i ai . (Tc/1.068Tc)^i        if T < TC
         //               chi_perp = b0 +1/4pi . Tc/(T-TC)                 if T > TC
         // -----------------------------------------------------------------------------------

         double chi = 0;     //1D vector the store the perpendicular susceptability for each cell

         //Fit Parameters
         double a0 = 0.00211549427182711;
         double a1 = 0.110224660661792;
         double a2 = -0.855153260915204;
         double a3 = 3.42088365387997;
         double a4 = -7.84821585896818;
         double a5 = 10.3247035469514;
         double a6 = -6.85608273303224;
         double a7 = 0.797453198330591;
         double a8 = 1.53787854178089;
         double a9 = -0.627128148404525;
         double PI= 3.14159;
         //normalisation when T=0
         double chi_0 = (a0+ a1+ a2+ a3 + a4 + a5 + a6 + a7+ a8+ a9);
         //conversion to SI
         chi_0 = chi_0*9.54393845712027;

         //calculates the perpendicular susceptability from the analytical expression for T<Tc
         if(T<(1.068*Tc)) chi = a0+ a1*((1.068*Tc-T)/(1.068*Tc))+ a2*pow((((1.068*Tc)-T)/(1.068*Tc)),2.)+ a3*pow((((1.068*Tc)-T)/(1.068*Tc)),3.)+ a4*pow((((1.068*Tc)-T)/(1.068*Tc)),4.) + a5*pow((((1.068*Tc)-T)/(1.068*Tc)),5.) + a6*pow((((1.068*Tc)-T)/(1.068*Tc)),6.) + a7*pow((((1.068*Tc)-T)/(1.068*Tc)),7.)+ a8*pow((((1.068*Tc)-T)/(1.068*Tc)),8.) + a9*pow((((1.068*Tc)-T)/(1.068*Tc)),9.);
         //calcualtes the susceptability for T > TC
         else chi = (0.8*1.4/660.*Tc)/(4*PI)/(T-Tc);
         //converts to SI
         chi = (chi*9.54393845712027); // (Tesla)
         //normalises with respect to the 0K susceptability.
         chi = chi/chi_0;

         //returns the anisotropy constant - the susceptability is only used for this
         chi = -ku/Ms*chi;

         return chi;            //returns the 1D vector containing the perpencicular susceptability for each cell,
      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
