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

      //calculates the temperature dependant parallel component of the susceptability from the analytical equation

      double calculate_chi_para(double T,int cell){

         // -----------------------------------------------------------------------------------
         //        chi_para = a0/4pi . Tc/(Tc-T) + sum_i ai (Tc -T)^i          if T < TC
         //        chi_para = b0 +1/4pi . Tc/(T-TC)                            if T > TC
         // -----------------------------------------------------------------------------------
         int shield = shield_number[cell];
         //Fit Parameter
         double a0 = 0.8;
         double a1 =-2.2e-07;
         double a2 = 1.95e-13;
         double a3 =-1.3e-17;
         double a4 =-4e-23;
         double a5 =-6.5076312364e-32;

         //cosntant to store the susepbability of the system

         double chi = 0;
         double PI= 3.14159;

         //calcualtes the susceptability for temerpature < Tc
         if(T<shield_Tc[shield]) chi =(a0/660.*shield_Tc[shield])/(4.*PI)/(shield_Tc[shield]-T)+a1*pow((shield_Tc[shield]-T),1.)+ a2*pow((shield_Tc[shield]-T),3.)+a3*pow((shield_Tc[shield]-T),4.)+ a4*pow((shield_Tc[shield]-T),6.)+ a5*pow((shield_Tc[shield]-T),9.);
         //calcualtes the susceptability for temperature greater than TC.
         else chi = (1.1*1.4/660.*shield_Tc[shield])/(4*PI)/(T-shield_Tc[shield]);

         //conversion from CGS to SI plus a small factor to stop the sum being less than 0
         chi = 1.0/(chi*9.54393845712027+0.308e-14 + 0.03);

         //returns the system suseptability
         return chi;
      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
