// Vampire headers
#include "environment.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include <math.h>
namespace environment{

   namespace internal{


     //calcuaktes the temperature dependant parallel component of the susceptability from the analytical equation

      double calculate_chi_para(double T){

         // -----------------------------------------------------------------------------------

         //        chi_para = a0/4pi . Tc/(Tc-T) + sum_i ai (Tc -T)^i          if T < TC
         //        chi_para = b0 +1/4pi . Tc/(T-TC)                            if T > TC

         // -----------------------------------------------------------------------------------


         //Fit Parameter
         double a0 = 0.8;
         double a1 =-2.2e-07;
         double a2 = 1.95e-13;
         double a3 =-1.3e-17;
         double a4 =-4e-23;
         double a5 =-6.5076312364e-32;

         double chi = 0.0;        //1D vector the store the parallel susceptability for each cell

         double PI= 3.14159;


        if(T<Tc) chi =(a0/660.*Tc)/(4.*PI)/(Tc-T)+a1*pow((Tc-T),1.)+ a2*pow((Tc-T),3.)+a3*pow((Tc-T),4.)+ a4*pow((Tc-T),6.)+ a5*pow((Tc-T),9.);
        else chi = (1.1*1.4/660.*Tc)/(4*PI)/(T-Tc);

        //conversion from CGS to SI plus a small factor to stop the sum being less than 0
        chi = 1.0/(chi*9.54393845712027+0.308e-14 + 0.03);

         return chi;            //returns the 1D vector for the susceptability,
      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
