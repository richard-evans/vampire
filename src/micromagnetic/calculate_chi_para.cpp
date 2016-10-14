// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include <math.h>
namespace micromagnetic
{

   namespace internal
   {
      std::vector<double> calculate_chi_para(int num_cells, double T)
      {
         //Fit Parameter
         double a0 = 0.8;
         double a1 =-2.2e-07;
         double a2 = 1.95e-13;
         double a3 =-1.3e-17;
         double a4 =-4e-23;
         double a5 =-6.5076312364e-32;
      //   double chi_0 = (a0+ a1+ a2+ a3 + a4 + a5)*9.54393845712027+0.308e-14 + 0.0011;
         std::vector<double> chi_CGS(num_cells,0.0);
         //double chi_SI  = 0.0;
         std::vector<double> chi(num_cells,0.0);
         double PI= 3.14159;
         //double factor =  0.75947907;
         for (int cell = 0; cell < num_cells; cell++)
         {
            if(T<Tc[cell]) chi_CGS[cell] =(a0/660.*Tc[cell])/(4.*PI)/(Tc[cell]-T)+a1*pow((Tc[cell]-T),1.)+ a2*pow((Tc[cell]-T),3.)+a3*pow((Tc[cell]-T),4.)+ a4*pow((Tc[cell]-T),6.)+ a5*pow((Tc[cell]-T),9.);
            else chi_CGS[cell] = (1.1*1.4/660.*Tc[cell])/(4*PI)/(T-Tc[cell]);

            //chi_SI = 4*PI*chi_CGS;     // CHI_SI = 4*PI Chi_CGS
            //chi = chi_SI*4*A_FePt*A_FePt*A_FePt/MU_S_FePt/MU_0;
            chi[cell] = chi_CGS[cell]*9.54393845712027+0.308e-14 + 0.0011; // (Tesla)
            chi[cell] = chi[cell];
         }
         return chi; // [T]
      }
   }
}
