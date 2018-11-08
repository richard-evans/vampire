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

      // -----------------------------------------------------------------------------------
      // calculates the temperature dependant parallel component of the susceptability from the analytical equation
      // -----------------------------------------------------------------------------------

      std::vector<double> calculate_chi_para(int num_local_cells,
                                             std::vector<int>local_cell_array,
                                              int num_cells,
                                             double T){        //temperature

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

         std::vector<double> chi(num_cells,0.0);        //1D vector the store the parallel susceptability for each cell

         double PI= 3.14159;
         //double factor =  0.75947907;
         for (int i = 0; i < num_local_cells; i++){

            int cell = local_cell_array[i];

            //analytical expression for the parallel susceptability,
            if(T<Tc[cell]) chi[cell] =(a0/660.*Tc[cell])/(4.*PI)/(Tc[cell]-T)+a1*pow((Tc[cell]-T),1.)+ a2*pow((Tc[cell]-T),3.)+a3*pow((Tc[cell]-T),4.)+ a4*pow((Tc[cell]-T),6.)+ a5*pow((Tc[cell]-T),9.);
            else chi[cell] = (1.1*1.4/660.*Tc[cell])/(4*PI)/(T-Tc[cell]);

            //conversion from CGS to SI plus a small factor to stop the sum being less than 0
            chi[cell] = 1.0/(chi[cell]*9.54393845712027+0.308e-14 + 0.03);
            if (chi[cell] < 0) {
               std::cout << "error" << '\t' << Tc[cell] << '\t' << chi[cell] <<std::endl;
               std::cin.get();
            }
         }
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &chi[0],     num_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
         #endif

         return chi;            //returns the 1D vector for the susceptability,
      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
