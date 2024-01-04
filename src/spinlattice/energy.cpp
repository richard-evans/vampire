//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//   implementation based on the paper Phys. Rev. B 103, 024429, (2021) M.Strungaru, M.O.A. Ellis et al
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "sld.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include "atoms.hpp"
#include "create.hpp"
#include "constants.hpp"
#include "material.hpp"






// sld module headers
#include "internal.hpp"


namespace sld{

   double compute_potential_energy(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& type_array){

//
                double pot=0.0;
                for (int at=start_index;at<end_index;at++){
                    pot += sld::internal::potential_eng[at];
                }

                pot /= (end_index-start_index);
            //   std::cout<<"spin temp param "<<mp::material[0].mu_s_SI /constants::kB <<"\t"<<mp::material[0].mu_s_SI<<"\t"<<constants::kB<<std::endl;
      return pot;

   }//end of potential energy
//
double compute_kinetic_energy(const int start_index, // first atom for exchange interactions to be calculated
            const int end_index,
            const std::vector<int>& type_array, // type for atom
            std::vector<double>& velo_array_x, // coord vectors for atoms
            std::vector<double>& velo_array_y,
            std::vector<double>& velo_array_z){

//
            double kinetic=0;
             for (int at=start_index;at<end_index;at++){
                 double vx = atoms::x_velo_array[at];
                 double vy = atoms::y_velo_array[at];
                 double vz = atoms::z_velo_array[at];
                 const unsigned int imat = atoms::type_array[at];
                
                 kinetic+=sld::internal::mp[imat].mass.get()*(vx*vx+vy*vy+vz*vz);
             }

            kinetic*=0.5/(end_index-start_index);
           // std::cout<<"kinetic en "<<kinetic<<std::endl;

   return kinetic;

}//end of kinetic energy

double compute_effective_J(const int start_index, // first atom for exchange interactions to be calculated
            const int end_index,
            std::vector<double>& sumJ){

             double sum1J=0.0;
             for (int at=start_index;at<end_index;at++){
                 sum1J += sumJ[at];
             }
             sum1J /= (end_index-start_index);
   return sum1J;

}//end of compute eff exchange

double compute_effective_C(const int start_index, // first atom for exchange interactions to be calculated
            const int end_index,
            std::vector<double>& sumC){

             double sum1C=0.0;
             for (int at=start_index;at<end_index;at++){
                 sum1C += sumC[at];
             }
             sum1C /= (end_index-start_index);
   return sum1C;

}//end of compute eff coupling

double compute_coupling_energy(const int start_index, // first atom for exchange interactions to be calculated
            const int end_index){


             double sumC=0.0;
             for (int at=start_index;at<end_index;at++){
                 sumC += sld::internal::coupl_eng[at];
             }

             sumC /= (end_index-start_index);
             sumC *= mp::material[0].mu_s_SI/1.602176634e-19; //in eV at the moment

   return sumC;

}//end of potential energy

double compute_exchange_energy(const int start_index, // first atom for exchange interactions to be calculated
            const int end_index){

//
             double sumJ=0.0;
             for (int at=start_index;at<end_index;at++){
                 sumJ += sld::internal::exch_eng[at];
             }

             sumJ /= (end_index-start_index);
             sumJ *= mp::material[0].mu_s_SI/1.602176634e-19;//in eV at the moment
   return sumJ;

}//end of potential energy


   } // end of sld namespace
