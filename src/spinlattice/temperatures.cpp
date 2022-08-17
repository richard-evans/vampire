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
#include "constants.hpp"
#include "material.hpp"

// sld module headers
#include "internal.hpp"


namespace sld{


   double compute_spin_temperature(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& type_array, // type for atom
               std::vector<double>& x_spin_array, // coord vectors for atoms
               std::vector<double>& y_spin_array,
               std::vector<double>& z_spin_array,
               std::vector<double>& fields_array_x, //  vectors for fields
               std::vector<double>& fields_array_y,
               std::vector<double>& fields_array_z){

                double SxH2=0.0;
                double SH=0.0;
                for (int at=start_index;at<end_index;at++){
                    double Sx = x_spin_array[at];
                    double Sy = y_spin_array[at];
                    double Sz = z_spin_array[at];
                    double Hx = fields_array_x[at];
                    double Hy = fields_array_y[at];
                    double Hz = fields_array_z[at];

                     double SxHx = Sy * Hz - Sz * Hy;
                     double SxHy = Sz * Hx - Sx * Hz;
                     double SxHz = Sx * Hy - Sy * Hx;
                     SxH2  = SxH2+ SxHx*SxHx + SxHy*SxHy + SxHz*SxHz;
                     SH  = SH +  Sx * Hx + Sy * Hy + Sz*Hz;

                }

               double T_spin=0.5* mp::material[0].mu_s_SI /constants::kB * SxH2 / SH;

      return T_spin;

      }//end of spin_temperature


double compute_lattice_temperature(const int start_index, // first atom for exchange interactions to be calculated
            const int end_index,
            const std::vector<int>& type_array, // type for atom
            std::vector<double>& velo_array_x, // coord vectors for atoms
            std::vector<double>& velo_array_y,
            std::vector<double>& velo_array_z){

            double kinetic=0;
             for (int at=start_index;at<end_index;at++){
                 double vx = velo_array_x[at];
                 double vy = velo_array_y[at];
                 double vz = velo_array_z[at];
                 kinetic+=vx*vx+vy*vy+vz*vz;
             }

            kinetic*=sld::internal::mp[0].mass.get()*0.5/(end_index-start_index);
            double T_lat=( 2.0 /(3.0*constants::kB_eV))*kinetic;

   return T_lat;

}//end of lattice_temperature

   } // end of sld namespace
