//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sergiu Ruta 2022. All rights reserved.
//
//   Email: sergiu.ruta@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "spinwaves.hpp"
#include "vmpi.hpp"
// sw module headers
#include "internal.hpp"
#include "iostream"
#include <iomanip>
#include <sstream>

//sergiu for SW
#include "unitcell.hpp"
#include "vector"
#include <cmath>
#include "fstream"
#include "atoms.hpp"


namespace spinwaves {

   // Local variables
   double sx;
   int atom;
   int spec;
   int index;
   int reduc_index;

   //----------------------------------------------------------------------------
   // Function to initialize sw module
   //----------------------------------------------------------------------------
   void fft_in_space( const std::vector<double>& rx,
                   const std::vector<double>& ry,
                   const std::vector<double>& rz,
                   const int time ){

      for(int k=0;k<internal::nk;k++){

         #ifdef MPICF

            //reset values
            for (int spec = 0; spec < internal::nspec; spec++){
               skx_r[k*internal::nspec+spec] = 0.0;
               skx_i[k*internal::nspec+spec] = 0.0;
            }

            // loop over atoms on each rank
            for(unsigned int j=0;j<internal::atom_mask.size();j++){

               // get atom, spectrum and spin component to fourier transform
               atom=internal::atom_mask[j];
               spec=internal::spec_mask[j];
               sx = (*internal::sw_array[spec])[atom];
               index = k*internal::nspec+spec;

               // calculate spacial fourier transorm
               if (internal::prefactor == true){
                  skx_r[index] += sx*internal::cos_k[k*(internal::atom_mask.size())+j];
                  skx_i[index] += sx*internal::sin_k[k*(internal::atom_mask.size())+j];
               }
               else if (internal::prefactor == false){
                  skx_r[index] += sx*cos(rx[atom]*internal::kx[k] + ry[atom]*internal::kx[k] + rz[atom]*internal::kx[k]);
                  skx_i[index] += sx*sin(rx[atom]*internal::kx[k] + ry[atom]*internal::ky[k] + rz[atom]*internal::ky[k]);
               }
            }


            // reduce kpoint to rank for fft in time.
            if (internal::reduc_ver == "direct-scatter"){
               reduc_index=(k % nk_per_rank)*internal::nt*internal::nspec+time*internal::nspec;
               MPI_Reduce(&skx_r[k*internal::nspec], &skx_r_scatter[reduc_index], internal::nspec, MPI_DOUBLE, MPI_SUM, k / nk_per_rank, MPI_COMM_WORLD);
               MPI_Reduce(&skx_i[k*internal::nspec], &skx_i_scatter[reduc_index], internal::nspec, MPI_DOUBLE, MPI_SUM, k / nk_per_rank, MPI_COMM_WORLD);
            }

         #else

            for(unsigned int j=0;j<internal::atom_mask.size();j++){

               // get atom, spectrum and spin component to fourier transform
               atom=internal::atom_mask[j];
               spec=internal::spec_mask[j];
               sx = (*internal::sw_array[spec])[atom];

               index = time*internal::nk*internal::nspec+k*internal::nspec+spec;

               // calculate spacial fourier transorm
               if (internal::prefactor == true){
                  spinwaves::skx_r_node[index] += sx*internal::cos_k[k*(internal::atom_mask.size())+j];
                  spinwaves::skx_i_node[index] += sx*internal::sin_k[k*(internal::atom_mask.size())+j];
               }
               else if (internal::prefactor == false){
                  skx_r_node[index] += sx*cos(rx[atom]*internal::kx[k] + ry[atom]*internal::kx[k] + rz[atom]*internal::kx[k]);
                  skx_i_node[index] += sx*sin(rx[atom]*internal::kx[k] + ry[atom]*internal::ky[k] + rz[atom]*internal::ky[k]);
               }
            }

         #endif

      }

      // reduce every kpoint to rank 0.
      #ifdef MPICF
         if (internal::reduc_ver == "rank0"){
            MPI_Reduce(&skx_r[0], &skx_r_node[time*internal::nk*internal::nspec], internal::nk*internal::nspec, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&skx_i[0], &skx_i_node[time*internal::nk*internal::nspec], internal::nk*internal::nspec, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
         }
      #endif

	   return;

   }

} // end of sw namespace
