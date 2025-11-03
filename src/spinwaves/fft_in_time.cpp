//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Joel Hirst 2023. All rights reserved.
//
//   Email: j.r.hirst@shu.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "spinwaves.hpp"

// sw module headers
#include "internal.hpp"
#include "stopwatch.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include "vio.hpp"

//sergiu for SW
#include "unitcell.hpp"
#include "vector"
#include "vmpi.hpp"
#include <cmath>
#include "fstream"
#include "atoms.hpp"

#ifdef FFT
#include <fftw3.h>
#endif

namespace spinwaves {



   void fft_in_time(){

      #ifdef FFT

      // Start time for time series fourier transform
      // Set timer for runtime
      std::cout     << "Starting calculation of time series Discrete Fourier Transform." << std::endl;
      zlog << zTs() << "Starting calculation of time series Discrete Fourier Transform." << std::endl;
      stopwatch_t fft_stopwatch;
      fft_stopwatch.start();

      // real and imag indexing
      const int real = 0;
      const int imag = 1;

      //------------------------------------------------------------------------------------
      //fftw_complex combined_real_imag[internal::nt]; //- can lead to stack overflow
      //fftw_complex combined_real_imag_fftd[internal::nt];
      // same using allocated memory
      fftw_complex *combined_real_imag = fftw_alloc_complex(internal::nt);
      fftw_complex *combined_real_imag_fftd = fftw_alloc_complex(internal::nt);
      fftw_plan fft_in_time = fftw_plan_dft_1d(internal::nt, &combined_real_imag[0], &combined_real_imag_fftd[0], FFTW_FORWARD, FFTW_MEASURE);

      #ifdef MPICF

         // // make a mask that determines which kpoints will be dft'd. This prevents dft of lots of 0s on last rank for parallel implementation
         std::vector<int> kmask;
         kmask.resize(nk_per_rank * vmpi::num_processors,0);

         for (int i = 0; i < internal::nk; i++){
            kmask[i] = 1;
         }
         std::cout     << "Created mask for k-points." << std::endl;
         zlog << zTs() << "Created mask for k-points." << std::endl;




         // split comm world for time series fft if nranks > nk
         MPI_Comm fft_com;
         int color=1;
         if (vmpi::my_rank < (internal::nk + nk_per_rank - 1)/nk_per_rank) color = 0;

         MPI_Comm_split(MPI_COMM_WORLD, color, vmpi::my_rank, &fft_com);
         // We are in the new communicator
         int new_rank, new_size;
         MPI_Comm_rank(fft_com, &new_rank);
         MPI_Comm_size(fft_com, &new_size);


         if (color == 0){
            if (internal::reduc_ver == "rank0"){

               // rearrange array for scatter - I think this bit needs to be done for serial and parallel
               if (vmpi::my_rank == 0){
                  for (int k = 0; k < internal::nk; k++){
                        for (int time=0; time < internal::nt; time++){
                           for (int spec  = 0; spec < internal::nspec; spec++){
                           // Fill FFTW arrays with values from spacial FFT.
                           skx_r_node_transposed[k*internal::nt*internal::nspec + time*internal::nspec+spec] = skx_r_node[time*internal::nk*internal::nspec + k*internal::nspec+spec];
                           skx_i_node_transposed[k*internal::nt*internal::nspec + time*internal::nspec+spec] = skx_i_node[time*internal::nk*internal::nspec + k*internal::nspec+spec];
                        }
                     }
                  }
               }
               std::cout << skx_r_node_transposed.size() << " " << scatterlength << std::endl;
               std::cout << skx_r_scatter.size() << std::endl;
               std::cout << new_size << std::endl;
               MPI_Scatter(&skx_r_node_transposed[0], scatterlength, MPI_DOUBLE, &skx_r_scatter[0], scatterlength, MPI_DOUBLE, 0, fft_com);
               MPI_Scatter(&skx_i_node_transposed[0], scatterlength, MPI_DOUBLE, &skx_i_scatter[0], scatterlength, MPI_DOUBLE, 0, fft_com);
               // for (int j = 0; j < skx_r_scatter.size(); j++){
               //    std::cout << new_rank << " " << new_size << " " << skx_r_scatter[j] << std::endl;

               // }
               std::cout     << "Time-series FFT distrubuted amongst each rank." << std::endl;
               zlog << zTs() << "Time-series FFT distrubuted amongst each rank." << std::endl;
            }

            // loop over the k-points on each rank
            for (int k = 0; k < nk_per_rank; k++){
               for (int spec  = 0; spec < internal::nspec; spec++){
               // if the mask is 1, the calculate the dft.
               if (kmask[k+new_rank*nk_per_rank] == 1){


                  // populate fftw_complex vector
                  for (int time=0; time < internal::nt; time++){
                     int index = k*internal::nt*internal::nspec + time*internal::nspec + spec;
                     combined_real_imag[time][real] = skx_r_scatter[index];
                     combined_real_imag[time][imag] = skx_i_scatter[index];
                  }

                  // exexcute the fft
                  fftw_execute(fft_in_time);

                  // write intermediate structure factor to file
                  if (internal::isf == true) {
                     spinwaves::internal::complex_magnitude(combined_real_imag);
                     spinwaves::internal::write_intermediate_to_file(combined_real_imag, k, spec);
                  }

                  // calculate one sided spectrum
                  if (internal::oss[spec] == true) spinwaves::internal::one_sided_spectrum(combined_real_imag_fftd);

                  // normalise to largest amplitude for each k-value
                  if (internal::cm[spec] == true) spinwaves::internal::complex_magnitude(combined_real_imag_fftd);

                  // normalise to largest amplitude for each k-value
                  if (internal::normk[spec] == true) spinwaves::internal::normalise_each_kpoint(combined_real_imag_fftd);

                  // write each k-value to file
                  spinwaves::internal::write_to_file(combined_real_imag_fftd, k, spec);

               }
            }
         }
      }

      MPI_Comm_free(&fft_com);

      #else

         // rearrange array for scatter - I think this bit needs to be done for serial and parallel
         for (int k = 0; k < internal::nk; k++){
               for (int time=0; time < internal::nt; time++){
                  for (int spec  = 0; spec < internal::nspec; spec++){

                  // Fill FFTW arrays with values from spacial FFT.
                  skx_r_node_transposed[k*internal::nt*internal::nspec + time*internal::nspec+spec] = skx_r_node[time*internal::nk*internal::nspec + k*internal::nspec+spec];
                  skx_i_node_transposed[k*internal::nt*internal::nspec + time*internal::nspec+spec] = skx_i_node[time*internal::nk*internal::nspec + k*internal::nspec+spec];

               }
            }
         }

         for (int k = 0; k < internal::nk; k++){
            for (int spec  = 0; spec < internal::nspec; spec++){
               for (int time=0; time < internal::nt; time++){

                  // Fill FFTW arrays with values from spacial FFT.
                  int index = k*internal::nt*internal::nspec + time*internal::nspec + spec;
                  combined_real_imag[time][real] = skx_r_node_transposed[index];
                  combined_real_imag[time][imag] = skx_i_node_transposed[index];

               }

               // exexcute the fft
               fftw_execute(fft_in_time);

               // write intermediate structure factor to file
               if (internal::isf == true) {
                  spinwaves::internal::complex_magnitude(combined_real_imag);
                  spinwaves::internal::write_intermediate_to_file(combined_real_imag, k, spec);
               }

               // calculate one sided spectrum
               if (internal::oss[spec] == true) spinwaves::internal::one_sided_spectrum(combined_real_imag_fftd);

               // normalise to largest amplitude for each k-value
               if (internal::cm[spec] == true) spinwaves::internal::complex_magnitude(combined_real_imag_fftd);

               // normalise to largest amplitude for each k-value
               if (internal::normk[spec] == true) spinwaves::internal::normalise_each_kpoint(combined_real_imag_fftd);

               // write each k-value to file
               spinwaves::internal::write_to_file(combined_real_imag_fftd, k, spec);
            }
         }


      #endif

      std::cout     << "Total duration of time series discrete fourier transform [s]: " << fft_stopwatch.elapsed_seconds() << std::endl;
      zlog << zTs() << "Total duration of time series discrete fourier transform [s]: " << fft_stopwatch.elapsed_seconds() << std::endl;

      // destroy fftw3 plan
      fftw_destroy_plan(fft_in_time);

      // free complex alocated arrays
      fftw_free(combined_real_imag);
      fftw_free(combined_real_imag_fftd);

      #endif // end of FFT macro

   }

} // end of spinwaves namespace
