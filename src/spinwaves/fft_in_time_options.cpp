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
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "spinwaves.hpp"
#include "unitcell.hpp"
#include "vmpi.hpp"

// sw module headers
#include "internal.hpp"

#ifdef FFT
   #include <fftw3.h>
#endif

namespace spinwaves {

    namespace internal {

        std::ofstream file_K_time_real;
        std::ofstream file_K_time_imag;
        int j2;

        #ifdef FFT

        const int real = 0;
        const int imag = 1;

        // =============================================================================================================================================
        // Calculate one sided spectrum ================================================================================================================
        // =============================================================================================================================================
        void one_sided_spectrum(fftw_complex *in){

            for (int j1 = 0; j1 < internal::nt/2; j1++){
                j2 = internal::nt-j1-1;
                in[j1][real] = in[j1][real] + in[j2][real];
                in[j1][imag] = in[j1][imag] + in[j2][imag];

            }
        }

        void complex_magnitude(fftw_complex *os){

            for (int j1 = 0; j1 < internal::nt; j1++){
                os[j1][real] = os[j1][real] * os[j1][real] + os[j1][imag] * os[j1][imag];
                os[j1][imag] = 0.0;

            }

        }


        void normalise_each_kpoint(fftw_complex *os){

            double largest = std::abs(os[0][real]);
            //double index;

            // Find largest value in k_z array - REAL
            for (int j1 = 1; j1 < internal::nt; j1++){

                if (largest < std::abs(os[j1][real])){
                    largest = std::abs(os[j1][real]);
                    //index = j1;
                }
                }
                // normlise each value
                for (int j1 = 0; j1 < internal::nt; j1++){
                os[j1][real] /= largest;
            }

            largest = std::abs(os[0][imag]);

            // Find largest value in k_z array - COMPLEX
            for (int j1 = 1; j1 < internal::nt; j1++){
                if (largest < std::abs(os[j1][imag])){
                    largest = std::abs(os[j1][imag]);
                    //index = j1;
                }
                }
                // normlise each value
                for (int j1 = 0; j1 < internal::nt; j1++){
                os[j1][imag] /= largest;
            }


        }

        void write_to_file(fftw_complex *os, int k, int spec){

            std::stringstream sstr_real;
            std::stringstream sstr_imag;

            if (cm[spec] == true){

                // output filename
                #ifdef MPICF
                    sstr_real << "sw_spec_" << std::setw(2) << std::setfill('0') << spec+1 << "_real_" << std::setw(4) << std::setfill('0') << std::to_string(k+vmpi::my_rank*spinwaves::nk_per_rank) << ".dat";
                #else
                    sstr_real << "sw_spec_" << std::setw(2) << std::setfill('0') << spec+1 << "_real_" << std::setw(4) << std::setfill('0') << std::to_string(k) << ".dat";
                #endif

                // open file
                file_K_time_real.open(sstr_real.str());

                 if (oss[spec] == true){
                    for (int time=0; time < internal::nt/2; time++){
                        file_K_time_real << os[time][real] << "\n";
                    }
                }
                else if (oss[spec] == false){
                    for (int time=0; time < internal::nt; time++){
                    file_K_time_real << os[time][real] << "\n";
                    }
                }

                // close file
                file_K_time_real.close();

            }
            else if (cm[spec] == false){

                // output filenames
                #ifdef MPICF
                    sstr_real << "sw_spec_" << std::setw(2) << std::setfill('0') << spec+1 << "_real_" << std::setw(4) << std::setfill('0') << std::to_string(k+vmpi::my_rank*spinwaves::nk_per_rank) << ".dat";
                    sstr_imag << "sw_spec_" << std::setw(2) << std::setfill('0') << spec+1 << "_imag_" << std::setw(4) << std::setfill('0') << std::to_string(k+vmpi::my_rank*spinwaves::nk_per_rank) << ".dat";
                #else
                    sstr_real << "sw_spec_" << std::setw(2) << std::setfill('0') << spec+1 << "_real_" << std::setw(4) << std::setfill('0') << std::to_string(k) << ".dat";
                    sstr_imag << "sw_spec_" << std::setw(2) << std::setfill('0') << spec+1 << "_imag_" << std::setw(4) << std::setfill('0') << std::to_string(k) << ".dat";
                #endif

                // open files
                file_K_time_real.open(sstr_real.str());
                file_K_time_imag.open(sstr_imag.str());

                if (oss[spec] == true){
                    for (int time=0; time < internal::nt/2; time++){
                        file_K_time_real << os[time][real] << "\n";
                        file_K_time_imag << os[time][imag] << "\n";
                    }
                }
                else if (oss[spec] == false){
                    for (int time=0; time < internal::nt; time++){
                        file_K_time_real << os[time][real] << "\n";
                        file_K_time_imag << os[time][imag] << "\n";
                    }
                }
            }

            file_K_time_real.close();
            file_K_time_imag.close();
        }

        void write_intermediate_to_file(fftw_complex *os, int k, int spec){

            std::stringstream sstr_real;
            // std::stringstream sstr_imag;

                // output filenames
                #ifdef MPICF
                    sstr_real << "sw_intermediate_spec_" << std::setw(2) << std::setfill('0') << spec+1 << "_real_" << std::setw(4) << std::setfill('0') << std::to_string(k+vmpi::my_rank*spinwaves::nk_per_rank) << ".dat";
                    // sstr_imag << "k_intermediate_imag_" << std::setw(4) << std::setfill('0') << std::to_string(k+vmpi::my_rank*spinwaves::nk_per_rank) << ".dat";
                #else
                    sstr_real << "sw_intermediate_spec_" << std::setw(2) << std::setfill('0') << spec+1 << "_real_" << std::setw(4) << std::setfill('0') << std::to_string(k) << ".dat";
                    // sstr_imag << "k_intermediate_imag_" << std::setw(4) << std::setfill('0') << std::to_string(k) << ".dat";
                #endif

                // open files
                file_K_time_real.open(sstr_real.str());
                // file_K_time_imag.open(sstr_imag.str());

                for (int time=0; time < internal::nt; time++){
                    file_K_time_real << os[time][real] << "\n";
                    // file_K_time_imag << os[time][imag] << "\n";
                }

            file_K_time_real.close();
            // file_K_time_imag.close();
        }

        #endif // end of FFT macro
    }

}
