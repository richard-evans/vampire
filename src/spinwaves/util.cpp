//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sergiu Ruta 2022. All rights reserved.
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
#include <string>
#include <sstream>
#include <vector>

// Vampire headers
#include "spinwaves.hpp"
#include "vmpi.hpp"
#include "unitcell.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "atoms.hpp"
#include "program.hpp"

// sw module headers
#include "internal.hpp"

namespace spinwaves {

    namespace internal {

        void save_frequencies(){

            #ifdef MPICF
                if (vmpi::my_rank==0) {

                    // open files
                    std::ofstream freq_file;
                    freq_file.open("frequencies.dat");
                    for (int i=0; i < internal::nt; i++){
                        freq_file << i/(sim::partial_time*mp::dt/1.76e11)/internal::nt << "\n";
                    }
                    freq_file.close();
                }
            #else
                std::ofstream freq_file;
                freq_file.open("frequencies.dat");
                for (int i=0; i < internal::nt; i++){
                    freq_file << i/(sim::partial_time*mp::dt/1.76e11)/internal::nt << "\n";
                }
                freq_file.close();
            #endif

        }

        void check_numbering_of_spectrums(){
            std::sort(internal::super_index_values.begin(), internal::super_index_values.end());
            internal::super_index_values.erase(std::unique(internal::super_index_values.begin(), internal::super_index_values.end()), internal::super_index_values.end());

            //
            if (static_cast<unsigned>(internal::nspec) != internal::super_index_values.size()){
                terminaltextcolor(RED);
                std::cerr << "The value of spinwaves:number-of-spectrums does not agree with the values of \"spinwaves[x]:\" found in input file" << std::endl;
                std::cerr << internal::nspec << " " << internal::super_index_values.size() << std::endl;
                err::vexit();
            }

            for (int i = 0; i < internal::nspec; i++){
                if (internal::super_index_values[i] != i){
                    terminaltextcolor(RED);
                    std::cerr << "Difference identified between spinwaves::num-spectrums and ordering of spectrum specific parameters." << std::endl;
                    std::cerr << "Number of spectrums specified with spinwaves::num-spectrums: " << internal::nspec << std::endl;
                    std::cerr << "Additional/Missing parameters for spectrum \"spinwaves[" << internal::super_index_values[i] <<  "]:\" found in input file" << std::endl;
                    err::vexit();
                }
            }
        }

        void determine_spin_component(){

            internal::sw_array.resize(internal::nspec);
            // Set spin_array based on the component we want to calculate spinwave from
            for (int spec = 0; spec < internal::nspec; spec++){
                if (internal::component[spec] == "sx") {
                    internal::sw_array[spec] = &atoms::x_spin_array;
                } else if (internal::component[spec] == "sy") {
                    internal::sw_array[spec] = &atoms::y_spin_array;
                } else if (internal::component[spec] == "sz") {
                    internal::sw_array[spec] = &atoms::z_spin_array;
                }
            }
        }

        // Determine the mask for calculations of SW spectrum for specific materials
        void calculate_material_mask(){

            #ifdef MPICF

                // loop over number of materials in EACH SPECTRUM
                for (unsigned int i = 0; i < mat.size(); i++){



                    // check mat value does not exceed number of materials
                    if (mat[i] >= mp::num_materials){
                        terminaltextcolor(RED);
                        std::cerr << "Error: Numeric value of materials in spectrum " << internal::mat_in_spec[i] << " exceeds the total number of materials." << std::endl;
                        err::vexit();
                    }

                    for (int atom=0;atom<vmpi::num_core_atoms+vmpi::num_bdry_atoms;atom++){

                        if (atoms::type_array[atom] == internal::mat[i]){
                            spec_mask.push_back(internal::mat_in_spec[i]);
                            atom_mask.push_back(atom);
                        }
                    }

                }

            #else

                // loop over number of materials in EACH SPECTRUM
                for (unsigned int i = 0; i < mat.size(); i++){

                    // check mat value does not exceed number of materials
                    if (mat[i] >= mp::num_materials){
                        terminaltextcolor(RED);
                        std::cerr << "Error: Numeric value of materials in spectrum " << internal::mat_in_spec[i] << " exceeds the total number of materials." << std::endl;
                        err::vexit();
                    }

                    for (int atom=0;atom<atoms::num_atoms;atom++){

                        if (atoms::type_array[atom] == internal::mat[i]){
                            spec_mask.push_back(internal::mat_in_spec[i]);
                            atom_mask.push_back(atom);
                        }
                    }
                }
            #endif
        }

         void calculate_fourier_prefactor(const std::vector<double>& rx, const std::vector<double>& ry, const std::vector<double>& rz){

            double arg;
            int nk=spinwaves::internal::kx.size();

            #ifdef MPICF

                // calculate memory requirements
                if (vmpi::my_rank == 0){
                    double mem = 0.0;
                    mem = (2.0 * vmpi::num_processors * nk*(internal::atom_mask.size()) * sizeof(double)) / 1.0e6;
                    std::cout << "Spinwave module using " << mem << " MB of RAM for real-space fourier prefactor." << std::endl;

                }
                // initialise array size for prefactor
                spinwaves::internal::cos_k.resize(nk*(internal::atom_mask.size()),0.0);
                spinwaves::internal::sin_k.resize(nk*(internal::atom_mask.size()),0.0);

                for (int k=0; k < nk; k++){

                    double kx = spinwaves::internal::kx[k];
                    double ky = spinwaves::internal::ky[k];
                    double kz = spinwaves::internal::kz[k];

                    for(unsigned int j=0;j<internal::atom_mask.size();j++){
                        int atom=internal::atom_mask[j];
                        arg=-kx*rx[atom]-ky*ry[atom]-kz*rz[atom];
                        spinwaves::internal::cos_k[k*(internal::atom_mask.size())+j] = cos(arg);
                        spinwaves::internal::sin_k[k*(internal::atom_mask.size())+j] = sin(arg);
                    }
                }

            #else

                // calculate memory requirements
                double mem = 0.0;
                mem = (2.0 *  nk*internal::atom_mask.size() * sizeof(double)) / 1.0e6;
                std::cout << "Spinwave module using " << mem << " MB of RAM for real-space fourier prefactor." << std::endl;

                spinwaves::internal::cos_k.resize(nk*internal::atom_mask.size(),0.0);
                spinwaves::internal::sin_k.resize(nk*internal::atom_mask.size(),0.0);


                for (int k=0; k < nk; k++){

                double kx = spinwaves::internal::kx[k];
                double ky = spinwaves::internal::ky[k];
                double kz = spinwaves::internal::kz[k];

                    for(unsigned int j=0;j<internal::atom_mask.size();j++){
                        int atom=internal::atom_mask[j];
                        arg=-kx*rx[atom]-ky*ry[atom]-kz*rz[atom];
                        spinwaves::internal::cos_k[k*internal::atom_mask.size()+j] = cos(arg);
                        spinwaves::internal::sin_k[k*internal::atom_mask.size()+j] = sin(arg);
                    }
                }

            #endif

        }






    }
}
