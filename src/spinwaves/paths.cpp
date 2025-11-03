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
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

// Vampire headers
#include "spinwaves.hpp"
#include "unitcell.hpp"
#include "internal.hpp"
#include "errors.hpp"
#include "atoms.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
#include "sim.hpp"

namespace spinwaves{

    namespace internal {

        void initialise_arrays(){

            // define number of time and kpoitns
            internal::nk = spinwaves::internal::kx.size();
            internal::nt = (sim::total_time / sim::partial_time);
            std::cout << "total number of spinwave kpoints "    << internal::nk << std::endl;
            std::cout << "total number of spinwave timepoints " << internal::nt << std::endl;


            // Spinwaves have to be calculated for every k on every rank
            spinwaves::skx_r.resize(internal::nk*internal::nspec,0.0);
            spinwaves::skx_i.resize(internal::nk*internal::nspec,0.0);

            // calculate memory requirements
            double mem = (4.0 * internal::nt * internal::nk * internal::nspec * sizeof(double)) / 1.0e6;
            std::cout << "Spinwave module using  " << mem << " MB of RAM per node for k vs time arrays for a simulation time of " << sim::total_time;
            std::cout << " timesteps with the dynamic structure factor calculated every " << sim::partial_time << " timesteps." << std::endl;

            zlog << zTs() << "Spinwave module using  " << mem << " MB of RAM per node for k vs time arrays for a simulation time of " << sim::total_time
            << " timesteps with the dynamic structure factor calculated every " << sim::partial_time << " timesteps." << std::endl;


            // the array containing the the values of k for every timepoint only needs to be stored on rank0 for each node. It also needs to be transposed
            #ifdef MPICF

                nk_per_rank = std::ceil(static_cast<double>(internal::nk) / static_cast<double>(vmpi::num_processors));
                std::cout << "nk_per_rank " << nk_per_rank << std::endl;
                scatterlength = nk_per_rank * internal::nt * internal::nspec;
                std::cout << "scatterlength " << scatterlength << std::endl;
                skx_r_scatter.resize(scatterlength,0.0);
                skx_i_scatter.resize(scatterlength,0.0);

                // resize arrays
                if (vmpi::my_rank == 0){

                    int lennode = internal::nt * internal::nspec * internal::nk;
                    spinwaves::skx_r_node.resize(lennode,0.0);
                    spinwaves::skx_i_node.resize(lennode,0.0);

                    // need to pad for scatter of time-series data to each ranks.
                    int lentrans = internal::nt * internal::nspec * ((internal::nk + nk_per_rank - 1)/nk_per_rank) * nk_per_rank;
                    spinwaves::skx_r_node_transposed.resize(lentrans,0.0);
                    spinwaves::skx_i_node_transposed.resize(lentrans,0.0);
                }

            #else
                spinwaves::skx_r_node.resize(internal::nt*internal::nk*internal::nspec,0.0);
                spinwaves::skx_i_node.resize(internal::nt*internal::nk*internal::nspec,0.0);
                spinwaves::skx_r_node_transposed.resize(internal::nt*internal::nk*internal::nspec,0.0);
                spinwaves::skx_i_node_transposed.resize(internal::nt*internal::nk*internal::nspec,0.0);
            #endif
        }

        void determine_path(){

            std::ifstream file_read_K_path;
            file_read_K_path.open(spinwaves::internal::filename);


            // else read the path from the file specified in interface.cpp
            if (file_read_K_path.is_open()){

                std::cout << "Spinwaves path file " << spinwaves::internal::filename << " succesfully opened." << std::endl;
                zlog << zTs()  << "Spinwaves path file " << spinwaves::internal::filename << " succesfully opened." << std::endl;

                // read contents of path file
                int linecount = 0;
                std::string line;
                while (std::getline(file_read_K_path, line)) {
                    std::istringstream iss(line);
                    linecount++;
                    double val1, val2, val3;
                    if (iss >> val1 >> val2 >> val3) {
                        spinwaves::internal::pathx.push_back(val1);
                        spinwaves::internal::pathy.push_back(val2);
                        spinwaves::internal::pathz.push_back(val3);
                    }
                    else {
                        std::cerr << "Error reading line: " << line << std::endl;
                    }
                }

                std::cout << "Spinwaves path file " << spinwaves::internal::filename << " contains " << linecount << " lines." << std::endl;
                zlog << zTs()  << "Spinwaves path file" << spinwaves::internal::filename << " contains " << linecount << " lines." << std::endl;

                // close file and check
                file_read_K_path.close();
                std::cout << "Spinwaves path file " << spinwaves::internal::filename << " has been closed." << std::endl;
                zlog << zTs()  << "Spinwaves path file " << spinwaves::internal::filename << " has been closed." << std::endl;

            }
            // if filename has not been found in spiwnaves/interface.cpp use a built in path
            else if (!file_read_K_path.is_open() && uc::sw_crystal_structure!=""){
                std::cout << "sw_crystal_structure = " << uc::sw_crystal_structure<< std::endl;

                // Determine which path to take depending on crystal structure
                if(uc::sw_crystal_structure      == "sc"                ) spinwaves::internal::path_sc();
                else if(uc::sw_crystal_structure == "bcc"               ) spinwaves::internal::path_bcc();
                else if(uc::sw_crystal_structure == "bcc-110"           ) spinwaves::internal::path_bcc();
                else if(uc::sw_crystal_structure == "fcc"               ) spinwaves::internal::path_fcc();
                else if(uc::sw_crystal_structure == "fcc-111"           ) spinwaves::internal::path_fcc();
                else if(uc::sw_crystal_structure == "hcp"               ) spinwaves::internal::path_hcp();
                else if(uc::sw_crystal_structure == "heusler"           ) spinwaves::internal::path_fcc();
                // else if(uc::sw_crystal_structure == "honeycomb"      ) spinwaves::internal::path_honeycomb();
                // else if(uc::sw_crystal_structure == "alpha-honeycomb") spinwaves::internal::path_honeycomb_alpha();
                // else if(uc::sw_crystal_structure == "beta-honeycomb" ) spinwaves::internal::path_honeycomb_beta();
                // else if(uc::sw_crystal_structure == "kagome"         ) spinwaves::internal::path_kagome();
                else if(uc::sw_crystal_structure == "mn2au"             ) spinwaves::internal::path_mn2au();
                // else if(uc::sw_crystal_structure == "NdFeB"          ) spinwaves::internal::path_NdFeB();
                else if(uc::sw_crystal_structure == "rocksalt"          ) spinwaves::internal::path_rocksalt();
                // else if(uc::sw_crystal_structure == "spinel"         ) spinwaves::internal::path_spinel();
                // else if(uc::sw_crystal_structure == "spinel-layered" ) spinwaves::internal::path_spinel_layered();
                // else if(uc::sw_crystal_structure == "SmFeN"          ) spinwaves::internal::path_SmFeN();
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error: Unknown spinwaves crystal_type "<< uc::sw_crystal_structure << " found during spinwave path initialisation. Exiting." << std::endl;
                    terminaltextcolor(WHITE);
                    zlog << zTs() << "Error: Unknown spinwaves crystal_type "<< uc::sw_crystal_structure << " found during spinwave path initialisation. Exiting." << std::endl;
                    err::vexit();
                }
            }
            else {
                terminaltextcolor(RED);
                std::cerr << "Error: Cannot find spinwaves path file or crystal structure. Exiting." << std::endl;
                std::cerr << "filename for kpath has been specified as: \"" << spinwaves::internal::filename << "\"." << std::endl;
                std::cerr << "Crystal structure has been defined as: \"" << uc::sw_crystal_structure << "\"." << std::endl;
                terminaltextcolor(WHITE);
                zlog << zTs() << "Error: Cannot find spinwaves path file \""<< spinwaves::internal::filename << "\"." << std::endl;
                zlog << zTs() << "Error: Cannot find spinwaves path file or crystal structure. Exiting." << std::endl;
                zlog << zTs() << "filename for kpath has been specified as: \"" << spinwaves::internal::filename << std::endl;
                zlog << zTs() << "Crystal structure has been defined as: \"" << uc::sw_crystal_structure << "\"." << std::endl;
                err::vexit();
            }
        }



         void determine_kpoints_from_user_specific_k(const double system_dimensions_x,
                    const double system_dimensions_y,
                    const double system_dimensions_z,
					const double unit_cell_size_x,
					const double unit_cell_size_y,
					const double unit_cell_size_z){


            // kfile
            std::ofstream kfile;
            kfile.open("kvectors.out");


            // number of kpoints
            int len=spinwaves::internal::pathx.size();

            // reciprocal lattice vectors
            double b[3];
            b[0] = 2.0 * M_PI * (unit_cell_size_y * unit_cell_size_z) / (unit_cell_size_x * (unit_cell_size_y * unit_cell_size_z));
            b[1] = 2.0 * M_PI * (unit_cell_size_z * unit_cell_size_x) / (unit_cell_size_x * (unit_cell_size_y * unit_cell_size_z));
            b[2] = 2.0 * M_PI * (unit_cell_size_x * unit_cell_size_y) / (unit_cell_size_x * (unit_cell_size_y * unit_cell_size_z));

            std::cout << "Determining k-points from user specific k-values..." << std::endl;

            // loop over number of kpoints
            for (int k = 0; k < len; k++){

                // convert from user defined values in units of 2pi/a to m^{-1} and push back to array
                double kxtemp = b[0] * pathx[k];
                double kytemp = b[1] * pathy[k];
                double kztemp = b[2] * pathz[k];

                spinwaves::internal::kx.push_back(kxtemp);
                spinwaves::internal::ky.push_back(kytemp);
                spinwaves::internal::kz.push_back(kztemp);

                // save to kfile
                kfile << kxtemp << " " << kytemp << " " << kztemp << "\n";

            }

            kfile.close();

        }

        void determine_kpoints_from_user_high_sym_path(const double system_dimensions_x,
                    const double system_dimensions_y,
                    const double system_dimensions_z,
					const double unit_cell_size_x,
					const double unit_cell_size_y,
					const double unit_cell_size_z){

            // kfile
            std::ofstream kfile;
            kfile.open("kvectors.out");

            int len=spinwaves::internal::pathx.size();

            // check path file contains an even number of high symm points.
            if ( len % 2 != 0){
                terminaltextcolor(RED);
                std::cerr <<     "Error: k-path file contains an odd number of high symmetry points." << std::endl;
                terminaltextcolor(WHITE);
                zlog << zTs() << "Error: k-path file contains an odd number of high symmetry points." << std::endl;
                err::vexit();
            }

            std::cout << unit_cell_size_x << " " << unit_cell_size_y << " " << unit_cell_size_z << std::endl;


            // Get reciprocal lattice vectors from cubic unit cell JRH 26/10/23
            // Based on equations through link: http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/fourier/reciprocal_lattice.php
            // The assumption is made that the unit cell is cubic
            double b[3];
            b[0] = 2.0 * M_PI * (unit_cell_size_y * unit_cell_size_z) / (unit_cell_size_x * (unit_cell_size_y * unit_cell_size_z));
            b[1] = 2.0 * M_PI * (unit_cell_size_z * unit_cell_size_x) / (unit_cell_size_x * (unit_cell_size_y * unit_cell_size_z));
            b[2] = 2.0 * M_PI * (unit_cell_size_x * unit_cell_size_y) / (unit_cell_size_x * (unit_cell_size_y * unit_cell_size_z));

            zlog << zTs()  << "Reciprocal Lattice Vectors:" << std::endl;
            zlog << zTs()  << b[0] << " 0.0 0.0" << std::endl;
            zlog << zTs()  << "0.0 " << b[1] << " 0.0" << std::endl;
            zlog << zTs()  << "0.0 0.0 " << b[2] << std::endl;


            // Generate k-points
            std::cout << "Determining k-points from user specified high-symmetry path..." << std::endl;
            zlog << zTs()  << "Determining k-points from user specified high-symmetry path..." << std::endl;

            double kx, ky, kz;

            for (int row=0; row<len; row+=2){

                double cellx = system_dimensions_x/unit_cell_size_x;
                double celly = system_dimensions_y/unit_cell_size_y;
                double cellz = system_dimensions_z/unit_cell_size_z;

                double kx0=spinwaves::internal::pathx[row];
                double ky0=spinwaves::internal::pathy[row];
                double kz0=spinwaves::internal::pathz[row];
                double kx1=spinwaves::internal::pathx[row+1];
                double ky1=spinwaves::internal::pathy[row+1];
                double kz1=spinwaves::internal::pathz[row+1];

                // std::cout << "kx " << kx0 << " " << kx1 << std::endl;
                // std::cout << "ky " << ky0 << " " << ky1 << std::endl;
                // std::cout << "kz " << kz0 << " " << kz1 << std::endl;
                double distancex = cellx * (kx1-kx0);
                double distancey = celly * (ky1-ky0);
                double distancez = cellz * (kz1-kz0);


                int res1 = spinwaves::internal::gcd(static_cast<int>(round(distancex)), static_cast<int>(round(distancey)));
                int common_denom = spinwaves::internal::gcd(res1, static_cast<int>(round(distancez)));

                std::cout << "Spinwaves module has found " << common_denom << " k-points between locations ";
                std::cout << "[" << kx0 << ", " << ky0 << ", " << kz0 << "] and ";
                std::cout << "[" << kx1 << ", " << ky1 << ", " << kz1 << "]" << std::endl;
                zlog << zTs() << "Spinwaves module has found " << common_denom << " k-points between locations ";
                zlog << "[" << kx0 << ", " << ky0 << ", " << kz0 << "] and ";
                zlog << "[" << kx1 << ", " << ky1 << ", " << kz1 << "]" << std::endl;

                kx=kx0;
                ky=ky0;
                kz=kz0;

                // make sure you always finish at the points specified by user.

                // loop over points between two locations
                for (int k = 0; k < common_denom; k++){

                    kfile << b[0] * kx << " " << b[1] * ky << " " << b[2] * kz << std::endl;

                    // make sure to convert to units of 2pi/latconst
                    spinwaves::internal::kx.push_back(b[0] * kx);
                    spinwaves::internal::ky.push_back(b[1] * ky);
                    spinwaves::internal::kz.push_back(b[2] * kz);

                    kx = kx + distancex/static_cast<double>(common_denom)/cellx;
                    ky = ky + distancey/static_cast<double>(common_denom)/celly;
                    kz = kz + distancez/static_cast<double>(common_denom)/cellz;

                }

            }

            // push back final row
            kfile << b[0] * kx << " " << b[1] * ky << " " << b[2] * kz << std::endl;
            spinwaves::internal::kx.push_back(b[0] * kx);
            spinwaves::internal::ky.push_back(b[1] * ky);
            spinwaves::internal::kz.push_back(b[2] * kz);
        }
   }
}
