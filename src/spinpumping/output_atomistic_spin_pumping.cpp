//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard Evans 2022. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <iostream>

// Vampire headers
#include "errors.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "spinpumping.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"

// internal module headers
#include "internal.hpp"

namespace spin_pumping{

   namespace internal{

      //------------------------------------------------------------------------
      // file-local storage variables
      //------------------------------------------------------------------------
      bool output_atomistic_spin_pumping_initialised = false; // bool to determine if storage has been allocated
      bool output_cells_spin_pumping_initialised = false; // bool to determine if storage has been allocated

      std::vector<double> data_from_all_atoms(0);   // 3N array to store atomic position/spin_current from all processors
      std::vector<double> data_from_local_atoms(0); // 3N array to store atomic position/spin_current on local processor
      std::vector<int> counts(0);                   // number of components received from each processor
      std::vector<int> displacements(0);            // offsets for data from remote processors

		
      //------------------------------------------------------------------------
      // Function to write meta data for coordinate and atomic data
      //------------------------------------------------------------------------
      void write_atomistic_coordinates_meta(){

         std::ofstream scmf; // spin coordinate meta file
         scmf.open("atoms-spin-pumping-coords.meta");

         // determine file format
         std::string format_string;
         format_string = "text";

         // Get system date
         time_t rawtime = time(NULL);
         struct tm * timeinfo = localtime(&rawtime);

         // get number of atoms (local CPU)
         const uint64_t local_num_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
         // get total number of atoms on master process
         uint64_t total_num_atoms = vmpi::reduce_sum(local_num_atoms);

         scmf << "#----------------------------------------------------------"<< std::endl;
         scmf << "# Atomistic coordinates configuration file for vampire V5+"<< std::endl;
         scmf << "#----------------------------------------------------------"<< std::endl;
         scmf << "# Date: "<< asctime(timeinfo);
         scmf << "#--------------------------------------------"<< std::endl;
         scmf << "Format: "<< format_string << std::endl;
         scmf << "#--------------------------------------------"<< std::endl;
         scmf << "Number of atoms: "<< total_num_atoms << std::endl;
         scmf << "#--------------------------------------------" << std::endl;
         scmf << "Number of materials: " << mp::num_materials << std::endl;
         for(int mat=0;mat<mp::num_materials;mat++){
            scmf << mat << "\t" << mp::material[mat].mu_s_SI/9.274e-24 << "\t" << mp::material[mat].element << "\t" <<
            mp::material[mat].name << std::endl;
         }
         scmf << "#--------------------------------------------" << std::endl;
         scmf << "Number of coord files: " << "1" << std::endl;
         scmf << "atoms-coords.data" << std::endl;
         scmf << std::flush;

         // close file
         scmf.close();

      }

      //------------------------------------------------------------------------
      // Function to output calculated atomistic coords on root process
      //------------------------------------------------------------------------
      void write_atomistic_coordinates_data(const std::vector<double>& x_coord_array, // atomic coordinates (angstroms)
                                        const std::vector<double>& y_coord_array,
                                        const std::vector<double>& z_coord_array){

         // get number of atoms (local CPU)
         const uint64_t local_num_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;

         // get total number of atoms on master process
         uint64_t total_num_atoms = vmpi::reduce_sum(local_num_atoms);

         // resize local arrays on all processes
         data_from_local_atoms.resize(3*local_num_atoms);

         // resize global array on root process
         if(vmpi::master) data_from_all_atoms.resize(3*total_num_atoms);

         // copy local position data to local array
         for(uint64_t atom = 0; atom < local_num_atoms; atom++){
            data_from_local_atoms[ 3 * atom + 0 ] = x_coord_array[atom];
            data_from_local_atoms[ 3 * atom + 1 ] = y_coord_array[atom];
            data_from_local_atoms[ 3 * atom + 2 ] = z_coord_array[atom];
            // data_from_local_atoms[ 4 * atom + 3 ] = moments_array[atom];
         }

         // set up counts and displacements
         vmpi::counts_and_displacements(data_from_local_atoms, data_from_all_atoms, counts, displacements);

         // collate global position data on master
         vmpi::fast_collate(data_from_local_atoms, data_from_all_atoms, counts, displacements);

         // compute average coordinates
         double mean_x = 0.0;
         double mean_y = 0.0;
         double mean_z = 0.0;
         for(uint64_t atom = 0; atom < total_num_atoms; atom++){
            mean_x += data_from_all_atoms[ 3 * atom + 0 ];
            mean_y += data_from_all_atoms[ 3 * atom + 1 ];
            mean_z += data_from_all_atoms[ 3 * atom + 2 ];
         }
         mean_x = mean_x / double (total_num_atoms);
         mean_y = mean_y / double (total_num_atoms);
         mean_z = mean_z / double (total_num_atoms);

         if(vmpi::master){
            std::ofstream ofile;
            std::stringstream filename;
            filename << "atoms-spin-pumping-coords.data";
            zlog << zTs() << "Outputting atomistic coordinates for spin pumping to file " << filename.str() << std::endl;
            ofile.open(std::string(filename.str()).c_str());

            ofile << total_num_atoms << "\n";
            for(uint64_t atom = 0; atom < total_num_atoms; atom++){
               ofile << data_from_all_atoms[ 3 * atom + 0 ] - mean_x << "\t" <<
                        data_from_all_atoms[ 3 * atom + 1 ] - mean_y << "\t" <<
                        data_from_all_atoms[ 3 * atom + 2 ] - mean_z << "\t\n";
            }
            ofile.close();
         }

         // Every process wait here
         vmpi::barrier();

         // re-calculate counts and displacements for 3-vector data
         data_from_local_atoms.resize(4*local_num_atoms);
         if(vmpi::master) data_from_all_atoms.resize(4*total_num_atoms);
         vmpi::counts_and_displacements(data_from_local_atoms, data_from_all_atoms, counts, displacements); // set up counts and displacements

         // set initialised variable to true
         output_atomistic_spin_pumping_initialised = true;

         return;

      }

      //------------------------------------------------------------------------
      // Function to call output of coordinate of atomic data
      //------------------------------------------------------------------------
      void output_atomistic_coordinates(const std::vector<double>& x_coord_array, // atomic coordinates (angstroms)
                                        const std::vector<double>& y_coord_array,
                                        const std::vector<double>& z_coord_array){ 

         // write spin coordinates info to "meta" file
         write_atomistic_coordinates_meta();
         // Write spin coordinates to "data" file
         write_atomistic_coordinates_data(x_coord_array, y_coord_array, z_coord_array);
      }

      //---------------------------------------------------------------------
      // Function to write meta data for each configuration
      //---------------------------------------------------------------------
      void write_atomistic_spin_pumping_meta(const uint64_t config_file_counter){

         // determine file name
         std::stringstream filename;
         filename << "atoms-spin-pumping-";
         filename << std::setfill ('0') << std::setw (8) << config_file_counter;
         filename << ".meta";

         // Get system date
         time_t rawtime = time(NULL);
         struct tm * timeinfo = localtime(&rawtime);

         // calculate real time
         const double real_time = double(sim::time) * mp::dt_SI;

         // Declare and open output file
         std::ofstream ofile;
         ofile.open (filename.str());

         ofile << "#------------------------------------------------------"<< "\n";
         ofile << "# Spin pumping configuration file for vampire v5+"<< "\n";
         ofile << "#------------------------------------------------------"<< "\n";
         ofile << "# Date: "<< asctime(timeinfo);
         ofile << "#------------------------------------------------------"<< "\n";
         ofile << "Time: " << real_time << "\n";
         ofile << "Temperature: "<< sim::temperature << "\n";
         ofile << "#------------------------------------------------------" << "\n";
         ofile << "Number of spin files: " << "1" << "\n"; 
         ofile << "atoms-spin-pumping-" << std::setfill ('0') << std::setw (8) << config_file_counter << ".data" << std::endl;
         // flush data to disk
         ofile << std::flush;
         ofile << "#------------------------------------------------------"<< "\n";

         return;

      }

      //------------------------------------------------------------------------
      // Function to output calculated atomistic coords on root process
      //------------------------------------------------------------------------
      void write_atomistic_spin_pumping_data(const uint64_t config_file_counter,  // counter of output files
                                             const std::vector<double>& moments_array){ // atomistic magnetic moments (bohr magnetons)

         // get number of atoms (local CPU)
         const uint64_t local_num_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;

         // get total number of atoms on master process
         uint64_t total_num_atoms = vmpi::reduce_sum(local_num_atoms);

         // resize local arrays on all processes
         data_from_local_atoms.resize(4*local_num_atoms);

         // copy local  data to local array
         for(uint64_t atom = 0; atom < local_num_atoms; atom++){
            data_from_local_atoms[ 4 * atom + 0 ] = spin_pumping::internal::x_atom_spin_pumping_array[atom];
            data_from_local_atoms[ 4 * atom + 1 ] = spin_pumping::internal::y_atom_spin_pumping_array[atom];
            data_from_local_atoms[ 4 * atom + 2 ] = spin_pumping::internal::z_atom_spin_pumping_array[atom];
            data_from_local_atoms[ 4 * atom + 3 ] = moments_array[atom];
         }

         // collate global data on master
         vmpi::fast_collate(data_from_local_atoms, data_from_all_atoms, counts, displacements);

         if(vmpi::master){
            std::ofstream ofile;
            std::stringstream filename;
            filename << "atoms-spin-pumping-";
            filename << std::setfill ('0') << std::setw (8) << config_file_counter;
            filename << ".data";
            zlog << zTs() << "Outputting spin pumping data to file " << filename.str() << std::endl;
            ofile.open(std::string(filename.str()).c_str());

            ofile << total_num_atoms << "\n";
            for(int atom = 0; atom < total_num_atoms; atom++){
               ofile << data_from_all_atoms[ 4 * atom + 0 ] << "\t" <<
                        data_from_all_atoms[ 4 * atom + 1 ] << "\t" <<
                        data_from_all_atoms[ 4 * atom + 2 ] << "\t" <<
                        data_from_all_atoms[ 4 * atom + 3 ] << "\t\n";
            }
            ofile.close();
         }

         return;

      } // end function output_atomistic_spin_pumping

      //------------------------------------------------------------------------
      // Function to call output of spin pumping of atomic data
      //------------------------------------------------------------------------
      void output_atomistic_spin_pumping(const uint64_t config_file_counter,  // counter of output files
                                       const std::vector<double>& moments_array){ // atomistic magnetic moments (bohr magnetons)

         // check for initialised data structures
         if(!output_atomistic_spin_pumping_initialised){
            std::cerr << "Programmer error : spin_pumping::output_atomistic_spin_pumping_initialised() called before initialisation" << std::endl;
            err::vexit();
         }

         // write spin configurations info to "meta" file
         write_atomistic_spin_pumping_meta(config_file_counter);
         // Write spin configurations to "data" file
         write_atomistic_spin_pumping_data(config_file_counter,moments_array);

      }

      //------------------------------------------------------------------------
      // Function to output calculated atomistic coords on root process
      //------------------------------------------------------------------------
      void output_cells_spin_pumping(const uint64_t config_file_counter){
      } // end function output_cells_spin_pumping

} // end of namespace internal
} // end of namespace spin_current
