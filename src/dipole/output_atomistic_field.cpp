//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <iostream>

// Vampire headers
#include "dipole.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// dipole module headers
#include "internal.hpp"

// alias interal dipole namespace for brevity
namespace dp = dipole::internal;

namespace dipole{

   namespace internal{

      //------------------------------------------------------------------------
      // file-local storage variables
      //------------------------------------------------------------------------
      bool output_atomistic_field_initialised = false; // bool to determine if storage has been allocated

      std::vector<double> data_from_all_atoms(0);   // 3N array to store atomic dipole position/field from all processors
      std::vector<double> data_from_local_atoms(0); // 3N array to store atomic dipole position/field on local processor
      std::vector<int> counts(0);                   // number of components received from each processor
      std::vector<int> displacements(0);            // offsets for data from remote processors

      //------------------------------------------------------------------------
      // Function to output calculated atomistic dipole coords on root process
      //------------------------------------------------------------------------
      void output_atomistic_coordinates(int num_atoms,                      // number of atoms (only correct in serial)
                                        std::vector<double>& x_coord_array, // atomic coordinates (angstroms)
                                        std::vector<double>& y_coord_array,
                                        std::vector<double>& z_coord_array,
                                        std::vector<double>& moments_array){ // atomistic magnetic moments (bohr magnetons)

         // inform user that atomistic fields are being outputted
         zlog << zTs() << "Outputting atomistic dipole fields to file" << std::endl;

         // get number of atoms (local CPU)
         const uint64_t local_num_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;

         // get total number of atoms on master process
         uint64_t total_num_atoms = vmpi::reduce_sum(local_num_atoms);

         // resize local arrays on all processes
         data_from_local_atoms.resize(4*local_num_atoms);

         // resize global array on root process
         if(vmpi::master) data_from_all_atoms.resize(4*total_num_atoms);

         // copy local position data to local array
         for(uint64_t atom = 0; atom < local_num_atoms; atom++){
            data_from_local_atoms[ 4 * atom + 0 ] = x_coord_array[atom];
            data_from_local_atoms[ 4 * atom + 1 ] = y_coord_array[atom];
            data_from_local_atoms[ 4 * atom + 2 ] = z_coord_array[atom];
            data_from_local_atoms[ 4 * atom + 3 ] = moments_array[atom];
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
            mean_x += data_from_all_atoms[ 4 * atom + 0 ];
            mean_y += data_from_all_atoms[ 4 * atom + 1 ];
            mean_z += data_from_all_atoms[ 4 * atom + 2 ];
         }
         mean_x = mean_x / double (total_num_atoms);
         mean_y = mean_y / double (total_num_atoms);
         mean_z = mean_z / double (total_num_atoms);

         // output data to file atomistic_dipole_positions.txt
         if(vmpi::master){
            std::ofstream ofile;
            ofile.open("atomistic_dipole_positions.txt");
            for(uint64_t atom = 0; atom < total_num_atoms; atom++){
               ofile << data_from_all_atoms[ 4 * atom + 0 ] - mean_x << "\t" <<
                        data_from_all_atoms[ 4 * atom + 1 ] - mean_y << "\t" <<
                        data_from_all_atoms[ 4 * atom + 2 ] - mean_z << "\t" <<
                        data_from_all_atoms[ 4 * atom + 3 ] << "\t\n";
            }
            ofile.close();
         }

         // Every process wait here
         vmpi::barrier();

         // re-calculate counts and displacements for 3-vector data
         data_from_local_atoms.resize(3*local_num_atoms);
         if(vmpi::master) data_from_all_atoms.resize(3*total_num_atoms);
         vmpi::counts_and_displacements(data_from_local_atoms, data_from_all_atoms, counts, displacements); // set up counts and displacements

         // set initialised variable to true
         output_atomistic_field_initialised = true;

         return;

      }


      //------------------------------------------------------------------------
      // Function to output calculated atomistic dipole coords on root process
      //------------------------------------------------------------------------
      void output_atomistic_dipole_fields(){

         // check for initialised data structures
         if(!output_atomistic_field_initialised){
            std::cerr << "Programmer error : dipole::output_atomistic_dipole_fields() called before initialisation" << std::endl;
            err::vexit();
         }

         // get number of atoms (local CPU)
         const uint64_t local_num_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;

         // copy local field data to local array
         for(uint64_t atom = 0; atom < local_num_atoms; atom++){
            data_from_local_atoms[ 3 * atom + 0 ] = dipole::atom_dipolar_field_array_x[atom];
            data_from_local_atoms[ 3 * atom + 1 ] = dipole::atom_dipolar_field_array_y[atom];
            data_from_local_atoms[ 3 * atom + 2 ] = dipole::atom_dipolar_field_array_z[atom];
         }

         // collate global field data on master
         vmpi::fast_collate(data_from_local_atoms, data_from_all_atoms, counts, displacements);

         // output data to file atomistic_dipole_field.txt
         if(vmpi::master){
            std::ofstream ofile;
            ofile.open("atomistic_dipole_field.txt");
            for(int atom = 0; atom < total_num_atoms; atom++){
               ofile << data_from_all_atoms[ 3 * atom + 0 ] << "\t" <<
                        data_from_all_atoms[ 3 * atom + 1 ] << "\t" <<
                        data_from_all_atoms[ 3 * atom + 2 ] << "\t\n";
            }
            ofile.close();
         }

         return;

      }

} // end of namespace internal
} // end of namespace dipole
