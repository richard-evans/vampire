//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <iostream>

// C library headers
#include <fenv.h>
#include <signal.h>

// Vampire headers
#include "cells.hpp" // needed for cells::cell_id_array but to be removed
#include "dipole.hpp"
#include "vio.hpp"
#include "vutil.hpp"

// dipole module headers
#include "internal.hpp"

namespace dipole{

   namespace internal{

      //------------------------------------------------------------------------
      // Function to initialise dipole tensors with default scheme.
      //
      // The tensors between local cells with the cutoff range are calculated
      // explictly from the atomistic coordinates. Longer range tensors assume
      // the dipole-dipole form.
      //------------------------------------------------------------------------
      void initialize_tensor_solver(const int cells_num_atoms_in_unit_cell,
                                    int cells_num_cells, /// number of macrocells
                                    int cells_num_local_cells, /// number of local macrocells
                                    const double cells_macro_cell_size_x,
                                    const double cells_macro_cell_size_y,
                                    const double cells_macro_cell_size_z,
                                    std::vector <int>& cells_local_cell_array,
                                    std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                                    std::vector <int>& cells_num_atoms_in_cell_global, /// number of atoms in each cell
                                    std::vector < std::vector <int> >& cells_index_atoms_array,
                                    std::vector<double>& cells_volume_array,
                                    std::vector<double>& cells_pos_and_mom_array, // array to store positions and moment of cells
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,
                                    std::vector<int>& atom_type_array,
                                    std::vector<int>& atom_cell_id_array,
                                    std::vector<double>& atom_coords_x, //atomic coordinates
                                    std::vector<double>& atom_coords_y,
                                    std::vector<double>& atom_coords_z,
                                    int num_atoms){

         //------------------------------------------------------
         // Collate atom coordinates for local cells
         //------------------------------------------------------
         #ifdef MPICF

            const int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;

            std::vector<double> atom_pos_x(num_local_atoms,0.0);
            std::vector<double> atom_pos_y(num_local_atoms,0.0);
            std::vector<double> atom_pos_z(num_local_atoms,0.0);

            for(int atom=0; atom<num_local_atoms; atom++){
               atom_pos_x[atom]=atom_coords_x[atom];
               atom_pos_y[atom]=atom_coords_y[atom];
               atom_pos_z[atom]=atom_coords_z[atom];
            }

            for(int lc=0; lc<cells_num_cells; lc++){
               // resize arrays
               cells_atom_in_cell_coords_array_x[lc].resize(cells_num_atoms_in_cell[lc]);
               cells_atom_in_cell_coords_array_y[lc].resize(cells_num_atoms_in_cell[lc]);
               cells_atom_in_cell_coords_array_z[lc].resize(cells_num_atoms_in_cell[lc]);
               cells_index_atoms_array[lc].resize(cells_num_atoms_in_cell[lc]);
            }

            // Call parallelisation function
            // Exchange cells data
            dipole::internal::send_recv_cells_data(dipole::internal::proc_cell_index_array1D,
                                                   cells_atom_in_cell_coords_array_x,
                                                   cells_atom_in_cell_coords_array_y,
                                                   cells_atom_in_cell_coords_array_z,
                                                   cells_index_atoms_array,
                                                   cells_pos_and_mom_array,
                                                   cells_num_atoms_in_cell,
                                                   cells::cell_id_array,
                                                   cells_local_cell_array,
                                                   cells_num_local_cells,
                                                   cells_num_cells);

            // Exchange atoms data
            dipole::internal::send_recv_atoms_data(dipole::internal::proc_cell_index_array1D,
                                                   cells::cell_id_array,
                                                   cells_local_cell_array,
                                                   atom_pos_x,
                                                   atom_pos_y,
                                                   atom_pos_z,
                                                   atom_type_array, // atomic moments (from dipole;:internal::atom_type_array)
                                                   cells_atom_in_cell_coords_array_x,
                                                   cells_atom_in_cell_coords_array_y,
                                                   cells_atom_in_cell_coords_array_z,
                                                   cells_index_atoms_array,
                                                   dipole::internal::cells_pos_and_mom_array,
                                                   dipole::internal::cells_num_atoms_in_cell,
                                                   dipole::internal::cells_num_local_cells,
                                                   dipole::internal::cells_num_cells,
                                                   cells_macro_cell_size_x,cells_macro_cell_size_y,cells_macro_cell_size_z );

            // Reorder data structure
            dipole::internal::sort_data(dipole::internal::proc_cell_index_array1D,
                                        cells::cell_id_array,
                                        cells_atom_in_cell_coords_array_x,
                                        cells_atom_in_cell_coords_array_y,
                                        cells_atom_in_cell_coords_array_z,
                                        cells_index_atoms_array,
                                        cells_pos_and_mom_array,
                                        cells_num_atoms_in_cell,
                                        cells_num_local_cells,
                                        cells_num_cells);

            // After transferring the data across cores, assign value cells_num_atoms_in_cell[] from cells_num_atoms_in_cell_global[]
            for(unsigned int i=0; i<cells_num_atoms_in_cell_global.size(); i++){
               if(cells_num_atoms_in_cell_global[i]>0 && cells_num_atoms_in_cell[i]==0){
                  cells_num_atoms_in_cell[i] = cells_num_atoms_in_cell_global[i];
               }
            }

            // Clear memory
            cells_num_atoms_in_cell_global.clear();

            // Clear atom_pos_x,y,z
            atom_pos_x.clear();
            atom_pos_y.clear();
            atom_pos_z.clear();

         #endif

         // Assign updated value of cells_num_atoms_in_cell to dipole::dipole_cells_num_atoms_in_cell. It is needed to print the config file. The actual value cells::num_atoms_in_cell is not changed instead
         dipole::dipole_cells_num_atoms_in_cell=cells_num_atoms_in_cell;

         // calculate matrix prefactors
         zlog << zTs() << "Precalculating rij matrix for dipole calculation using tensor solver... " << std::endl;
         std::cout     << "Precalculating rij matrix for dipole calculation using tensor solver"     << std::flush;

         // instantiate timer
         vutil::vtimer_t timer;

         // start timer
         timer.start();

         // loop over local cells
         for(int lc=0;lc<cells_num_local_cells;lc++){

            // print out progress to screen
            if(fmod(ceil(lc),ceil(cells_num_local_cells)/10.0) == 0) std::cout << "." << std::flush;

            // reference global cell ID
            //int i = cells_local_cell_array[lc];
            int i = cells::cell_id_array[lc];
            if(cells_num_atoms_in_cell[i]>0){

            	// Loop over all other cells to calculate contribution to local cell
               for(int j=0;j<cells_num_cells;j++){

                  /*==========================================================*/
                  /* Calculation of inter part of dipolar tensor              */
                  /*==========================================================*/
                	if(i!=j && cells_num_atoms_in_cell[j]>0){

                     // calculate inter term of dipolar tensor
                     compute_inter_tensor(cells_macro_cell_size_x,cells_macro_cell_size_y,cells_macro_cell_size_z, i,j,lc,cells_num_atoms_in_cell,cells_atom_in_cell_coords_array_x,cells_atom_in_cell_coords_array_y,cells_atom_in_cell_coords_array_z);

                  } // End of Inter part

                  /*==========================================================*/
                  /* Calculation of intra part of dipolar tensor              */
                  /*==========================================================*/
                  else if( i==j && cells_num_atoms_in_cell[j]>0){

                     //Compute inter component of dipolar tensor
                     compute_intra_tensor(i,j,lc,cells_num_atoms_in_cell,cells_atom_in_cell_coords_array_x,cells_atom_in_cell_coords_array_y,cells_atom_in_cell_coords_array_z);

                  } // End of Intra part
               }
   			}
   		}

         // hold parallel calculation until all processors have completed the dipole calculation
         vmpi::barrier();

         // stop timer
         timer.stop();

         std::cout << "done! [ " << timer.elapsed_time() << " s ]" << std::endl;
         zlog << zTs() << "Precalculation of rij matrix for dipole calculation complete. Time taken: " << timer.elapsed_time() << " s"<< std::endl;

         return;

      }

   } // end of namespace internal

} // end of namespace dipole
