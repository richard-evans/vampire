//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C standard library headers
#include <cmath>
#include <cstdlib>

// Vampire headers
#include "cells.hpp" // needed for cells::macrocell_size but to be removed
#include "dipole.hpp"
#include "vio.hpp"
#include "vutil.hpp"
#include "atoms.hpp" // needed fro m spin array but to be removed

// dipole module headers
#include "internal.hpp"

namespace dipole{

   namespace internal{

      //------------------------------------------------------------------------
      // Function to compute the tensor row of cell celli: its interactions with
      // every cell, stored in rij_tensor_*[row]. Used for every local cell when
      // the tensors are stored, and at every field update with
      // dipole:on-the-fly-tensors (one row, reused for each local cell).
      //------------------------------------------------------------------------
      void compute_tensor_row(const int celli,
                              const int row,
                              const double real_cutoff,
                              const std::vector<int>& cells_num_atoms_in_cell_global,
                              const std::vector<double>& cells_pos_and_mom_array,
                              const std::vector<int>& list_of_atoms_with_cells,
                              const std::vector< std::vector<double> >& atoms_in_cells_array){

            // a cell with no magnetic atoms (probe cell) still gets a row; the tensor functions guard the empty cell
            //std::cout << i << '\t' << "interaction" << "\t"  << cells_num_cells << "\t" << cells_num_atoms_in_cell[i] <<  std::endl;

            // Loop over all other cells to calculate contribution to local cell
            for( int cellj = 0; cellj < cells_num_cells; cellj++ ){

               //--------------------------------------------------------------
               // Calculation of inter part of dipolar tensor
               //--------------------------------------------------------------
               if ( cells_num_atoms_in_cell_global[cellj] > 0 ){ // only calculate interaction if there are atoms in remote cell

                  // check if cell is not the same
                  if( celli != cellj ){
                     // calculate inter term of dipolar tensor
                     compute_inter_tensor(celli,
                                          cellj,
                                          row,
                                          real_cutoff,
                                          cells_num_atoms_in_cell_global,
                                          cells_pos_and_mom_array,
                                          list_of_atoms_with_cells,
                                          atoms_in_cells_array);

                  } // End of Inter part

                  //--------------------------------------------------------------
                  // Calculation of intra part of dipolar tensor
                  //--------------------------------------------------------------
                  else if( celli == cellj ){

                     //Compute inter component of dipolar tensor
                     compute_intra_tensor(celli,
                                          cellj,
                                          row,
                                          cells_num_atoms_in_cell_global,
                                          list_of_atoms_with_cells,
                                          atoms_in_cells_array);

                  }
                  // End of Intra part
               //   std::cout << i<< '\t' <<j << '\t' <<   dipole::internal::rij_tensor_xx[row][j]<< '\t' <<  dipole::internal::rij_tensor_xy[row][j]<< '\t' <<  dipole::internal::rij_tensor_xz[row][j]<< '\t' <<  dipole::internal::rij_tensor_yy[row][j] << '\t' <<  dipole::internal::rij_tensor_yz[row][j] << '\t' <<  dipole::internal::rij_tensor_zz[row][j] << std::endl;

               }
            }

         return;

      }

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
                                    const double cells_macro_cell_size,
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

            // all processors wait here before starting the timer
            vmpi::barrier();

            // get realspace cutoff distance. On a regular cell grid many cell pairs lie
            // exactly at the cutoff, and rounding of the cell centres (which depends on
            // the number of processors) would put each on either side of it; the
            // relative tolerance of 1e-6 treats them all atomistically
            const double real_cutoff = dipole::cutoff * cells::macro_cell_size * (1.0 + 1.0e-6);

            std::vector< std::vector<double> > atoms_in_cells_array; // 2D list of [cell][atom] for local cells needed for computing dipole tensor
            std::vector<int> list_of_atoms_with_cells; // list of cell IDs to enable parsing of atomistic data

            // distribute atomistic data to enable tensor dipole calculation
            initialise_atomistic_cell_data(cells_num_cells,
                                           cells_num_local_cells,
                                           real_cutoff,                     // cutoff range for dipole tensor construction (Angstroms)
                                           cells_num_atoms_in_cell,         // number of atoms in each cell (local CPU)
                                           cells_local_cell_array,          // numerical list of cells containing atoms on local processor
                                           cells_num_atoms_in_cell_global,  // number of atoms in each cell
                                           cells_pos_and_mom_array,         // array of positions and cell moments
                                           cells_index_atoms_array,         // 2D array of [cells][atomID]
                                           atom_coords_x,                   // input arrays of atom coordinates
                                           atom_coords_y,                   //
                                           atom_coords_z,                   //
                                           atoms::m_spin_array,             // input array of atom moments (Bohr magnetons)
                                           list_of_atoms_with_cells,        // 2D list of [cell][atom] for local cells needed for computing dipole tensor (output)
                                           atoms_in_cells_array             // list of cell IDs to enable parsing of atomistic data (output)
                                          );

         // After transferring the data across cores, assign value cells_num_atoms_in_cell[] from cells_num_atoms_in_cell_global[]
         for(unsigned int i=0; i<cells_num_atoms_in_cell_global.size(); i++){
            //if(cells_num_atoms_in_cell_global[i]>0 && cells_num_atoms_in_cell[i]==0){
               dipole::internal::cells_num_atoms_in_cell[i] = cells_num_atoms_in_cell_global[i];
            //}
         }

         // print informative message to user
         zlog << zTs() << "Precalculating rij matrix for dipole calculation using tensor solver... " << std::endl;
         std::cout     << "Precalculating rij matrix for dipole calculation using tensor solver"     << std::flush;

         // instantiate timer
         vutil::vtimer_t timer;

         // start timer
         timer.start();

         //--------------------------------------------------------------------------------------------
         // Compute the dipole tensor
         //--------------------------------------------------------------------------------------------

         // with dipole:on-the-fly-tensors the rows are computed at each field update instead
         // (dipole::internal::update_field), from the data kept here
         if(dipole::internal::on_the_fly_tensors){
            dipole::internal::otf_cutoff = real_cutoff;
            dipole::internal::otf_global_atoms_in_cell_count = cells_num_atoms_in_cell_global;
            dipole::internal::otf_cells_pos_and_mom_array = cells_pos_and_mom_array;
            dipole::internal::otf_list_of_cells_with_atoms.swap(list_of_atoms_with_cells);
            dipole::internal::otf_atoms_in_cells_array.swap(atoms_in_cells_array);
         }
         else{
         // loop over local cells
         for( int lc = 0; lc < cells_num_local_cells; lc++){

            // print out progress to screen
            if(fmod(ceil(lc),ceil(cells_num_local_cells)/10) == 0) std::cout << "." << std::flush;

               // get global cell ID of source cell and compute its tensor row
               compute_tensor_row(cells_local_cell_array[lc], lc, real_cutoff, cells_num_atoms_in_cell_global,
                                  cells_pos_and_mom_array, list_of_atoms_with_cells, atoms_in_cells_array);

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
