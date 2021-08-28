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

            // get realspace cutoff distance
            const double real_cutoff = dipole::cutoff * cells::macro_cell_size;

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

         // Assign updated value of cells_num_atoms_in_cell to dipole::dipole_cells_num_atoms_in_cell. It is needed to print the config file. The actual value cells::num_atoms_in_cell is not changed instead
         dipole::dipole_cells_num_atoms_in_cell = cells_num_atoms_in_cell;

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

         // loop over local cells
         for( int lc = 0; lc < cells_num_local_cells; lc++){

            // print out progress to screen
            if(fmod(ceil(lc),ceil(cells_num_local_cells)/10) == 0) std::cout << "." << std::flush;

            // get global cell ID of source cell
            int celli = cells_local_cell_array[lc];

            // check that local cell contains some local atoms (if not we don't need the tensor)
            if( cells_num_atoms_in_cell[celli] > 0 ){
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
                                             lc,
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
                                             lc,
                                             cells_num_atoms_in_cell_global,
                                             list_of_atoms_with_cells,
                                             atoms_in_cells_array);

                     }
                     // End of Intra part

                     // check for close to zero value tensors and round down to zero
                     if (dipole::internal::rij_tensor_xx[lc][cellj]*dipole::internal::rij_tensor_xx[lc][cellj] < 1e-15) dipole::internal::rij_tensor_xx[lc][cellj] = 0.0;
                     if (dipole::internal::rij_tensor_xy[lc][cellj]*dipole::internal::rij_tensor_xy[lc][cellj] < 1e-15) dipole::internal::rij_tensor_xy[lc][cellj] = 0.0;
                     if (dipole::internal::rij_tensor_xz[lc][cellj]*dipole::internal::rij_tensor_xz[lc][cellj] < 1e-15) dipole::internal::rij_tensor_xz[lc][cellj] = 0.0;
                     if (dipole::internal::rij_tensor_yy[lc][cellj]*dipole::internal::rij_tensor_yy[lc][cellj] < 1e-15) dipole::internal::rij_tensor_yy[lc][cellj] = 0.0;
                     if (dipole::internal::rij_tensor_yz[lc][cellj]*dipole::internal::rij_tensor_yz[lc][cellj] < 1e-15) dipole::internal::rij_tensor_yz[lc][cellj] = 0.0;
                     if (dipole::internal::rij_tensor_zz[lc][cellj]*dipole::internal::rij_tensor_zz[lc][cellj] < 1e-15) dipole::internal::rij_tensor_zz[lc][cellj] = 0.0;
               //   std::cout << i<< '\t' <<j << '\t' <<   dipole::internal::rij_tensor_xx[lc][j]<< '\t' <<  dipole::internal::rij_tensor_xy[lc][j]<< '\t' <<  dipole::internal::rij_tensor_xz[lc][j]<< '\t' <<  dipole::internal::rij_tensor_yy[lc][j] << '\t' <<  dipole::internal::rij_tensor_yz[lc][j] << '\t' <<  dipole::internal::rij_tensor_zz[lc][j] << std::endl;
               }
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
