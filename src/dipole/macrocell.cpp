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
      // Function to initialise dipole tensors with bare macrocell solver
      //------------------------------------------------------------------------
      void initialize_macrocell_solver(const int cells_num_atoms_in_unit_cell,
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
                                                   cells_pos_and_mom_array,
                                                   cells_num_atoms_in_cell,
                                                   cells_num_local_cells,
                                                   cells_num_cells,
                                                   cells_macro_cell_size);

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
         zlog << zTs() << "Precalculating rij matrix for dipole calculation using macrocell solver... " << std::endl;
         std::cout     << "Precalculating rij matrix for dipole calculation using macrocell solver"     << std::flush;

         // instantiate timer
         vutil::vtimer_t timer;

         // start timer
         timer.start();

         // loop over local cells
         for(int lc=0;lc<cells_num_local_cells;lc++){

            // print out progress to screen
            //if(lc % (dipole::internal::cells_num_local_cells/10) == 0) std::cout << "." << std::flush;
            if(fmod(ceil(lc),ceil(cells_num_local_cells)/10.0) == 0) std::cout << "." << std::flush;

            // reference global cell ID
            //int i = dipole::internal::cells_local_cell_array[lc];
            int i = cells::cell_id_array[lc];

            // check that the cell constains at least one atom
            if(cells_num_atoms_in_cell[i]>0){

            	// Loop over all other cells to calculate contribution to local cell
               for(int j=0;j<cells_num_cells;j++){

                  /*==========================================================*/
                  /* Calculation of inter part of dipolar tensor              */
                  /*==========================================================*/
                	if(i!=j && cells_num_atoms_in_cell[j]>0){

                     // Calculate distance vectors between cells
                     double rx = cells_pos_and_mom_array[4*j+0] - cells_pos_and_mom_array[4*i+0];
                     double ry = cells_pos_and_mom_array[4*j+1] - cells_pos_and_mom_array[4*i+1];
                     double rz = cells_pos_and_mom_array[4*j+2] - cells_pos_and_mom_array[4*i+2];

                    	double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance

                     // define unitarian distance vectors
                  	const double ex = rx*rij;
                  	const double ey = ry*rij;
                  	const double ez = rz*rij;

                  	const double rij3 = (rij*rij*rij); // Angstroms

                     // calculate dipolar matrix for 6 entries because of symmetry
                  	dipole::internal::rij_tensor_xx[lc][j] = ((3.0*ex*ex - 1.0)*rij3);
                  	dipole::internal::rij_tensor_xy[lc][j] = ( 3.0*ex*ey      )*rij3 ;
                  	dipole::internal::rij_tensor_xz[lc][j] = ( 3.0*ex*ez      )*rij3 ;

                  	dipole::internal::rij_tensor_yy[lc][j] = ((3.0*ey*ey - 1.0)*rij3);
                  	dipole::internal::rij_tensor_yz[lc][j] = ( 3.0*ey*ez      )*rij3 ;
                  	dipole::internal::rij_tensor_zz[lc][j] = ((3.0*ez*ez - 1.0)*rij3);
                  //   std::cout << dipole::internal::rij_tensor_xx[lc][j] << std::endl;
                  }

                  /*==========================================================*/
                  /* Calculation of intra part of dipolar tensor              */
                  /*==========================================================*/

                  else if( i==j && dipole::internal::cells_num_atoms_in_cell[j]>0){

                     double tmp_rij_intra_xx = 0.0;
                     double tmp_rij_intra_xy = 0.0;
                     double tmp_rij_intra_xz = 0.0;

                     double tmp_rij_intra_yy = 0.0;
                     double tmp_rij_intra_yz = 0.0;
                     double tmp_rij_intra_zz = 0.0;

                     // const int mmax = cells_num_atoms_in_cell[i];
                     //
                     // for(int pi=0; pi<mmax; pi++){
                     //
                     //      const double cix = cells_atom_in_cell_coords_array_x[i][pi];
                     //      const double ciy = cells_atom_in_cell_coords_array_y[i][pi];
                     //      const double ciz = cells_atom_in_cell_coords_array_z[i][pi];
                     //
                     //      // use double loops to avoid if pi != qj statement
                     //      for(int qj=0; qj<pi; qj++){
                     //
                     //         const double rx = cells_atom_in_cell_coords_array_x[j][qj] - cix; //cells_atom_in_cell_coords_array_x[i][pi];
                     //         const double ry = cells_atom_in_cell_coords_array_y[j][qj] - ciy; //cells_atom_in_cell_coords_array_y[i][pi];
                     //         const double rz = cells_atom_in_cell_coords_array_z[j][qj] - ciz; //cells_atom_in_cell_coords_array_z[i][pi];
                     //
                     //         const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance
                     //
                     //         const double ex = rx*rij;
                     //         const double ey = ry*rij;
                     //         const double ez = rz*rij;
                     //
                     //         const double rij3 = (rij*rij*rij); // Angstroms
                     //
                     //         tmp_rij_intra_xx += ((3.0*ex*ex - 1.0)*rij3);
                     //         tmp_rij_intra_xy += ((3.0*ex*ey      )*rij3);
                     //         tmp_rij_intra_xz += ((3.0*ex*ez      )*rij3);
                     //
                     //         tmp_rij_intra_yy += ((3.0*ey*ey - 1.0)*rij3);
                     //         tmp_rij_intra_yz += ((3.0*ey*ez      )*rij3);
                     //         tmp_rij_intra_zz += ((3.0*ez*ez - 1.0)*rij3);
                     //
                     //      }
                     //      for(int qj=pi+1; qj<mmax; qj++){
                     //
                     //         const double rx = cells_atom_in_cell_coords_array_x[j][qj] - cix; //cells_atom_in_cell_coords_array_x[i][pi];
                     //         const double ry = cells_atom_in_cell_coords_array_y[j][qj] - ciy; //cells_atom_in_cell_coords_array_y[i][pi];
                     //         const double rz = cells_atom_in_cell_coords_array_z[j][qj] - ciz; //cells_atom_in_cell_coords_array_z[i][pi];
                     //
                     //         const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance
                     //
                     //         const double ex = rx*rij;
                     //         const double ey = ry*rij;
                     //         const double ez = rz*rij;
                     //
                     //         const double rij3 = (rij*rij*rij); // Angstroms
                     //
                     //         tmp_rij_intra_xx += ((3.0*ex*ex - 1.0)*rij3);
                     //         tmp_rij_intra_xy += ((3.0*ex*ey      )*rij3);
                     //         tmp_rij_intra_xz += ((3.0*ex*ez      )*rij3);
                     //
                     //         tmp_rij_intra_yy += ((3.0*ey*ey - 1.0)*rij3);
                     //         tmp_rij_intra_yz += ((3.0*ey*ez      )*rij3);
                     //         tmp_rij_intra_zz += ((3.0*ez*ez - 1.0)*rij3);
                     //
                     //    }
                     // }

                     dipole::internal::rij_tensor_xx[lc][i] =  (tmp_rij_intra_xx);
                     dipole::internal::rij_tensor_xy[lc][i] =  (tmp_rij_intra_xy);
                     dipole::internal::rij_tensor_xz[lc][i] =  (tmp_rij_intra_xz);

                     dipole::internal::rij_tensor_yy[lc][i] =  (tmp_rij_intra_yy);
                     dipole::internal::rij_tensor_yz[lc][i] =  (tmp_rij_intra_yz);
                     dipole::internal::rij_tensor_zz[lc][i] =  (tmp_rij_intra_zz);

                     dipole::internal::rij_tensor_xx[lc][i] = dipole::internal::rij_tensor_xx[lc][i]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                     dipole::internal::rij_tensor_xy[lc][i] = dipole::internal::rij_tensor_xy[lc][i]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                     dipole::internal::rij_tensor_xz[lc][i] = dipole::internal::rij_tensor_xz[lc][i]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));

                     dipole::internal::rij_tensor_yy[lc][i] = dipole::internal::rij_tensor_yy[lc][i]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                     dipole::internal::rij_tensor_yz[lc][i] = dipole::internal::rij_tensor_yz[lc][i]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                     dipole::internal::rij_tensor_zz[lc][i] = dipole::internal::rij_tensor_zz[lc][i]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));


                  } // End of Intra part
               }
   			}
   		}
         //std::cin.get();
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
