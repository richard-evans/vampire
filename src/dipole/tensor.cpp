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
                                    const double cells_macro_cell_size,
                                    std::vector <int>& cells_local_cell_array,
                                    std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                                    std::vector <int>& cells_num_atoms_in_cell_global, /// number of atoms in each cell
                                    std::vector < std::vector <int> >& cells_index_atoms_array,
                                    const std::vector<double>& cells_volume_array,
                                    std::vector<double>& cells_pos_and_mom_array, // array to store positions and moment of cells
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,
                                    const std::vector<int>& atom_type_array,
                                    const std::vector<int>& atom_cell_id_array,
                                    const std::vector<double>& atom_coords_x, //atomic coordinates
                                    const std::vector<double>& atom_coords_y,
                                    const std::vector<double>& atom_coords_z,
                                    const int num_atoms){

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

            for(int lc=0; lc<dipole::internal::cells_num_cells; lc++){
               // resize arrays
               cells_atom_in_cell_coords_array_x[lc].resize(dipole::internal::cells_num_atoms_in_cell[lc]);
               cells_atom_in_cell_coords_array_y[lc].resize(dipole::internal::cells_num_atoms_in_cell[lc]);
               cells_atom_in_cell_coords_array_z[lc].resize(dipole::internal::cells_num_atoms_in_cell[lc]);
               cells_index_atoms_array[lc].resize(dipole::internal::cells_num_atoms_in_cell[lc]);
            }

            // Call parallelisation function
            // Exchange cells data
            dipole::internal::send_recv_cells_data(dipole::internal::proc_cell_index_array1D,
                                                   cells_atom_in_cell_coords_array_x,
                                                   cells_atom_in_cell_coords_array_y,
                                                   cells_atom_in_cell_coords_array_z,
                                                   cells_index_atoms_array,
                                                   dipole::internal::cells_pos_and_mom_array,
                                                   dipole::internal::cells_num_atoms_in_cell,
                                                   cells::cell_id_array,
                                                   dipole::internal::cells_local_cell_array,
                                                   dipole::internal::cells_num_local_cells,
                                                   dipole::internal::cells_num_cells);

            // Exchange atoms data
            dipole::internal::send_recv_atoms_data(dipole::internal::proc_cell_index_array1D,
                                                   cells::cell_id_array,
                                                   dipole::internal::cells_local_cell_array,
                                                   atom_pos_x,
                                                   atom_pos_y,
                                                   atom_pos_z,
                                                   dipole::internal::atom_type_array, // atomic moments (from dipole;:internal::atom_type_array)
                                                   cells_atom_in_cell_coords_array_x,
                                                   cells_atom_in_cell_coords_array_y,
                                                   cells_atom_in_cell_coords_array_z,
                                                   cells_index_atoms_array,
                                                   dipole::internal::cells_pos_and_mom_array,
                                                   dipole::internal::cells_num_atoms_in_cell,
                                                   dipole::internal::cells_num_local_cells,
                                                   dipole::internal::cells_num_cells,
                                                   cells_macro_cell_size);

            // Reorder data structure
            dipole::internal::sort_data(dipole::internal::proc_cell_index_array1D,
                                        cells::cell_id_array,
                                        cells_atom_in_cell_coords_array_x,
                                        cells_atom_in_cell_coords_array_y,
                                        cells_atom_in_cell_coords_array_z,
                                        cells_index_atoms_array,
                                        dipole::internal::cells_pos_and_mom_array,
                                        dipole::internal::cells_num_atoms_in_cell,
                                        dipole::internal::cells_num_local_cells,
                                        dipole::internal::cells_num_cells);

            // After transferring the data across cores, assign value dipole::internal::cells_num_atoms_in_cell[] from cells_num_atoms_in_cell_global[]
            for(unsigned int i=0; i<cells_num_atoms_in_cell_global.size(); i++){
               if(cells_num_atoms_in_cell_global[i]>0 && dipole::internal::cells_num_atoms_in_cell[i]==0){
                  dipole::internal::cells_num_atoms_in_cell[i] = cells_num_atoms_in_cell_global[i];
               }
            }

            // Clear memory
            cells_num_atoms_in_cell_global.clear();

            // Clear atom_pos_x,y,z
            atom_pos_x.clear();
            atom_pos_y.clear();
            atom_pos_z.clear();

         #endif

         // calculate matrix prefactors
         zlog << zTs() << "Precalculating rij matrix for dipole calculation using tensor solver... " << std::endl;
         std::cout     << "Precalculating rij matrix for dipole calculation using tensor solver"     << std::flush;

         // instantiate timer
         vutil::vtimer_t timer;

         // start timer
         timer.start();

         // loop over local cells
         for(int lc=0;lc<dipole::internal::cells_num_local_cells;lc++){

            // print out progress to screen
            if(fmod(ceil(lc),ceil(dipole::internal::cells_num_local_cells)/10.0) == 0) std::cout << "." << std::flush;

            // reference global cell ID
            //int i = dipole::internal::cells_local_cell_array[lc];
            int i = cells::cell_id_array[lc];
            if(dipole::internal::cells_num_atoms_in_cell[i]>0){

            	// Loop over all other cells to calculate contribution to local cell
               for(int j=0;j<dipole::internal::cells_num_cells;j++){

                  /*==========================================================*/
                  /* Calculation of inter part of dipolar tensor              */
                  /*==========================================================*/
                	if(i!=j && dipole::internal::cells_num_atoms_in_cell[j]>0){
                     // calculate inter term of dipolar tensor
                     compute_inter_tensor(cells_macro_cell_size,i,j,lc,cells_num_atoms_in_cell,cells_atom_in_cell_coords_array_x,cells_atom_in_cell_coords_array_y,cells_atom_in_cell_coords_array_z);
/*
                     // create temporary variable to store components of tensor
                    	double tmp_rij_inter_xx = 0.0;
                    	double tmp_rij_inter_xy = 0.0;
                    	double tmp_rij_inter_xz = 0.0;

                    	double tmp_rij_inter_yy = 0.0;
                    	double tmp_rij_inter_yz = 0.0;
                    	double tmp_rij_inter_zz = 0.0;

                     // Calculate distance vectors between cells
                     double rx = dipole::internal::cells_pos_and_mom_array[4*j+0] - dipole::internal::cells_pos_and_mom_array[4*i+0];
                     double ry = dipole::internal::cells_pos_and_mom_array[4*j+1] - dipole::internal::cells_pos_and_mom_array[4*i+1];
                     double rz = dipole::internal::cells_pos_and_mom_array[4*j+2] - dipole::internal::cells_pos_and_mom_array[4*i+2];

                    	double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance
                    	double rij_1 = 1.0/rij;

                     // If distance between macro-cells > cutoff nm => continuum approach (bare macro-cell method)
                     if( (rij_1)/cells_macro_cell_size > dipole::cutoff){
                        // define unitarian distance vectors
                     	const double ex = rx*rij;
                     	const double ey = ry*rij;
                     	const double ez = rz*rij;

                     	const double rij3 = (rij*rij*rij); // Angstroms

                        // calculate dipolar matrix for 6 entries because of symmetry
                     	dipole::internal::rij_inter_xx[lc][j] = ((3.0*ex*ex - 1.0)*rij3);
                     	dipole::internal::rij_inter_xy[lc][j] = ( 3.0*ex*ey      )*rij3 ;
                     	dipole::internal::rij_inter_xz[lc][j] = ( 3.0*ex*ez      )*rij3 ;

                     	dipole::internal::rij_inter_yy[lc][j] = ((3.0*ey*ey - 1.0)*rij3);
                     	dipole::internal::rij_inter_yz[lc][j] = ( 3.0*ey*ez      )*rij3 ;
                     	dipole::internal::rij_inter_zz[lc][j] = ((3.0*ez*ez - 1.0)*rij3);

                     }

                     //--------------------------------------------------------------------------
                     // If distance between macro-cells < cutoff ==> apply inter-intra method
                     //--------------------------------------------------------------------------
                     else if( (1.0/rij)/cells_macro_cell_size <= dipole::cutoff){

                        for(int pi=0; pi<dipole::internal::cells_num_atoms_in_cell[i]; pi++){

                           const double cix = cells_atom_in_cell_coords_array_x[i][pi];
                           const double ciy = cells_atom_in_cell_coords_array_y[i][pi];
                           const double ciz = cells_atom_in_cell_coords_array_z[i][pi];

                           for(int qj=0; qj<dipole::internal::cells_num_atoms_in_cell[j]; qj++){

                              const double rx = cells_atom_in_cell_coords_array_x[j][qj] - cix;
                              const double ry = cells_atom_in_cell_coords_array_y[j][qj] - ciy;
                              const double rz = cells_atom_in_cell_coords_array_z[j][qj] - ciz;

                              rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz);  //Reciprocal of the distance


                              const double ex = rx*rij;
                              const double ey = ry*rij;
                              const double ez = rz*rij;

                              const double rij3 = (rij*rij*rij); // Angstroms

                              tmp_rij_inter_xx += ((3.0*ex*ex - 1.0)*rij3);
                              tmp_rij_inter_xy += ((3.0*ex*ey      )*rij3);
                              tmp_rij_inter_xz += ((3.0*ex*ez      )*rij3);

                              tmp_rij_inter_yy += ((3.0*ey*ey - 1.0)*rij3);
                              tmp_rij_inter_yz += ((3.0*ey*ez      )*rij3);
                              tmp_rij_inter_zz += ((3.0*ez*ez - 1.0)*rij3);

                           }
                        }

                        dipole::internal::rij_inter_xx[lc][j] =  (tmp_rij_inter_xx);
                        dipole::internal::rij_inter_xy[lc][j] =  (tmp_rij_inter_xy);
                        dipole::internal::rij_inter_xz[lc][j] =  (tmp_rij_inter_xz);

                        dipole::internal::rij_inter_yy[lc][j] =  (tmp_rij_inter_yy);
                        dipole::internal::rij_inter_yz[lc][j] =  (tmp_rij_inter_yz);
                        dipole::internal::rij_inter_zz[lc][j] =  (tmp_rij_inter_zz);

                        // Normalisation by the number of atoms in the cell. This is required for the correct evaluation of the field in the update.cpp routine
                        dipole::internal::rij_inter_xx[lc][j] = dipole::internal::rij_inter_xx[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                        dipole::internal::rij_inter_xy[lc][j] = dipole::internal::rij_inter_xy[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                        dipole::internal::rij_inter_xz[lc][j] = dipole::internal::rij_inter_xz[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));

                        dipole::internal::rij_inter_yy[lc][j] = dipole::internal::rij_inter_yy[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                        dipole::internal::rij_inter_yz[lc][j] = dipole::internal::rij_inter_yz[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                        dipole::internal::rij_inter_zz[lc][j] = dipole::internal::rij_inter_zz[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));

                     }  // End of Inter part calculated atomicstically */

                  } // End of Inter part

                  /*==========================================================*/
                  /* Calculation of intra part of dipolar tensor              */
                  /*==========================================================*/
                  else if( i==j && dipole::internal::cells_num_atoms_in_cell[j]>0){
                     //Compute inter component of dipolar tensor
                     compute_intra_tensor(i,j,lc,cells_num_atoms_in_cell,cells_atom_in_cell_coords_array_x,cells_atom_in_cell_coords_array_y,cells_atom_in_cell_coords_array_z);
                     /*
                     // initialise temp vectors
                     double tmp_rij_intra_xx = 0.0;
                     double tmp_rij_intra_xy = 0.0;
                     double tmp_rij_intra_xz = 0.0;

                     double tmp_rij_intra_yy = 0.0;
                     double tmp_rij_intra_yz = 0.0;
                     double tmp_rij_intra_zz = 0.0;

                     const int mmax = dipole::internal::cells_num_atoms_in_cell[i];

                  	for(int pi=0; pi<mmax; pi++){

                        const double cix = cells_atom_in_cell_coords_array_x[i][pi];
                        const double ciy = cells_atom_in_cell_coords_array_y[i][pi];
                        const double ciz = cells_atom_in_cell_coords_array_z[i][pi];

                        // use double loops to avoid if pi != qj statement
                        for(int qj=0; qj<pi; qj++){

                           const double rx = cells_atom_in_cell_coords_array_x[j][qj] - cix; //cells_atom_in_cell_coords_array_x[i][pi];
                           const double ry = cells_atom_in_cell_coords_array_y[j][qj] - ciy; //cells_atom_in_cell_coords_array_y[i][pi];
                           const double rz = cells_atom_in_cell_coords_array_z[j][qj] - ciz; //cells_atom_in_cell_coords_array_z[i][pi];

                           const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance

                           const double ex = rx*rij;
                           const double ey = ry*rij;
                           const double ez = rz*rij;

                           const double rij3 = (rij*rij*rij); // Angstroms

                           tmp_rij_intra_xx += ((3.0*ex*ex - 1.0)*rij3);
                           tmp_rij_intra_xy += ((3.0*ex*ey      )*rij3);
                           tmp_rij_intra_xz += ((3.0*ex*ez      )*rij3);

                           tmp_rij_intra_yy += ((3.0*ey*ey - 1.0)*rij3);
                           tmp_rij_intra_yz += ((3.0*ey*ez      )*rij3);
                           tmp_rij_intra_zz += ((3.0*ez*ez - 1.0)*rij3);

                        }
                        for(int qj=pi+1; qj<mmax; qj++){

                           const double rx = cells_atom_in_cell_coords_array_x[j][qj] - cix; //cells_atom_in_cell_coords_array_x[i][pi];
                           const double ry = cells_atom_in_cell_coords_array_y[j][qj] - ciy; //cells_atom_in_cell_coords_array_y[i][pi];
                           const double rz = cells_atom_in_cell_coords_array_z[j][qj] - ciz; //cells_atom_in_cell_coords_array_z[i][pi];

                           const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance

                           const double ex = rx*rij;
                           const double ey = ry*rij;
                           const double ez = rz*rij;

                           const double rij3 = (rij*rij*rij); // Angstroms

                           tmp_rij_intra_xx += ((3.0*ex*ex - 1.0)*rij3);
                           tmp_rij_intra_xy += ((3.0*ex*ey      )*rij3);
                           tmp_rij_intra_xz += ((3.0*ex*ez      )*rij3);

                           tmp_rij_intra_yy += ((3.0*ey*ey - 1.0)*rij3);
                           tmp_rij_intra_yz += ((3.0*ey*ez      )*rij3);
                           tmp_rij_intra_zz += ((3.0*ez*ez - 1.0)*rij3);

                    		}
                   	}

                   	dipole::internal::rij_inter_xx[lc][i] =  (tmp_rij_intra_xx);
                   	dipole::internal::rij_inter_xy[lc][i] =  (tmp_rij_intra_xy);
                   	dipole::internal::rij_inter_xz[lc][i] =  (tmp_rij_intra_xz);

                   	dipole::internal::rij_inter_yy[lc][i] =  (tmp_rij_intra_yy);
                   	dipole::internal::rij_inter_yz[lc][i] =  (tmp_rij_intra_yz);
                   	dipole::internal::rij_inter_zz[lc][i] =  (tmp_rij_intra_zz);

                   	dipole::internal::rij_inter_xx[lc][i] = dipole::internal::rij_inter_xx[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                   	dipole::internal::rij_inter_xy[lc][i] = dipole::internal::rij_inter_xy[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                   	dipole::internal::rij_inter_xz[lc][i] = dipole::internal::rij_inter_xz[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));

                   	dipole::internal::rij_inter_yy[lc][i] = dipole::internal::rij_inter_yy[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                   	dipole::internal::rij_inter_yz[lc][i] = dipole::internal::rij_inter_yz[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                   	dipole::internal::rij_inter_zz[lc][i] = dipole::internal::rij_inter_zz[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));  */

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
