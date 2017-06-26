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
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Vampire headers
#include "cells.hpp"
#include "dipole.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"

#include <fenv.h>
#include <signal.h>

// dipole module headers
#include "internal.hpp"

namespace dipole{

   //----------------------------------------------------------------------------
   // Function to initialize dipole module
   //----------------------------------------------------------------------------
   void initialize(const int cells_num_atoms_in_unit_cell,
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

                  const int num_atoms
				){

	//-------------------------------------------------------------------------------------
	// Check for dipole calculation enabled, if not do nothing
	//-------------------------------------------------------------------------------------
      if(!dipole::activated) return;

      // output informative message
      std::cout << "Initialising dipole field calculation" << std::endl;
		zlog << zTs() << "Initialising dipole field calculation" << std::endl;

		// check for prior initialisation
		if(dipole::internal::initialised){
      	zlog << zTs() << "Warning:  Dipole field calculation already initialised. Continuing." << std::endl;
      	return;
		}

		//-------------------------------------------------------------------------------------
		// Set simulation constants
		//-------------------------------------------------------------------------------------

      //-------------------------------------------------------------------------------------
      // Set const for functions
      //-------------------------------------------------------------------------------------

      dipole::internal::num_atoms                  = num_atoms;
      dipole::internal::atom_type_array            = atom_type_array;
      dipole::internal::atom_cell_id_array         = atom_cell_id_array;

      dipole::internal::cells_num_cells            = cells_num_cells;
      dipole::internal::cells_num_local_cells      = cells_num_local_cells;
      dipole::internal::cells_local_cell_array     = cells_local_cell_array;
      dipole::internal::cells_num_atoms_in_cell    = cells_num_atoms_in_cell;
      dipole::internal::cells_volume_array         = cells_volume_array;

      dipole::internal::cells_pos_and_mom_array    = cells_pos_and_mom_array;

		//-------------------------------------------------------------------------------------
		// Starting calculation of dipolar field
		//-------------------------------------------------------------------------------------

      // instantiate timer
      vutil::vtimer_t timer;

      // start timer
      timer.start();

       // Check memory requirements and print to screen
       zlog << zTs() << "Fast dipole field calculation has been enabled and requires " << double(dipole::internal::cells_num_cells)*double(dipole::internal::cells_num_local_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;
       std::cout << "Fast dipole field calculation has been enabled and requires " << double(dipole::internal::cells_num_cells)*double(dipole::internal::cells_num_local_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;

       // For MPI version, only add local atoms
       #ifdef MPICF
          int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
       #else
          int num_local_atoms = num_atoms;
       #endif


      #ifdef MPICF
         //------------------------------------------------------
         //        parallel section
         //------------------------------------------------------
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
         dipole::send_recv_cells_data(dipole::internal::proc_cell_index_array1D,
                                    cells_atom_in_cell_coords_array_x,
                                    cells_atom_in_cell_coords_array_y,
                                    cells_atom_in_cell_coords_array_z,
                                    cells_index_atoms_array,
                                    dipole::internal::cells_pos_and_mom_array,
                                    dipole::internal::cells_num_atoms_in_cell,
                                    cells::cell_id_array,
                                    dipole::internal::cells_local_cell_array,
                                    dipole::internal::cells_num_local_cells,
                                    dipole::internal::cells_num_cells
                                    );

         // Exchange atoms data
         dipole::send_recv_atoms_data(dipole::internal::proc_cell_index_array1D,
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
                                    cells_macro_cell_size
                                    );

         // Reorder data structure
      	sort_data(dipole::internal::proc_cell_index_array1D,
                  cells::cell_id_array,
                  cells_atom_in_cell_coords_array_x,
                  cells_atom_in_cell_coords_array_y,
                  cells_atom_in_cell_coords_array_z,
                  cells_index_atoms_array,
                  dipole::internal::cells_pos_and_mom_array,
                  dipole::internal::cells_num_atoms_in_cell,
                  dipole::internal::cells_num_local_cells,
                  dipole::internal::cells_num_cells
                  );

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

      // allocate arrays to store data [nloccell x ncells]
      for(int lc=0;lc<dipole::internal::cells_num_local_cells; lc++){

         dipole::internal::rij_inter_xx.push_back(std::vector<double>());
         dipole::internal::rij_inter_xx[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_inter_xy.push_back(std::vector<double>());
         dipole::internal::rij_inter_xy[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_inter_xz.push_back(std::vector<double>());
         dipole::internal::rij_inter_xz[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_inter_yy.push_back(std::vector<double>());
         dipole::internal::rij_inter_yy[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_inter_yz.push_back(std::vector<double>());
         dipole::internal::rij_inter_yz[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_inter_zz.push_back(std::vector<double>());
         dipole::internal::rij_inter_zz[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_intra_xx.push_back(std::vector<double>());
         dipole::internal::rij_intra_xx[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_intra_xy.push_back(std::vector<double>());
         dipole::internal::rij_intra_xy[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_intra_xz.push_back(std::vector<double>());
         dipole::internal::rij_intra_xz[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_intra_yy.push_back(std::vector<double>());
         dipole::internal::rij_intra_yy[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_intra_yz.push_back(std::vector<double>());
         dipole::internal::rij_intra_yz[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_intra_zz.push_back(std::vector<double>());
         dipole::internal::rij_intra_zz[lc].resize(dipole::internal::cells_num_cells,0.0);
         // resize B-field cells array
         dipole::cells_field_array_x.resize(dipole::internal::cells_num_cells,0.0);
         dipole::cells_field_array_y.resize(dipole::internal::cells_num_cells,0.0);
         dipole::cells_field_array_z.resize(dipole::internal::cells_num_cells,0.0);
         // resize mu_0*Hd-field cells array
         dipole::cells_mu0Hd_field_array_x.resize(dipole::internal::cells_num_cells,0.0);
         dipole::cells_mu0Hd_field_array_y.resize(dipole::internal::cells_num_cells,0.0);
         dipole::cells_mu0Hd_field_array_z.resize(dipole::internal::cells_num_cells,0.0);
      }

      zlog << zTs() << "Number of local cells for dipole calculation = " << dipole::internal::cells_num_local_cells << std::endl;
      zlog << zTs() << "Number of total cells for dipole calculation = " << dipole::internal::cells_num_cells << std::endl;

      // calculate matrix prefactors
      zlog << zTs() << "Precalculating rij matrix for dipole calculation... " << std::endl;
      std::cout     << "Precalculating rij matrix for dipole calculation"     << std::flush;

      //==========================================================
      //----------------------------------------------------------
      // Calculation of dipolar tensor
      //----------------------------------------------------------
      //==========================================================

      // loop over local cells
      for(int lc=0;lc<dipole::internal::cells_num_local_cells;lc++){

         // print out progress to screen
         if(lc % (dipole::internal::cells_num_local_cells/10) == 0) std::cout << "." << std::flush;

         // reference global cell ID
         //int i = dipole::internal::cells_local_cell_array[lc];
         int i = cells::cell_id_array[lc];
         if(dipole::internal::cells_num_atoms_in_cell[i]>0){

         	// Loop over all other cells to calculate contribution to local cell
            for(int j=0;j<dipole::internal::cells_num_cells;j++){

               /*==========================================================*/
               /* Calculation of intra part of dipolar tensor              */
               /*==========================================================*/
             	if(i!=j && dipole::internal::cells_num_atoms_in_cell[j]>0){
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
                  	double ex = rx*rij;
                  	double ey = ry*rij;
                  	double ez = rz*rij;

                  	double rij3 = (rij*rij*rij); // Angstroms
                     // calculate dipolar matrix for 6 entries because of symmetry
                  	dipole::internal::rij_inter_xx[lc][j] = ((3.0*ex*ex - 1.0)*rij3);
                  	dipole::internal::rij_inter_xy[lc][j] = ( 3.0*ex*ey      )*rij3 ;
                  	dipole::internal::rij_inter_xz[lc][j] = ( 3.0*ex*ez      )*rij3 ;

                  	dipole::internal::rij_inter_yy[lc][j] = ((3.0*ey*ey - 1.0)*rij3);
                  	dipole::internal::rij_inter_yz[lc][j] = ( 3.0*ey*ez      )*rij3 ;
                  	dipole::internal::rij_inter_zz[lc][j] = ((3.0*ez*ez - 1.0)*rij3);
                  }
                  // If distance between macro-cells < cutoff ==> apply inter-intra method
                  else if( (1.0/rij)/cells_macro_cell_size <= dipole::cutoff){
                     for(int pi=0; pi<dipole::internal::cells_num_atoms_in_cell[i]; pi++){
                        for(int qj=0; qj<dipole::internal::cells_num_atoms_in_cell[j]; qj++){

                           rx = cells_atom_in_cell_coords_array_x[j][qj] - cells_atom_in_cell_coords_array_x[i][pi];
                           ry = cells_atom_in_cell_coords_array_y[j][qj] - cells_atom_in_cell_coords_array_y[i][pi];
                           rz = cells_atom_in_cell_coords_array_z[j][qj] - cells_atom_in_cell_coords_array_z[i][pi];

                           rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz);  //Reciprocal of the distance
                           rij_1 = 1.0/rij;
                           // Check if there are not prolems with distances, otherwise print out error message
                           if( rij_1==0.0 ) {
                              fprintf(stderr,">>>>> (Inter)  Warning: atoms are overlapping in cells i=%d\tand j=%d\ton rank=%d\t<<<<<\n",i,j,vmpi::my_rank);
                              std::cout << ">>>>> (Inter) Warning: atoms are overlapping in cells i,j " << i << "\t" << j <<"<<<<<" << std::endl;
                              std::cout << "Cell and atomic coordinates used for calculation of Dipolar matrix" << std::endl;
                              std::cout << " xj= "  << dipole::internal::cells_pos_and_mom_array[4*j+0]  << " yj= "  << dipole::internal::cells_pos_and_mom_array[4*j+1]  << " zj= "  << dipole::internal::cells_pos_and_mom_array[4*j+2] << std::endl;
                              std::cout << " xi= "  << dipole::internal::cells_pos_and_mom_array[4*i+0]  << " yi= "  << dipole::internal::cells_pos_and_mom_array[4*i+1]  << " zi= "  << dipole::internal::cells_pos_and_mom_array[4*i+2] << std::endl;
                              std::cout << " pi= "<<pi << " xpi= " << cells_atom_in_cell_coords_array_x[i][pi] << " ypi= " << cells_atom_in_cell_coords_array_y[i][pi] << " zpi= " << cells_atom_in_cell_coords_array_z[i][pi] << std::endl;
                              std::cout << " qj= "<<qj << " xqj= " << cells_atom_in_cell_coords_array_x[j][qj] << " yqj= " << cells_atom_in_cell_coords_array_y[j][qj] << " zqj= " << cells_atom_in_cell_coords_array_z[j][qj] << std::endl;
                              std::cout << " rx= "  << rx << " ry= " << ry << " rz= " << rz << std::endl;
                              std::cout << " rij= " << 1.0/rij << " = " << (1.0/rij)/((cells_macro_cell_size*sqrt(3)+0.01)) << " macro-cells" << std::endl;
                              break;
                           }

                           const double ex = rx*rij;
                           const double ey = ry*rij;
                           const double ez = rz*rij;

                           double rij3 = (rij*rij*rij); // Angstroms

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
                  }  // End of Inter part calculated atomicstically
               } // End of Inter part

               /*==========================================================*/
               /* Calculation of intra part of dipolar tensor              */
               /*==========================================================*/
               else if( i==j && dipole::internal::cells_num_atoms_in_cell[j]>0){
                  // initialise temp vectors
                  double tmp_rij_intra_xx = 0.0;
                  double tmp_rij_intra_xy = 0.0;
                  double tmp_rij_intra_xz = 0.0;

                  double tmp_rij_intra_yy = 0.0;
                  double tmp_rij_intra_yz = 0.0;
                  double tmp_rij_intra_zz = 0.0;

               	for(int pi=0; pi<dipole::internal::cells_num_atoms_in_cell[i]; pi++){
                 		for(int qj=0; qj<dipole::internal::cells_num_atoms_in_cell[i]; qj++){
                  		if( pi!=qj ){

                   			double rx = cells_atom_in_cell_coords_array_x[j][qj] - cells_atom_in_cell_coords_array_x[i][pi];
                   			double ry = cells_atom_in_cell_coords_array_y[j][qj] - cells_atom_in_cell_coords_array_y[i][pi];
                   			double rz = cells_atom_in_cell_coords_array_z[j][qj] - cells_atom_in_cell_coords_array_z[i][pi];

                   			const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance

                   			if( 1.0/rij==0.0 ){
                              fprintf(stderr,">>>>> (Intra)  Warning: atoms are overlapping in cells i=%d\tand j=%d\ton rank=%d\t<<<<<\n",i,j,vmpi::my_rank);
                     			std::cout << ">>>>> (Intra)  Warning: atoms are overlapping in cells i=\t" << i << "\t and j=\t" << j <<"\t<<<<<" << std::endl;
                     			std::cout << "Cell and atomic coordinates used for calculation of Dipolar matrix" << std::endl;
                     			std::cout << " qj= "<<qj << " xqj= " << cells_atom_in_cell_coords_array_x[j][qj] << " yqj= " << cells_atom_in_cell_coords_array_y[j][qj] << " zqj= " << cells_atom_in_cell_coords_array_z[j][qj] << std::endl;
                     			std::cout << " pi= "<<pi << " xpi= " << cells_atom_in_cell_coords_array_x[i][pi] << " ypi= " << cells_atom_in_cell_coords_array_y[i][pi] << " zpi= " << cells_atom_in_cell_coords_array_z[i][pi] << std::endl;
                              std::cout << " xj= "  << dipole::internal::cells_pos_and_mom_array[4*j+0]  << " yj= "  << dipole::internal::cells_pos_and_mom_array[4*j+1]  << " zj= "  << dipole::internal::cells_pos_and_mom_array[4*j+2] << std::endl;
                              std::cout << " xi= "  << dipole::internal::cells_pos_and_mom_array[4*i+0]  << " yi= "  << dipole::internal::cells_pos_and_mom_array[4*i+1]  << " zi= "  << dipole::internal::cells_pos_and_mom_array[4*i+2] << std::endl;
                     			std::cout << " rx= "  << rx << " ry= " << ry << " rz= " << rz << std::endl;
                     			std::cout << " rij= " << 1.0/rij << " = " << (1.0/rij)/((cells_macro_cell_size*sqrt(3)+0.01)) << " macro-cells" << std::endl;
                              break;
                   			}

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
                	}

                	dipole::internal::rij_intra_xx[lc][i] =  (tmp_rij_intra_xx);
                	dipole::internal::rij_intra_xy[lc][i] =  (tmp_rij_intra_xy);
                	dipole::internal::rij_intra_xz[lc][i] =  (tmp_rij_intra_xz);

                	dipole::internal::rij_intra_yy[lc][i] =  (tmp_rij_intra_yy);
                	dipole::internal::rij_intra_yz[lc][i] =  (tmp_rij_intra_yz);
                	dipole::internal::rij_intra_zz[lc][i] =  (tmp_rij_intra_zz);

                	dipole::internal::rij_intra_xx[lc][i] = dipole::internal::rij_intra_xx[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                	dipole::internal::rij_intra_xy[lc][i] = dipole::internal::rij_intra_xy[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                	dipole::internal::rij_intra_xz[lc][i] = dipole::internal::rij_intra_xz[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));

                	dipole::internal::rij_intra_yy[lc][i] = dipole::internal::rij_intra_yy[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                	dipole::internal::rij_intra_yz[lc][i] = dipole::internal::rij_intra_yz[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                	dipole::internal::rij_intra_zz[lc][i] = dipole::internal::rij_intra_zz[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
               } // End of Intra part
            }
			}
		}

      cells::num_cells = dipole::internal::cells_num_cells;
      cells::num_atoms_in_cell = dipole::internal::cells_num_atoms_in_cell;

      // hold parallel calculation until all processors have completed the dipole calculation
      vmpi::barrier();

      // stop timer
      timer.stop();

      std::cout << "Done! [ " << timer.elapsed_time() << " s ]" << std::endl;
      zlog << zTs() << "Precalculation of rij matrix for dipole calculation complete. Time taken: " << timer.elapsed_time() << " s"<< std::endl;

      // Set initialised flag
      dipole::internal::initialised=true;

      // start timer
      timer.start();

      // now calculate fields
      dipole::calculate_field();

      // hold parallel calculation until all processors have completed the update
      vmpi::barrier();

      // stop timer
      timer.stop();

      zlog << zTs() << "Time required for dipole update: " << timer.elapsed_time() << " s." << std::endl;

	   //-------------------------------------------------------//
	   //------- CPUs OUTPUT Dij on different fiels ------------//
	   //-------------------------------------------------------//

      // Output informative message to log file
      zlog << zTs() << "Outputting dipole matrix " << std::endl;

      // Output Demag tensor only if first step of simulation since depending only on shape
      if(sim::time == 0){

         int num_atoms_magnetic = 0.0;   // Initialise tot num of magnetic atoms
         // Calculate number of magnetic atoms
         for(int lc=0; lc<dipole::internal::cells_num_local_cells; lc++){
            int i = cells::cell_id_array[lc];
            num_atoms_magnetic += dipole::internal::cells_num_atoms_in_cell[i];
         }

         // Define and initialise Demag factor N tensor components
         double Nxx = 0.0;
         double Nxy = 0.0;
         double Nxz = 0.0;
         double Nyy = 0.0;
         double Nyz = 0.0;
         double Nzz = 0.0;

         // Every cpus print to check dipolar martrix inter term
         for(int lc=0; lc<dipole::internal::cells_num_local_cells; lc++){
            int i = cells::cell_id_array[lc];
            if(dipole::internal::cells_num_atoms_in_cell[i]>0){

               for(unsigned int j=0; j<dipole::internal::rij_inter_xx[lc].size(); j++){
                  if(dipole::internal::cells_num_atoms_in_cell[j]>0){

                     // To obtain dipolar matrix free of units, multiply tensor by "factor"
                     const double Vatomic = dipole::internal::cells_volume_array[j]/double(dipole::internal::cells_num_atoms_in_cell[j]);
                     const double factor = Vatomic*double(dipole::internal::cells_num_atoms_in_cell[j]) * double(dipole::internal::cells_num_atoms_in_cell[i]);
                     // Sum over dipolar tensor to obtain total tensor
                     Nxx += factor*(dipole::internal::rij_intra_xx[lc][j]+dipole::internal::rij_inter_xx[lc][j]);
                     Nxy += factor*(dipole::internal::rij_intra_xy[lc][j]+dipole::internal::rij_inter_xy[lc][j]);
                     Nxz += factor*(dipole::internal::rij_intra_xz[lc][j]+dipole::internal::rij_inter_xz[lc][j]);
                     Nyy += factor*(dipole::internal::rij_intra_yy[lc][j]+dipole::internal::rij_inter_yy[lc][j]);
                     Nyz += factor*(dipole::internal::rij_intra_yz[lc][j]+dipole::internal::rij_inter_yz[lc][j]);
                     Nzz += factor*(dipole::internal::rij_intra_zz[lc][j]+dipole::internal::rij_inter_zz[lc][j]);
                  }
               }
            }
         }

         // Reduce values of num of magnetic atoms and demag factors on all procs
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &num_atoms_magnetic, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &Nxx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &Nxy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &Nxz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &Nyy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &Nyz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &Nzz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif

         // Compute demag factor tensor from dipolar matrix adding self term
         Nxx = ((Nxx /num_atoms_magnetic)-4.0*M_PI/3.0)/(-4.0*M_PI);
         Nxy =  (Nxy /num_atoms_magnetic              )/(-4.0*M_PI);
         Nxz =  (Nxz /num_atoms_magnetic              )/(-4.0*M_PI);
         Nyy = ((Nyy /num_atoms_magnetic)-4.0*M_PI/3.0)/(-4.0*M_PI);
         Nyz =  (Nyz /num_atoms_magnetic              )/(-4.0*M_PI);
         Nzz = ((Nzz /num_atoms_magnetic)-4.0*M_PI/3.0)/(-4.0*M_PI);

         // Write demag factor to log file zlog << zTs() <<
         zlog << zTs() << "Demagnetisation tensor in format Nxx\t\tNxy\t\tNxz\t\tNyx\t\tNyy\tNyz\t\tNzx\t\tNzy\t\tNzz :\n";
         zlog << zTs() << Nxx << "\t" << Nxy << "\t" << Nxz << "\t";
         zlog <<          Nxy << "\t" << Nyy << "\t" << Nyz << "\t";
         zlog <<          Nxz << "\t" << Nyz << "\t" << Nzz << "\n";

      } // close if loop for sim::time == 0
	   //--------------------------------------------------/
      // End of outptu dipolar tensor
	   //--------------------------------------------------/

    	return;
    }
} // end of dipole namespace
