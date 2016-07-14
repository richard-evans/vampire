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
#include <vector>
#include <string>

// Vampire headers
#include "cells.hpp"
#include "dipole.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// dipole module headers
#include "internal.hpp"

namespace dipole{

   //----------------------------------------------------------------------------
   // Function to initialize dipole module
   //----------------------------------------------------------------------------
   void initialize(const int cells_num_atoms_in_unit_cell,
                   const int cells_num_cells, /// number of macrocells
                   const int cells_num_local_cells, /// number of local macrocells
                   const double cells_macro_cell_size
                   const std::vector <int>& cells_local_cell_array,
                   const std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                   const std::vector < std::vector <int> >& cells_index_atoms_array,

                   const std::vector<double>& cells_volume_array,
                   const std::vector<double>& cells_cell_coords_array_x, /// arrays to store cells positions
                   const std::vector<double>& cells_cell_coords_array_y,
                   const std::vector<double>& cells_cell_coords_array_z,
                   const std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                   const std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                   const std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,
                   const std::vector<double>& cells_mag_array_x, /// arrays to store cells magnetisation
                   const std::vector<double>& cells_mag_array_y,
                   const std::vector<double>& cells_mag_array_z,
                   const std::vector<double>& cells_field_array_x, /// arrays to store cells field
                   const std::vector<double>& cells_field_array_y,
                   const std::vector<double>& cells_field_array_z,

                   const std::vector<int>& atom_type_array,
                   const std::vector<int>& atom_cell_array,
                   const int num_atoms,

                   const std::vector<double>& atom_dipolar_field_array_x, /// arrays to store atoms dipolar field
                   const std::vector<double>& atom_dipolar_field_array_y,
                   const std::vector<double>& atom_dipolar_field_array_z
				){

	//-------------------------------------------------------------------------------------
	// Check for dipole calculation enabled, if not do nothing
	//-------------------------------------------------------------------------------------
		if(!dipole::internal::enabled) return;

		// output informative message
		zlog << zTs() << "Initialising demagnetisation field calculation" << std::endl;

		// check for prior initialisation
		if(dipole::internal::initialised){
      		zlog << zTs() << "Warning:  Demagnetisation field calculation already initialised. Continuing." << std::endl;
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
      dipole::internal::atom_cell_array            = atom_cell_array;
      dipole::internal::atom_dipolar_field_array_x = atom_dipolar_field_array_x;
      dipole::internal::atom_dipolar_field_array_y = atom_dipolar_field_array_y;
      dipole::internal::atom_dipolar_field_array_z = atom_dipolar_field_array_z;

      dipole::internal::cells_mag_array_x          = cells_mag_array_x;
      dipole::internal::cells_mag_array_y          = cells_mag_array_y;
      dipole::internal::cells_mag_array_z          = cells_mag_array_z;
      dipole::internal::cells_field_array_x        = cells_field_array_x;
      dipole::internal::cells_field_array_y        = cells_field_array_y;
      dipole::internal::cells_field_array_z        = cells_field_array_z;
      dipole::internal::cells_volume_array         = cells_volume_array;

		//-------------------------------------------------------------------------------------
		// Starting calculation of dipolar field
		//-------------------------------------------------------------------------------------
       // timing function
       #ifdef MPICF
          double t1 = MPI_Wtime();
       #else
          time_t t1;
          t1 = time (NULL);
       #endif
         // Check memory requirements and print to screen
         zlog << zTs() << "Fast demagnetisation field calculation has been enabled and requires " << double(cells_num_cells)*double(cells_num_local_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;
         std::cout << "Fast demagnetisation field calculation has been enabled and requires " << double(cells_num_cells)*double(cells_num_local_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;

       // allocate arrays to store data [nloccell x ncells]
      for(int lc=0;lc<cells_num_local_cells; lc++){

         dipole::internal::rij_inter_xx.push_back(std::vector<double>());
         dipole::internal::rij_inter_xx[lc].resize(cells_num_cells,0.0);

         dipole::internal::rij_inter_xy.push_back(std::vector<double>());
         dipole::internal::rij_inter_xy[lc].resize(cells_num_cells,0.0);

         dipole::internal::rij_inter_xz.push_back(std::vector<double>());
         dipole::internal::rij_inter_xz[lc].resize(cells_num_cells,0.0);

         dipole::internal::rij_inter_yy.push_back(std::vector<double>());
         dipole::internal::rij_inter_yy[lc].resize(cells_num_cells,0.0);

         dipole::internal::rij_inter_yz.push_back(std::vector<double>());
         dipole::internal::rij_inter_yz[lc].resize(cells_num_cells,0.0);

         dipole::internal::rij_inter_zz.push_back(std::vector<double>());
         dipole::internal::rij_inter_zz[lc].resize(cells_num_cells,0.0);

         dipole::internal::rij_intra_xx.push_back(std::vector<double>());
         dipole::internal::rij_intra_xx[lc].resize(cells_num_cells,0.0);

         dipole::internal::rij_intra_xy.push_back(std::vector<double>());
         dipole::internal::rij_intra_xy[lc].resize(cells_num_cells,0.0);

         dipole::internal::rij_intra_xz.push_back(std::vector<double>());
         dipole::internal::rij_intra_xz[lc].resize(cells_num_cells,0.0);

         dipole::internal::rij_intra_yy.push_back(std::vector<double>());
         dipole::internal::rij_intra_yy[lc].resize(cells_num_cells,0.0);

         dipole::internal::rij_intra_yz.push_back(std::vector<double>());
         dipole::internal::rij_intra_yz[lc].resize(cells_num_cells,0.0);

         dipole::internal::rij_intra_zz.push_back(std::vector<double>());
         dipole::internal::rij_intra_zz[lc].resize(cells_num_cells,0.0);
      }

      // calculate matrix prefactors
      zlog << zTs() << "Precalculating rij matrix for demag calculation... " << std::endl;


      std::cout<< "Number of local cells= "<<cells_num_local_cells << std::endl;
      std::cout<< "Number of  cells= "<<cells_num_cells << std::endl;

      double cutoff=12.0; //after 12 macrocell of distance, the bare macrocell model gives the same result

      // loop over local cells
      for(int lc=0;lc<cells::num_local_cells;lc++){

         // reference global cell ID
         int i = cells_local_cell_array[lc];

         if(cells::num_atoms_in_cell[i]>0){

         	// Loop over all other cells to calculate contribution to local cell
          	for(int j=0;j<cells_num_cells;j++){

             	// Calculation of INTER TERM of Dipolar interaction
             	if(i!=j && cells_num_atoms_in_cell[j]>0){

                 	double tmp_rij_inter_xx = 0.0;
                 	double tmp_rij_inter_xy = 0.0;
                 	double tmp_rij_inter_xz = 0.0;

                 	double tmp_rij_inter_yy = 0.0;
                 	double tmp_rij_inter_yz = 0.0;
                 	double tmp_rij_inter_zz = 0.0;

                 	double rx = cells_coord_array_x[j] - cells_coord_array_x[i]; // Angstroms
                 	double ry = cells_coord_array_y[j] - cells_coord_array_y[i];
                 	double rz = cells_coord_array_z[j] - cells_coord_array_z[i];

                 	double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance
                 	double rij_1 = 1.0/rij;

                  // If distance between macro-cells > cutoff nm => continuum approach
                  if( (rij_1)/cells_macro_cell_size > cutoff){

                  	double ex = rx*rij;
                  	double ey = ry*rij;
                  	double ez = rz*rij;

                  	double rij3 = (rij*rij*rij); // Angstroms

                  	rij_inter_xx[lc][j] = ((3.0*ex*ex - 1.0)*rij3);
                  	rij_inter_xy[lc][j] = ( 3.0*ex*ey      )*rij3 ;
                  	rij_inter_xz[lc][j] = ( 3.0*ex*ez      )*rij3 ;

                  	rij_inter_yy[lc][j] = ((3.0*ey*ey - 1.0)*rij3);
                  	rij_inter_yz[lc][j] = ( 3.0*ey*ez      )*rij3 ;
                  	rij_inter_zz[lc][j] = ((3.0*ez*ez - 1.0)*rij3);
                  }
                  else if( (1.0/rij)/cells_macro_cell_size <= cutoff){
                     for(int pi=0; pi<cells_num_atoms_in_cell[i]; pi++){
                        for(int qj=0; qj<cells_num_atoms_in_cell[j]; qj++){

                           rx = cells_atom_in_cell_coords_array_x[j][qj] - cells_atom_in_cell_coords_array_x[i][pi];
                           ry = cells_atom_in_cell_coords_array_y[j][qj] - cells_atom_in_cell_coords_array_y[i][pi];
                           rz = cells_atom_in_cell_coords_array_z[j][qj] - cells_atom_in_cell_coords_array_z[i][pi];

                           rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz);  //Reciprocal of the distance
                           rij_1 = 1.0/rij;

                           if( rij_1==0.0 ) {
                              std::cout << ">>>>> (Inter) Warning: atoms are overlapping in cells i,j " << i << "\t" << j <<"<<<<<" << std::endl;
                              std::cout << "Cell and atomic coordinates used for calculation of Dipolar matrix" << std::endl;
                              std::cout << " xj= "  << cells_cell_coords_array_x[j]  << " yj= "  << cells_cell_coords_array_y[j]  << " zj= "  << cells_cell_coords_array_z[j] << std::endl;
                              std::cout << " xi= "  << cells_cell_coords_array_x[i]  << " yi= "  << cells_cell_coords_array_y[i]  << " zi= "  << cells_cell_coords_array_z[i] << std::endl;
//                            std::cout << " rx_cell= "  << rx_cell << " ry_cell= " << ry_cell << " rz_cell= " << rz_cell << std::endl;
                              std::cout << " xpi= " << cells_atom_in_cell_coords_array_x[i][pi] << " ypi= " << cells_atom_in_cell_coords_array_y[i][pi] << " zpi= " << cells_atom_in_cell_coords_array_z[i][pi] << std::endl;
                              std::cout << " xqj= " << cells_atom_in_cell_coords_array_x[j][qj] << " yqj= " << cells_atom_in_cell_coords_array_y[j][qj] << " zqj= " << cells_atom_in_cell_coords_array_z[j][qj] << std::endl;
                              std::cout << " rx= "  << rx << " ry= " << ry << " rz= " << rz << std::endl;
                              std::cout << " rij= " << 1.0/rij << " = " << (1.0/rij)/((cells_macro_cell_size*sqrt(3)+0.01)) << " macro-cells" << std::endl;
                              return;
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

                     rij_inter_xx[lc][j] =  (tmp_rij_inter_xx);
                     rij_inter_xy[lc][j] =  (tmp_rij_inter_xy);
                     rij_inter_xz[lc][j] =  (tmp_rij_inter_xz);

                     rij_inter_yy[lc][j] =  (tmp_rij_inter_yy);
                     rij_inter_yz[lc][j] =  (tmp_rij_inter_yz);
                     rij_inter_zz[lc][j] =  (tmp_rij_inter_zz);

                     rij_inter_xx[lc][j] = rij_inter_xx[lc][j]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                     rij_inter_xy[lc][j] = rij_inter_xy[lc][j]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                     rij_inter_xz[lc][j] = rij_inter_xz[lc][j]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));

                     rij_inter_yy[lc][j] = rij_inter_yy[lc][j]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                     rij_inter_yz[lc][j] = rij_inter_yz[lc][j]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                     rij_inter_zz[lc][j] = rij_inter_zz[lc][j]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                  }  // End of Inter part calculated atomicstically
               } // End of Inter part

               // Calculation of INTRA TERM of Dipolar interaction
               else if( i==j && cells::num_atoms_in_cell[j]>0){

               double tmp_rij_intra_xx = 0.0;
               double tmp_rij_intra_xy = 0.0;
               double tmp_rij_intra_xz = 0.0;

               double tmp_rij_intra_yy = 0.0;
               double tmp_rij_intra_yz = 0.0;
               double tmp_rij_intra_zz = 0.0;

               	for(int pi=0; pi<cells_num_atoms_in_cell[i]; pi++){
                 		for(int qj=0; qj<cells_num_atoms_in_cell[i]; qj++){
                  		if( pi!=qj ){

                   			double rx = cells_atom_in_cell_coords_array_x[i][qj] - cells_atom_in_cell_coords_array_x[i][pi];
                   			double ry = cells_atom_in_cell_coords_array_y[i][qj] - cells_atom_in_cell_coords_array_y[i][pi];
                   			double rz = cells_atom_in_cell_coords_array_z[i][qj] - cells_atom_in_cell_coords_array_z[i][pi];

                   			const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance

                   			if( 1.0/rij==0.0 ){
                     			std::cout << ">>>>> (Intra)  Warning: atoms are overlapping in cells i=\t" << i << "\t and j=\t" << j <<"\t<<<<<" << std::endl;
                     			std::cout << "Cell and atomic coordinates used for calculation of Dipolar matrix" << std::endl;
                     			std::cout << " xqj= " << cells_atom_in_cell_coords_array_x[j][qj] << " yqj= " << cells_atom_in_cell_coords_array_y[j][qj] << " zqj= " << cells_atom_in_cell_coords_arrayz[j][qj] << std::endl;
                     			std::cout << " xpi= " << cells_atom_in_cell_coords_array_x[i][pi] << " ypi= " << cells_atom_in_cell_coords_array_y[i][pi] << " zpi= " << cells_atom_in_cell_coords_arrayz[i][pi] << std::endl;
                     			std::cout << " xj= "  << cells_coord_atoms_array_x[j]  << " yj= "  << cells_coord_atoms_array_y[j]  << " zj= "  << cells_coord_atoms_array_y[j] << std::endl;
                     			std::cout << " xi= "  << cells_coord_atoms_array_x[i]  << " yi= "  << cells_coord_atoms_array_y[i]  << " zi= "  << cells_coord_atoms_array_y[i] << std::endl;
                     			std::cout << " rx= "  << rx << " ry= " << ry << " rz= " << rz << std::endl;
                     			std::cout << " rij= " << 1.0/rij << " = " << (1.0/rij)/((cells_macro_cell_size*sqrt(3)+0.01)) << " macro-cells" << std::endl;
                     			return;
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

                	rij_intra_xx[lc][i] =  (tmp_rij_intra_xx);
                	rij_intra_xy[lc][i] =  (tmp_rij_intra_xy);
                	rij_intra_xz[lc][i] =  (tmp_rij_intra_xz);

                	rij_intra_yy[lc][i] =  (tmp_rij_intra_yy);
                	rij_intra_yz[lc][i] =  (tmp_rij_intra_yz);
                	rij_intra_zz[lc][i] =  (tmp_rij_intra_zz);

                	rij_intra_xx[lc][j] = rij_intra_xx[lc][j]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                	rij_intra_xy[lc][j] = rij_intra_xy[lc][j]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                	rij_intra_xz[lc][j] = rij_intra_xz[lc][j]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));

                	rij_intra_yy[lc][j] = rij_intra_yy[lc][j]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                	rij_intra_yz[lc][j] = rij_intra_yz[lc][j]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
                	rij_intra_zz[lc][j] = rij_intra_zz[lc][j]/(double(cells_num_atoms_in_cell[i]) * double(cells_num_atoms_in_cell[j]));
               }
            }
			}
		}


      // timing function
      #ifdef MPICF
			double t1 = MPI_Wtime();
      #else
         time_t t1;
         t1 = time (NULL);
      #endif

      // now calculate fields
      dipole::internal::update_field();

      // timing function
      #ifdef MPICF
         double t2 = MPI_Wtime();
      #else
         time_t t2;
         t2 = time (NULL);
      #endif

      zlog << zTs() << "Time required for demag update: " << t2-t1 << "s." << std::endl;

    	return;

    }

} // end of dipole namespace

