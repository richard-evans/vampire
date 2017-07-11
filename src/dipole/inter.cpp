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
      // Function to calculate inter component of dipole tensors.
      //
      // The tensors between local cells with the cutoff range are calculated
      // explictly from the atomistic coordinates. Longer range tensors assume
      // the dipole-dipole form.
      //------------------------------------------------------------------------
      void compute_inter_tensor(const double cells_macro_cell_size,
                                const int i,
                                const int j,
                                const int lc,
                                std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                                //std::vector<double>& cells_pos_and_mom_array, // array to store positions and moment of cells
                                std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                                std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                                std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z){

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
	         dipole::internal::rij_tensor_xx[lc][j] = ((3.0*ex*ex - 1.0)*rij3);
	         dipole::internal::rij_tensor_xy[lc][j] = ( 3.0*ex*ey      )*rij3 ;
	         dipole::internal::rij_tensor_xz[lc][j] = ( 3.0*ex*ez      )*rij3 ;

	         dipole::internal::rij_tensor_yy[lc][j] = ((3.0*ey*ey - 1.0)*rij3);
	         dipole::internal::rij_tensor_yz[lc][j] = ( 3.0*ey*ez      )*rij3 ;
	         dipole::internal::rij_tensor_zz[lc][j] = ((3.0*ez*ez - 1.0)*rij3);

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

            dipole::internal::rij_tensor_xx[lc][j] =  (tmp_rij_inter_xx);
            dipole::internal::rij_tensor_xy[lc][j] =  (tmp_rij_inter_xy);
            dipole::internal::rij_tensor_xz[lc][j] =  (tmp_rij_inter_xz);

            dipole::internal::rij_tensor_yy[lc][j] =  (tmp_rij_inter_yy);
            dipole::internal::rij_tensor_yz[lc][j] =  (tmp_rij_inter_yz);
            dipole::internal::rij_tensor_zz[lc][j] =  (tmp_rij_inter_zz);

            // Normalisation by the number of atoms in the cell. This is required for the correct evaluation of the field in the update.cpp routine
            dipole::internal::rij_tensor_xx[lc][j] = dipole::internal::rij_tensor_xx[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
            dipole::internal::rij_tensor_xy[lc][j] = dipole::internal::rij_tensor_xy[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
            dipole::internal::rij_tensor_xz[lc][j] = dipole::internal::rij_tensor_xz[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));

            dipole::internal::rij_tensor_yy[lc][j] = dipole::internal::rij_tensor_yy[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
            dipole::internal::rij_tensor_yz[lc][j] = dipole::internal::rij_tensor_yz[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
            dipole::internal::rij_tensor_zz[lc][j] = dipole::internal::rij_tensor_zz[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
         }  // End of Inter part calculated atomicstically
      }  // End of funtion calculating inter component of dipole tensor

   } // End of namespace internal
} // End of namespace dipole
