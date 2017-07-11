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
      // Function to calculate Intra component of dipole tensors.
      //
      // The tensors between local cells with the cutoff range are calculated
      // explictly from the atomistic coordinates. Longer range tensors assume
      // the dipole-dipole form.
      //------------------------------------------------------------------------
      void compute_intra_tensor(const int i,
                                const int j,
                                const int lc,
                                std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                                std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                                std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                                std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z){

         // Here we start
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

         dipole::internal::rij_tensor_xx[lc][i] =  (tmp_rij_intra_xx);
         dipole::internal::rij_tensor_xy[lc][i] =  (tmp_rij_intra_xy);
         dipole::internal::rij_tensor_xz[lc][i] =  (tmp_rij_intra_xz);

         dipole::internal::rij_tensor_yy[lc][i] =  (tmp_rij_intra_yy);
      	dipole::internal::rij_tensor_yz[lc][i] =  (tmp_rij_intra_yz);
         dipole::internal::rij_tensor_zz[lc][i] =  (tmp_rij_intra_zz);

      	dipole::internal::rij_tensor_xx[lc][i] = dipole::internal::rij_tensor_xx[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
         dipole::internal::rij_tensor_xy[lc][i] = dipole::internal::rij_tensor_xy[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
         dipole::internal::rij_tensor_xz[lc][i] = dipole::internal::rij_tensor_xz[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));

         dipole::internal::rij_tensor_yy[lc][i] = dipole::internal::rij_tensor_yy[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
         dipole::internal::rij_tensor_yz[lc][i] = dipole::internal::rij_tensor_yz[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
         dipole::internal::rij_tensor_zz[lc][i] = dipole::internal::rij_tensor_zz[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));

      }  // End of funtion calculating Intra component of dipole tensor

   } // End of namespace internal
} // End of namespace dipole
