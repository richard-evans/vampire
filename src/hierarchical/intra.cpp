//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2018. All rights reserved.
//
//   Email: sarah.jenkins@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
// C++ standard library headers
#include <string>
#include <cmath>
#include <cstdlib>
#include <iostream>
// C library headers
#include <fenv.h>
#include <signal.h>
#include <math.h>
// Vampire headers
#include "cells.hpp" // needed for cells::cell_id_array but to be removed
#include "dipole.hpp"
#include "vio.hpp"
#include "vutil.hpp"
#include "atoms.hpp"
// Vampire headers
#include "hierarchical.hpp"


// hierarchical module headers
#include "internal.hpp"

namespace ha = hierarchical::internal;

namespace hierarchical{

  void ha::calc_intra(int cell_i, int cell_j, int interaction_no, std::vector < std::vector < double > >& cells_atom_in_cell_coords_array_x, std::vector < std::vector < double > >& cells_atom_in_cell_coords_array_y, std::vector < std::vector < double > >& cells_atom_in_cell_coords_array_z){


     const int mmax = ha::num_atoms_in_cell[cell_i];
     double rxx =0;
     double rxy =0;
     double rxz =0;
     double ryy =0;
     double ryz =0;
     double rzz =0;

       for(int pi=0; pi<mmax; pi++){

       const double cix = cells_atom_in_cell_coords_array_x[cell_i][pi];
       const double ciy = cells_atom_in_cell_coords_array_y[cell_i][pi];
       const double ciz = cells_atom_in_cell_coords_array_z[cell_i][pi];
       // use double loops to avoid if pi != qj statement
         //   std::cout << cell_j << "\t" << pi << '\t' << cix  << '\t' << ciy << '\t' << ciz << std::endl;
       for(int qj=0; qj<pi; qj++){

         const double rx = cells_atom_in_cell_coords_array_x[cell_j][qj] - cix;
         const double ry = cells_atom_in_cell_coords_array_y[cell_j][qj] - ciy;
         const double rz = cells_atom_in_cell_coords_array_z[cell_j][qj] - ciz;


         const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance



         const double ex = rx*rij;
         const double ey = ry*rij;
         const double ez = rz*rij;


         const double rij3 = (rij*rij*rij); // Angstroms



         rxx += ((3.0*ex*ex - 1.0)*rij3);
         rxy += ((3.0*ex*ey      )*rij3);
         rxz += ((3.0*ex*ez      )*rij3);

         ryy += ((3.0*ey*ey - 1.0)*rij3);
         ryz += ((3.0*ey*ez      )*rij3);
         rzz += ((3.0*ez*ez - 1.0)*rij3);

       }
       for(int qj=pi+1; qj<mmax; qj++){

          const double rx = cells_atom_in_cell_coords_array_x[cell_j][qj] - cix;
          const double ry = cells_atom_in_cell_coords_array_y[cell_j][qj] - ciy;
          const double rz = cells_atom_in_cell_coords_array_z[cell_j][qj] - ciz;

          const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance

          const double ex = rx*rij;
          const double ey = ry*rij;
          const double ez = rz*rij;

          const double rij3 = (rij*rij*rij); // Angstroms

          rxx += ((3.0*ex*ex - 1.0)*rij3);
          rxy += ((3.0*ex*ey      )*rij3);
          rxz += ((3.0*ex*ez      )*rij3);

          ryy += ((3.0*ey*ey - 1.0)*rij3);
          ryz += ((3.0*ey*ez      )*rij3);
          rzz += ((3.0*ez*ez - 1.0)*rij3);

       }
   }
   ha::rij_tensor_xx[interaction_no] += rxx/(double(ha::num_atoms_in_cell[cell_i]) * double(ha::num_atoms_in_cell[cell_i]));
   ha::rij_tensor_xy[interaction_no] += rxy/(double(ha::num_atoms_in_cell[cell_i]) * double(ha::num_atoms_in_cell[cell_i]));
   ha::rij_tensor_xz[interaction_no] += rxz/(double(ha::num_atoms_in_cell[cell_i]) * double(ha::num_atoms_in_cell[cell_i]));

   ha::rij_tensor_yy[interaction_no] += ryy/(double(ha::num_atoms_in_cell[cell_i]) * double(ha::num_atoms_in_cell[cell_i]));
   ha::rij_tensor_yz[interaction_no] += ryz/(double(ha::num_atoms_in_cell[cell_i]) * double(ha::num_atoms_in_cell[cell_i]));
   ha::rij_tensor_zz[interaction_no] += rzz/(double(ha::num_atoms_in_cell[cell_i]) * double(ha::num_atoms_in_cell[cell_i]));

   if (ha::rij_tensor_xx[interaction_no]*ha::rij_tensor_xx[interaction_no] < 1e-20) ha::rij_tensor_xx[interaction_no] = 0;
   if (ha::rij_tensor_xy[interaction_no]*ha::rij_tensor_xy[interaction_no] < 1e-20) ha::rij_tensor_xy[interaction_no] = 0;
   if (ha::rij_tensor_xz[interaction_no]*ha::rij_tensor_xz[interaction_no] < 1e-20) ha::rij_tensor_xz[interaction_no] = 0;

   //std::cout << "atom" <<  '\t' << cell_i <<'\t' << cell_j << "\t" << ha::rij_tensor_xx[interaction_no] << "\t" << ha::rij_tensor_xy[interaction_no] << '\t' << ha::rij_tensor_xz[interaction_no] << std::endl;


  return;

  } //function

} // end of hierarchical namespace
