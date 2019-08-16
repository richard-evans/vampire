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

  void ha::calc_inter(int cell_i, int cell_j, int interaction_no, std::vector < std::vector < double > >& cells_atom_in_cell_coords_array_x, std::vector < std::vector < double > >& cells_atom_in_cell_coords_array_y, std::vector < std::vector < double > >& cells_atom_in_cell_coords_array_z){


     // instantiate timer
//     vutil::vtimer_t timer;
   //std::cout << dipole::cutoff*cells::macro_cell_size << std::endl;
   //timer.start();
    // create temporary variable to store components of tensor
    double tmp_rij_inter_xx = 0.0;
    double tmp_rij_inter_xy = 0.0;
    double tmp_rij_inter_xz = 0.0;

    double tmp_rij_inter_yy = 0.0;
    double tmp_rij_inter_yz = 0.0;
    double tmp_rij_inter_zz = 0.0;




       // Calculate distance vectors between cells
       double rx2 = ha::cell_positions_mom[4*cell_j+0] - ha::cell_positions_mom[4*cell_i+0];
       double ry2 = ha::cell_positions_mom[4*cell_j+1] - ha::cell_positions_mom[4*cell_i+1];
       double rz2 = ha::cell_positions_mom[4*cell_j+2] - ha::cell_positions_mom[4*cell_i+2];
       double rij = sqrt(rx2*rx2+ry2*ry2+rz2*rz2); //Reciprocal of the distance
       double rij_1 = 1.0/rij;
      //if (cell_i == 0) std::cout <<rij << std::endl;
      if (rij > dipole::cutoff*cells::macro_cell_size){
      //if (cell_i ==0)   std::cout << cell_j << '\t'<< rij << '\t' << "macrocell" <<std::endl;

                  // define unitarian distance vectors

         const double ex = rx2*rij_1;
         const double ey = ry2*rij_1;
         const double ez = rz2*rij_1;

         const double rij3 = (rij_1*rij_1*rij_1); // Angstroms

         // calculate dipolar matrix for 6 entries because of symmetry
         ha::rij_tensor_xx[interaction_no] += ((3.0*ex*ex - 1.0)*rij3);
         ha::rij_tensor_xy[interaction_no] += ((3.0*ex*ey      )*rij3);
         ha::rij_tensor_xz[interaction_no] += ((3.0*ex*ez      )*rij3);

         ha::rij_tensor_yy[interaction_no] += ((3.0*ey*ey - 1.0)*rij3);
         ha::rij_tensor_yz[interaction_no] += ((3.0*ey*ez      )*rij3);
         ha::rij_tensor_zz[interaction_no] += ((3.0*ez*ez - 1.0)*rij3);
         // if (ha::rij_tensor_zz[interaction_no] != ha::rij_tensor_zz[interaction_no])  {std::cout << "zero" <<std::endl;
         // std::cin.get();}
         //std::cout << cell_i << '\t' << cell_j << '\t' << ha::cell_positions_mom[4*cell_j+0] << '\t'  << ha::cell_positions_mom[4*cell_i+0] << "\t" << ha::rij_tensor_xx[interaction_no] << std::endl;

      }


         //--------------------------------------------------------------------------
         // If distance between macro-cells < cutoff ==> apply inter-intra method
         //--------------------------------------------------------------------------
          else if(rij <= dipole::cutoff*cells::macro_cell_size){
         //    if (cell_i ==0)   std::cout << cell_j << '\t'<< rij << '\t' << "tensor" << "\t" << ha::num_atoms_in_cell[cell_i] << '\t' << ha::num_atoms_in_cell[cell_j] << std::endl;


           if  (cell_j > ha::num_zero_level_cells) std::cout << ha::num_zero_level_cells << '\t' << cell_j << '\t' << rij <<  "\t" << dipole::cutoff*cells::macro_cell_size << std::endl;
          //std::cout << "ATOM" << '\t' << ha::num_atoms_in_cell[cell_i] << std::endl;
            //if (cell_j > 120) std::cout << cell_i << "\t" << cell_j << std::endl;
             for(int pi=0; pi<ha::num_atoms_in_cell[cell_i]; pi++){

               const double cix = cells_atom_in_cell_coords_array_x[cell_i][pi];
               const double ciy = cells_atom_in_cell_coords_array_y[cell_i][pi];
               const double ciz = cells_atom_in_cell_coords_array_z[cell_i][pi];

                 for(int qj=0; qj<ha::num_atoms_in_cell[cell_j]; qj++){

                  const double dx = cells_atom_in_cell_coords_array_x[cell_j][qj] - cix;
                  const double dy = cells_atom_in_cell_coords_array_y[cell_j][qj] - ciy;
                  const double dz = cells_atom_in_cell_coords_array_z[cell_j][qj] - ciz;

                  rij = 1.0/sqrt(dx*dx+dy*dy+dz*dz);  //Reciprocal of the distance

                  const double ex = dx*rij;
                  const double ey = dy*rij;
                  const double ez = dz*rij;

                  const double rij3 = (rij*rij*rij); // Angstroms

                  tmp_rij_inter_xx += ((3.0*ex*ex - 1.0)*rij3);
                  tmp_rij_inter_xy += ((3.0*ex*ey      )*rij3);
                  tmp_rij_inter_xz += ((3.0*ex*ez      )*rij3);

                  tmp_rij_inter_yy += ((3.0*ey*ey - 1.0)*rij3);
                  tmp_rij_inter_yz += ((3.0*ey*ez      )*rij3);
                  tmp_rij_inter_zz += ((3.0*ez*ez - 1.0)*rij3);
//                   //std::cout<<pi << '\t' << qj << '\t' << ex << '\t' << ey << '\t' << ez << "\t" << tmp_rij_inter_xx <<std::endl;
//std::cout << cell_j << '\t' << rij3 << '\t'  << ha::num_atoms_in_cell[cell_i] << '\t' << ha::num_atoms_in_cell[cell_j] << '\t' << dx << '\t' << dy << '\t' << dz << "\t" << tmp_rij_inter_xx << '\t' << tmp_rij_inter_xy << '\t' << tmp_rij_inter_xz << '\t' <<std::endl;
                 }
              }
            //std::cout << interaction_no << '\t' << ha::rij_tensor_xx.size() <<std::endl;
            ha::rij_tensor_xx[interaction_no] =  (tmp_rij_inter_xx);
            ha::rij_tensor_xy[interaction_no] =  (tmp_rij_inter_xy);
            ha::rij_tensor_xz[interaction_no] =  (tmp_rij_inter_xz);

            ha::rij_tensor_yy[interaction_no] =  (tmp_rij_inter_yy);
            ha::rij_tensor_yz[interaction_no] =  (tmp_rij_inter_yz);
            ha::rij_tensor_zz[interaction_no] =  (tmp_rij_inter_zz);
//
// //   std::cout << "atom" <<  '\t' << cell_i <<'\t' << cell_j << "\t" << ha::rij_tensor_xx[interaction_no] << "\t" << ha::rij_tensor_xy[interaction_no] << '\t' << ha::rij_tensor_xz[interaction_no] << std::endl;
         //
         // //    // Normalisation by the number of atoms in the cell. This is required for the correct evaluation of the field in the update.cpp routine
            ha::rij_tensor_xx[interaction_no] = ha::rij_tensor_xx[interaction_no]/(double(ha::num_atoms_in_cell[cell_i]) * double(ha::num_atoms_in_cell[cell_j]));
            ha::rij_tensor_xy[interaction_no] = ha::rij_tensor_xy[interaction_no]/(double(ha::num_atoms_in_cell[cell_i]) * double(ha::num_atoms_in_cell[cell_j]));
            ha::rij_tensor_xz[interaction_no] = ha::rij_tensor_xz[interaction_no]/(double(ha::num_atoms_in_cell[cell_i]) * double(ha::num_atoms_in_cell[cell_j]));

            ha::rij_tensor_yy[interaction_no] = ha::rij_tensor_yy[interaction_no]/(double(ha::num_atoms_in_cell[cell_i]) * double(ha::num_atoms_in_cell[cell_j]));
            ha::rij_tensor_yz[interaction_no] = ha::rij_tensor_yz[interaction_no]/(double(ha::num_atoms_in_cell[cell_i]) * double(ha::num_atoms_in_cell[cell_j]));
            ha::rij_tensor_zz[interaction_no] = ha::rij_tensor_zz[interaction_no]/(double(ha::num_atoms_in_cell[cell_i]) * double(ha::num_atoms_in_cell[cell_j]));

            }  // End of Inter part calculated atomicstically
//  timer.stop();
// // std::cout <<"cell\t" << cell_i <<   "\tdone! [ " << timer.elapsed_time() << " s ]" << std::endl;
// zlog << zTs() << "cell\t" << cell_i <<   "\tPrecalculation of rij matrix for dipole calculation complete. Time taken: " << timer.elapsed_time() << " s"<< std::endl;

      return;

  } //function

} // end of hierarchical namespace
