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
      void initialize_macrocell_solver(){

         // calculate matrix prefactors
         zlog << zTs() << "Precalculating rij matrix for dipole calculation using macrocell solver... " << std::endl;
         std::cout     << "Precalculating rij matrix for dipole calculation using macrocell solver"     << std::flush;

         // instantiate timer
         vutil::vtimer_t timer;

         // start timer
         timer.start();

         // loop over local cells
         for(int lc=0;lc<dipole::internal::cells_num_local_cells;lc++){

            // print out progress to screen
            if(lc % (dipole::internal::cells_num_local_cells/10) == 0) std::cout << "." << std::flush;

            // reference global cell ID
            //int i = dipole::internal::cells_local_cell_array[lc];
            int i = cells::cell_id_array[lc];

            // check that the cell constains at least one atom
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

                  /*==========================================================*/
                  /* Calculation of intra part of dipolar tensor              */
                  /*==========================================================*/
                  // ** Need to fix this !!!! ** //
                  else if( i==j && dipole::internal::cells_num_atoms_in_cell[j]>0){

                     const double third = 1.0/3.0;

                   	dipole::internal::rij_intra_xx[lc][i] = third;
                   	dipole::internal::rij_intra_xy[lc][i] = 0.0;
                   	dipole::internal::rij_intra_xz[lc][i] = 0.0;

                   	dipole::internal::rij_intra_yy[lc][i] = third;
                   	dipole::internal::rij_intra_yz[lc][i] = 0.0;
                   	dipole::internal::rij_intra_zz[lc][i] = third;

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
