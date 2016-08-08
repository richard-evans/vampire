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
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>


// Vampire headers
#include "dipole.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "cells.hpp"
#include "errors.hpp"

// dipole module headers
#include "internal.hpp"

namespace dipole{

   //-----------------------------------------------------------------------------
   // Function for updating local temperature fields
   //-----------------------------------------------------------------------------


	void dipole::internal::update_field(){

      if(!dipole::activated) return;

		if(err::check==true){
			terminaltextcolor(RED);
			std::cerr << "dipole::update has been called " << vmpi::my_rank << std::endl;
			terminaltextcolor(WHITE);
		}

		// loop over local cells
    	for(int lc=0;lc<dipole::internal::cells_num_local_cells;lc++){

        	int i = dipole::internal::cells_local_cell_array[lc];

        		if(dipole::internal::cells_num_atoms_in_cell[i]>0){

         	const double eightPI_three_cell_volume = 8.0*M_PI/(3.0*dipole::internal::cells_volume_array[i]);
//         	const double self_demag = demag::prefactor*eightPI_three_cell_volume;
         	double self_demag = eightPI_three_cell_volume;

         	// Add self-demagnetisation as mu_0/4_PI * 8PI/3V
         	dipole::cells_field_array_x[i]=self_demag*(dipole::internal::cells_mag_array_x[i]/9.27400915e-24);
         	dipole::cells_field_array_y[i]=self_demag*(dipole::internal::cells_mag_array_y[i]/9.27400915e-24);
         	dipole::cells_field_array_z[i]=self_demag*(dipole::internal::cells_mag_array_z[i]/9.27400915e-24);

         	// Loop over all other cells to calculate contribution to local cell
         	for(int j=0;j<dipole::internal::cells_num_cells;j++){
           		if(dipole::internal::cells_num_atoms_in_cell[j]>0){

            		const double mx = dipole::internal::cells_mag_array_x[j]/9.27400915e-24;
            		const double my = dipole::internal::cells_mag_array_y[j]/9.27400915e-24;
            		const double mz = dipole::internal::cells_mag_array_z[j]/9.27400915e-24;

            		if( i!=j ){
                		dipole::cells_field_array_x[i]+=(mx*internal::rij_inter_xx[lc][j] + my*internal::rij_inter_xy[lc][j] + mz*internal::rij_inter_xz[lc][j]);
                		dipole::cells_field_array_y[i]+=(mx*internal::rij_inter_xy[lc][j] + my*internal::rij_inter_yy[lc][j] + mz*internal::rij_inter_yz[lc][j]);
                		dipole::cells_field_array_z[i]+=(mx*internal::rij_inter_xz[lc][j] + my*internal::rij_inter_yz[lc][j] + mz*internal::rij_inter_zz[lc][j]);
            		}
            		else{
                		dipole::cells_field_array_x[i]+=(mx*internal::rij_intra_xx[lc][i] + my*internal::rij_intra_xy[lc][i] + mz*internal::rij_intra_xz[lc][i]);
                		dipole::cells_field_array_y[i]+=(mx*internal::rij_intra_xy[lc][i] + my*internal::rij_intra_yy[lc][i] + mz*internal::rij_intra_yz[lc][i]);
                		dipole::cells_field_array_z[i]+=(mx*internal::rij_intra_xz[lc][i] + my*internal::rij_intra_yz[lc][i] + mz*internal::rij_intra_zz[lc][i]);
            		}
          		}
         	}
         	dipole::cells_field_array_x[i] *= 9.27400915e-01;
         	dipole::cells_field_array_y[i] *= 9.27400915e-01;
         	dipole::cells_field_array_z[i] *= 9.27400915e-01;
            /*std::cout << dipole::cells_field_array_x[i] << "\t";
            std::cout << dipole::cells_field_array_y[i] << "\t";
            std::cout << dipole::cells_field_array_z[i] << std::endl; */
     		}
    	}
	}

} // end of dipole namespace
