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
         //std::cout << std::endl << "dipole::internal::cells_local_cell_array[lc] = " << i << std::endl;
         fprintf(stderr,"lc = %d, i = %d, x = %f y = %f z = %f M = %e on rank = %d\n",lc,i,dipole::internal::cells_pos_and_mom_array[4*i+0],dipole::internal::cells_pos_and_mom_array[4*i+1],dipole::internal::cells_pos_and_mom_array[4*i+2],dipole::internal::cells_pos_and_mom_array[4*i+3],vmpi::my_rank);

        		if(dipole::internal::cells_num_atoms_in_cell[i]>0){

         	const double eightPI_three_cell_volume = 8.0*M_PI/(3.0*dipole::internal::cells_volume_array[i]);
//         	const double self_demag = demag::prefactor*eightPI_three_cell_volume;
         	double self_demag = eightPI_three_cell_volume;

         	// Add self-demagnetisation as mu_0/4_PI * 8PI/3V
         	dipole::cells_field_array_x[i]=self_demag*(cells::mag_array_x[i]/9.27400915e-24);
         	dipole::cells_field_array_y[i]=self_demag*(cells::mag_array_y[i]/9.27400915e-24);
         	dipole::cells_field_array_z[i]=self_demag*(cells::mag_array_z[i]/9.27400915e-24);

         	// Loop over all other cells to calculate contribution to local cell
         	for(int j=0;j<dipole::internal::cells_num_cells;j++){
           		if(dipole::internal::cells_num_atoms_in_cell[j]>0){

            		const double mx = cells::mag_array_x[j]/9.27400915e-24;
            		const double my = cells::mag_array_y[j]/9.27400915e-24;
            		const double mz = cells::mag_array_z[j]/9.27400915e-24;
                  fprintf(stderr,"lc %d i %d j %d mx(%d) %f my(%d) %f mz(%d) %f on rank = %d\n",lc,i,j,j,mx,j,my,j,mz,vmpi::my_rank);

                  //#ifdef MPICF
                  //if((dipole::internal::cells_pos_and_mom_array[4*i+0]!=dipole::internal::cells_pos_and_mom_array[4*j+0]) &&
                  //   (dipole::internal::cells_pos_and_mom_array[4*i+1]!=dipole::internal::cells_pos_and_mom_array[4*j+1]) &&
                  //   (dipole::internal::cells_pos_and_mom_array[4*i+2]!=dipole::internal::cells_pos_and_mom_array[4*j+2]) ){
                  //#else
            		if( i!=j ){
                  //#endif
                     fprintf(stderr," >> inter >> j = %d, x = %f y = %f z = %f M = %e on rank = %d\n",j,dipole::internal::cells_pos_and_mom_array[4*j+0],dipole::internal::cells_pos_and_mom_array[4*j+1],dipole::internal::cells_pos_and_mom_array[4*j+2],dipole::internal::cells_pos_and_mom_array[4*j+3],vmpi::my_rank);
                		dipole::cells_field_array_x[i]+=(mx*internal::rij_inter_xx[lc][j] + my*internal::rij_inter_xy[lc][j] + mz*internal::rij_inter_xz[lc][j]);
                		dipole::cells_field_array_y[i]+=(mx*internal::rij_inter_xy[lc][j] + my*internal::rij_inter_yy[lc][j] + mz*internal::rij_inter_yz[lc][j]);
                		dipole::cells_field_array_z[i]+=(mx*internal::rij_inter_xz[lc][j] + my*internal::rij_inter_yz[lc][j] + mz*internal::rij_inter_zz[lc][j]);
            		}
            		else{
                     fprintf(stderr," >> intra >> j = %d, x = %f y = %f z = %f M = %e on rank = %d\n",j,dipole::internal::cells_pos_and_mom_array[4*j+0],dipole::internal::cells_pos_and_mom_array[4*j+1],dipole::internal::cells_pos_and_mom_array[4*j+2],dipole::internal::cells_pos_and_mom_array[4*j+3],vmpi::my_rank);
                		dipole::cells_field_array_x[i]+=(mx*internal::rij_intra_xx[lc][i] + my*internal::rij_intra_xy[lc][i] + mz*internal::rij_intra_xz[lc][i]);
                		dipole::cells_field_array_y[i]+=(mx*internal::rij_intra_xy[lc][i] + my*internal::rij_intra_yy[lc][i] + mz*internal::rij_intra_yz[lc][i]);
                		dipole::cells_field_array_z[i]+=(mx*internal::rij_intra_xz[lc][i] + my*internal::rij_intra_yz[lc][i] + mz*internal::rij_intra_zz[lc][i]);
            		}
          		}
         	}
         	dipole::cells_field_array_x[i] = dipole::cells_field_array_x[i] * 9.27400915e-01;
         	dipole::cells_field_array_y[i] = dipole::cells_field_array_y[i] * 9.27400915e-01;
         	dipole::cells_field_array_z[i] = dipole::cells_field_array_z[i] * 9.27400915e-01;
     		}
    	}

    	for(int lc=0;lc<dipole::internal::cells_num_cells;lc++){
         //if(dipole::internal::cells_num_atoms_in_cell[lc]>0){
            //std::cout << lc << "\t";
            //std::cout << cells::mag_array_x[lc]<< "\t" << dipole::cells_field_array_x[lc] << "\t";
            //std::cout << cells::mag_array_y[lc]<< "\t" << dipole::cells_field_array_y[lc] << "\t";
            //std::cout << cells::mag_array_z[lc]<< "\t" << dipole::cells_field_array_z[lc] << std::endl;
            fprintf(stderr,"lc = %d\tmag_x = %e\tmag_y = %e\tmag_z = %e\tfield_x = %f\tfield_y = %f\tfield_z = %f on rank = %d\n",lc,cells::mag_array_x[lc],cells::mag_array_y[lc],cells::mag_array_z[lc],dipole::cells_field_array_x[lc],dipole::cells_field_array_y[lc],dipole::cells_field_array_z[lc],vmpi::my_rank);
         //}
      }

      MPI::COMM_WORLD.Barrier();

      fprintf(stderr,"\n\n  >>>>>>> Before reduction of cells field <<<<<< \n\n");

      //#ifdef MPICF
      //// Reduce field on all CPUs
      //MPI_Allreduce(MPI_IN_PLACE, &dipole::cells_field_array_x[0], dipole::cells_field_array_x.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      //MPI_Allreduce(MPI_IN_PLACE, &dipole::cells_field_array_y[0], dipole::cells_field_array_y.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      //MPI_Allreduce(MPI_IN_PLACE, &dipole::cells_field_array_z[0], dipole::cells_field_array_z.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      //////MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &dipole::cells_field_array_x[0], dipole::cells_field_array_x.size(), MPI_DOUBLE, MPI_SUM);
      //////MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &dipole::cells_field_array_y[0], dipole::cells_field_array_y.size(), MPI_DOUBLE, MPI_SUM);
      //////MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &dipole::cells_field_array_z[0], dipole::cells_field_array_z.size(), MPI_DOUBLE, MPI_SUM);
      //#endif

      fprintf(stderr,"\n\n  >>>>>>> After reduction of cells field <<<<<< \n\n");

      MPI::COMM_WORLD.Barrier();

    	for(int lc=0;lc<dipole::internal::cells_num_cells;lc++){
         fprintf(stderr,"lc = %d\tmag_x = %e\tmag_y = %e\tmag_z = %e\tfield_x = %f\tfield_y = %f\tfield_z = %f on rank = %d\n",lc,cells::mag_array_x[lc],cells::mag_array_y[lc],cells::mag_array_z[lc],dipole::cells_field_array_x[lc],dipole::cells_field_array_y[lc],dipole::cells_field_array_z[lc],vmpi::my_rank);
      }

      MPI::COMM_WORLD.Barrier();

	}

} // end of dipole namespace
