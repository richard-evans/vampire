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
#include "sim.hpp"

// dipole module headers
#include "internal.hpp"


namespace dp = dipole::internal;


namespace dipole{

   //-----------------------------------------------------------------------------
   // Function for updating local temperature fields
   //-----------------------------------------------------------------------------

#ifdef FFT

	void dipole::internal::update_field_fft(){

      if(err::check==true){
   		terminaltextcolor(RED);
   		std::cerr << "demag::fft_update has been called " << vmpi::my_rank << std::endl;
   		terminaltextcolor(WHITE);
   	}


		for (int i=0 ; i<dp::num_macro_cells_x ; i++){
         for (int j=0 ; j<dp::num_macro_cells_y ; j++){
             for (int k=0 ; k<dp::num_macro_cells_z ; k++){
						int id = (i*dp::num_macro_cells_x+j)*dp::num_macro_cells_y+k;
						Mx_in[id][0]=0;
						Mx_in[id][1]=0;
						My_in[id][0]=0;
						My_in[id][1]=0;
						Mz_in[id][0]=0;
						Mz_in[id][1]=0;


						Mx_out[id][0]=0;
						Mx_out[id][1]=0;
						My_out[id][0]=0;
						My_out[id][1]=0;
						Mz_out[id][0]=0;
						Mz_out[id][1]=0;


						Hx_in[id][0]=0;
						Hx_in[id][1]=0;
						Hy_in[id][0]=0;
						Hy_in[id][1]=0;
						Hz_in[id][0]=0;
						Hz_in[id][1]=0;


						Hx_out[id][0]=0;
						Hx_out[id][1]=0;
						Hy_out[id][0]=0;
						Hy_out[id][1]=0;
						Hz_out[id][0]=0;
						Hz_out[id][1]=0;

				 }
			 }
		}


   	//		std::cout << "a" <<std::endl;
   		int cell = 0;
   		for (int i=0 ; i<dp::num_macro_cells_x; i++){
   			for (int j=0 ; j<dp::num_macro_cells_y; j++){
   				for (int k=0 ; k<dp::num_macro_cells_z; k++){
						int id = (i*dp::num_macro_cells_x+j)*dp::num_macro_cells_y+k;
   					Mx_in[id][0] = cells::mag_array_x[cell]/9.27400915e-24;
   					My_in[id][0] = cells::mag_array_y[cell]/9.27400915e-24;
   					Mz_in[id][0] = cells::mag_array_z[cell]/9.27400915e-24;
                  //if (cell == 0) std::cout << "A" << Mx(i,j,k)[0]  << "\t" <<   My(i,j,k)[0]  << "\t" <<   Mz(i,j,k)[0]  << "\t" <<  std::endl;

					//	std::cout << cells::mag_array_x[cell] << '\t' << cells::mag_array_y[cell] << '\t' << cells::mag_array_z[cell] << '\t' << Mx(i,j,k)[0] << '\t' << My(i,j,k)[0] << '\t' << Mz(i,j,k)[0] <<std::endl;
					   cell ++;

   				}
   			 }
   		}


   	   fftw_plan MxP,MyP,MzP;

   		//std::cout << 'g' <<std::endl;
         MxP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::Mx_in,dp::Mx_out,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(MxP);
         MyP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::My_in,dp::My_out,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(MyP);
         MzP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::Mz_in,dp::Mz_out,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(MzP);
   		//std::cout << 'h' <<std::endl;

cell = 0;
         // performs the converlusion between Nk and Mk
      for (int i=0 ; i<2*dp::num_macro_cells_x ; i++){
         for (int j=0 ; j<2*dp::num_macro_cells_y ; j++){
             for (int k=0 ; k<2*dp::num_macro_cells_z ; k++){

					 int id = (i*dp::num_macro_cells_x+j)*dp::num_macro_cells_y+k;
					 Hx_in[id][0] = dp::N2xx[id][0]*dp::Mx_out[id][0] + dp::N2xy[id][0]*dp::My_out[id][0] + dp::N2xz[id][0]*dp::Mz_out[id][0]; //summing the real part
					 Hx_in[id][0] -= (dp::N2xx[id][1]*dp::Mx_out[id][1] + dp::N2xy[id][1]*dp::My_out[id][1] + dp::N2xz[id][1]*dp::Mz_out[id][1]);

					 Hx_in[id][1] = dp::N2xx[id][0]*dp::Mx_out[id][1] + dp::N2xy[id][0]*dp::My_out[id][1] + dp::N2xz[id][0]*dp::Mz_out[id][1];
					 Hx_in[id][1] += (dp::N2xx[id][1]*dp::Mx_out[id][0] + dp::N2xy[id][1]*dp::My_out[id][0] + dp::N2xz[id][1]*dp::Mz_out[id][0]);

					 Hy_in[id][0] = dp::N2yx[id][0]*dp::Mx_out[id][0] + dp::N2yy[id][0]*dp::My_out[id][0] + dp::N2yz[id][0]*dp::Mz_out[id][0];
					 Hy_in[id][0] -= (dp::N2yx[id][1]*dp::Mx_out[id][1] + dp::N2yy[id][1]*dp::My_out[id][1] + dp::N2yz[id][1]*dp::Mz_out[id][1]);

					 Hy_in[id][1] = dp::N2yx[id][0]*dp::Mx_out[id][1] + dp::N2yy[id][0]*dp::My_out[id][1] + dp::N2yz[id][0]*dp::Mz_out[id][1];
					 Hy_in[id][1] += (dp::N2yx[id][1]*dp::Mx_out[id][0] + dp::N2yy[id][1]*dp::My_out[id][0] + dp::N2yz[id][1]*dp::Mz_out[id][0]);

					 Hz_in[id][0] = dp::N2zx[id][0]*dp::Mx_out[id][0] + dp::N2zy[id][0]*dp::My_out[id][0] + dp::N2zz[id][0]*dp::Mz_out[id][0]; //summing the real part
					 Hz_in[id][0] -= (dp::N2zx[id][1]*dp::Mx_out[id][1] + dp::N2zy[id][1]*dp::My_out[id][1] + dp::N2zz[id][1]*dp::Mz_out[id][1]);

					 Hz_in[id][1] = dp::N2zx[id][0]*dp::Mx_out[id][1] + dp::N2zy[id][0]*dp::My_out[id][1] + dp::N2zz[id][0]*dp::Mz_out[id][1];
					 Hz_in[id][1] += (dp::N2zx[id][1]*dp::Mx_out[id][0] + dp::N2zy[id][1]*dp::My_out[id][0] + dp::N2zz[id][1]*dp::Mz_out[id][0]);
					 cell++;

             }
         }
      }

//std::cin.get();

      // performs the backward transform to give the dipole field, Hx, Hy, Hz
      fftw_plan HxP,HyP,HzP;

      HxP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::Hx_in,dp::Hx_out,FFTW_BACKWARD,FFTW_ESTIMATE);
      fftw_execute(HxP);
      HyP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::Hy_in,dp::Hy_out,FFTW_BACKWARD,FFTW_ESTIMATE);
      fftw_execute(HyP);
      HzP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::Hz_in,dp::Hz_out,FFTW_BACKWARD,FFTW_ESTIMATE);
      fftw_execute(HzP);


      for(int lc=0;lc<dipole::internal::cells_num_local_cells;lc++){

         //int i = dipole::internal::cells_local_cell_array[lc];
         int i = cells::cell_id_array[lc];
         //std::cout << std::endl << "dipole::internal::cells_local_cell_array[lc] = " << i << std::endl;
         //fprintf(stderr,"lc = %d, i = %d, x = %f y = %f z = %f M = %e on rank = %d\n",lc,i,dipole::internal::cells_pos_and_mom_array[4*i+0],dipole::internal::cells_pos_and_mom_array[4*i+1],dipole::internal::cells_pos_and_mom_array[4*i+2],dipole::internal::cells_pos_and_mom_array[4*i+3],vmpi::my_rank);

            if(dipole::internal::cells_num_atoms_in_cell[i]>0){

            const double eightPI_three_cell_volume = 8.0*M_PI/(3.0*dipole::internal::cells_volume_array[i]);
      //         	const double self_demag = demag::prefactor*eightPI_three_cell_volume;
            double self_demag = eightPI_three_cell_volume;
            //fprintf(stderr,"  $$$$$$$$$ dipole::internal::cells_volume_array[%d] = %f self_demag = %e on rank = %d\n",i,dipole::internal::cells_volume_array[i],self_demag,vmpi::my_rank);

            // Add self-demagnetisation as mu_0/4_PI * 8PI/3V
            dipole::cells_field_array_x[i]=self_demag*(cells::mag_array_x[i]/9.27400915e-24);
            dipole::cells_field_array_y[i]=self_demag*(cells::mag_array_y[i]/9.27400915e-24);
            dipole::cells_field_array_z[i]=self_demag*(cells::mag_array_z[i]/9.27400915e-24);

         }
      }

   	cell = 0;
   	for (int i=0 ; i<dp::num_macro_cells_x ; i++){
   		for (int j=0 ; j<dp::num_macro_cells_y ; j++){
   			for (int k=0 ; k<dp::num_macro_cells_z ; k++){
   			//	if (cell ==0 && sim::time % 1000 == 0) std::cout << "cell field" << '\t' << sim::temperature << '\t'<< Hx(i,j,k)[0]/eight_num_cells << '\t' << Hy(i,j,k)[0]/eight_num_cells << '\t' << Hz(i,j,k)[0]/eight_num_cells << '\t' << std::endl;
				int id = (i*dp::num_macro_cells_x+j)*dp::num_macro_cells_y+k;
					dipole::cells_field_array_x[cell] += Hx_out[id][0]/dp::eight_num_cells;
   				dipole::cells_field_array_y[cell] += Hy_out[id][0]/dp::eight_num_cells;
   				dipole::cells_field_array_z[cell] += Hz_out[id][0]/dp::eight_num_cells;
               dipole::cells_field_array_x[cell] *= 9.27400915e-01;
               dipole::cells_field_array_y[cell] *= 9.27400915e-01;
               dipole::cells_field_array_z[cell] *= 9.27400915e-01;
            //   if (cell ==0) std::cout << "cell field" << '\t' << sim::temperature << '\t'<< dipole::cells_field_array_x[cell] << '\t' << dipole::cells_field_array_y[cell] << '\t' << dipole::cells_field_array_z[cell] << '\t' << std::endl;

   		//	if (cell == 0)	std::cout << "fft" <<   '\t' << cell << '\t' << cells::num_cells << '\t' << dipole::cells_field_array_x[cell] <<  '\t' << dipole::cells_field_array_y[cell] << '\t' << dipole::cells_field_array_z[cell] << std::endl;

   				cell++;
   			}
   		}
   	}

   }
#endif
}
