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



   void dipole::internal::initialize_fft_solver(){
      #ifdef FFT
      // allocate arrays to store data [nloccell x ncells]
      if(dipole::fft==true) {

         // calculate matrix prefactors
         zlog << zTs() << "Precalculating rij matrix for dipole calculation... " << std::endl;



      // determine number of cells in each direction (with small shift to prevent the fence post problem)
      dp::num_macro_cells_x = static_cast<unsigned int>(ceil((cs::system_dimensions[0]+0.01)/cells::macro_cell_size[0]));
      dp::num_macro_cells_y = static_cast<unsigned int>(ceil((cs::system_dimensions[1]+0.01)/cells::macro_cell_size[1]));
      dp::num_macro_cells_z = static_cast<unsigned int>(ceil((cs::system_dimensions[2]+0.01)/cells::macro_cell_size[2]));

      dp::eight_num_cells = 8*dp::num_macro_cells_x*dp::num_macro_cells_y*dp::num_macro_cells_z;


      dp::N2xx =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2xy =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2xz =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2yx =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2yy =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2yz =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2zx =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2zy =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2zz =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);

      dp::N2xx0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2xy0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2xz0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2yx0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2yy0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2yz0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2zx0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2zy0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::N2zz0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);

      dp::Mx_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::My_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::Mz_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);

      dp::Mx_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::My_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::Mz_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);

      dp::Hx_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::Hy_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::Hz_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);

      dp::Hx_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::Hy_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);
      dp::Hz_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dp::eight_num_cells);

      //   std::cout << "b" <<std::endl;

      for(unsigned int i = 0 ; i < dp::num_macro_cells_x ; i++)
      {
          for(unsigned int j = 0 ; j < dp::num_macro_cells_y; j++)
          {
            for(unsigned int k = 0 ; k < dp::num_macro_cells_z ; k++)
            {
                int id = (i*dp::num_macro_cells_x+j)*dp::num_macro_cells_y+k;
                  dp::N2xx0[id][0]=0;
                  dp::N2xx0[id][1]=0;
                  dp::N2xx[id][0] =0;
                  dp::N2xx[id][1] =0;
                  dp::N2xy0[id][0]=0;
                  dp::N2xy0[id][1]=0;
                  dp::N2xy[id][0] =0;
                  dp::N2xy[id][1] =0;
                  dp::N2xz0[id][0]=0;
                  dp::N2xz0[id][1]=0;

                  dp::N2yx0[id][0]=0;
                  dp::N2yx0[id][1]=0;
                  dp::N2yx[id][0] =0;
                  dp::N2yx[id][1] =0;
                  dp::N2yy0[id][0]=0;
                  dp::N2yy0[id][1]=0;
                  dp::N2yy[id][0] =0;
                  dp::N2yy[id][1] =0;
                  dp::N2yz0[id][0]=0;
                  dp::N2yz0[id][1]=0;

                  dp::N2zx0[id][0]=0;
                  dp::N2zx0[id][1]=0;
                  dp::N2zx[id][0] =0;
                  dp::N2zx[id][1] =0;
                  dp::N2zy0[id][0]=0;
                  dp::N2zy0[id][1]=0;
                  dp::N2zy[id][0] =0;
                  dp::N2zy[id][1] =0;
                  dp::N2zz0[id][0]=0;
                  dp::N2zz0[id][1]=0;
            }
          }
      }

      double ii,jj,kk;
         for(int i=0;i<dp::num_macro_cells_x*2;i++){
            if (i >= dp::num_macro_cells_x) ii = i - 2*dp::num_macro_cells_x;
            else ii = i;
            for(int j=0;j<dp::num_macro_cells_y*2;j++){
               if (j >= dp::num_macro_cells_y) jj = j - 2*dp::num_macro_cells_y;
               else jj = j;
               for(int k=0;k<dp::num_macro_cells_z*2;k++){
                  if (k>= dp::num_macro_cells_z) kk = k - 2*dp::num_macro_cells_z;
                  else kk = k;
                  if((ii!=jj) && (jj != kk)){

                     const double rx = ii*cells::macro_cell_size[0]; // Angstroms
                     const double ry = jj*cells::macro_cell_size[1];
                     const double rz = kk*cells::macro_cell_size[2];
                  //   std::cout << "r" << rx << '\t' << ry << '\t' << rz << "\t" << cells::macro_cell_size << '\t' << rij << '\t' << rij3 << '\t' << ex << '\t' << ey << '\t' << ez <<std::endl;
                     const double rij = 1.0/pow(rx*rx+ry*ry+rz*rz,0.5);

                     const double ex = rx*rij;
                     const double ey = ry*rij;
                     const double ez = rz*rij;

                     const double rij3 = rij*rij*rij; // Angstroms
                  // std::cout << "r" << rx << '\t' << ry << '\t' << rz << "\t" << cell_size[0] << '\t' << rij << '\t' << rij3 << '\t' << ex << '\t' << ey << '\t' << ez <<std::endl;
                     int id = (i*dp::num_macro_cells_x+j)*dp::num_macro_cells_y+k;
                     dp::N2xx0[id][0] = (3.0*ex*ex - 1.0)*rij3;
                     dp::N2xx0[id][0] = (3.0*ex*ex - 1.0)*rij3;
                     dp::N2xy0[id][0] = (3.0*ex*ey      )*rij3;
                     dp::N2xz0[id][0] = (3.0*ex*ez      )*rij3;

                     dp::N2yx0[id][0] = (3.0*ey*ex - 1.0)*rij3;
                     dp::N2yy0[id][0] = (3.0*ey*ey      )*rij3;
                     dp::N2yz0[id][0] = (3.0*ey*ez      )*rij3;

                     dp::N2zx0[id][0] = (3.0*ez*ex - 1.0)*rij3;
                     dp::N2zy0[id][0] = (3.0*ez*ey      )*rij3;
                     dp::N2zz0[id][0] = (3.0*ez*ez      )*rij3;


            //   std::cout << 	i << '\t' << j << "\t" << k << '\t' << Nxx0[id][0] << '\t' << Nxy0[id][0] << '\t' << Nxz0[id][0] << '\t' << Nyy[id][0] << '\t' << Nyz0[id][0] << '\t' << Nzz0[id][0] <<std::endl;

                  }
               }
            }
         }


         // fft calculations
         fftw_plan NxxP,NxyP,NxzP,NyxP,NyyP,NyzP,NzxP,NzyP,NzzP;


          //deterines the forward transform for the N arrays
         NxxP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::N2xx0,dp::N2xx,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NxxP);
         NyxP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::N2yx0,dp::N2yx,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NyxP);
         NzxP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::N2zx0,dp::N2zx,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NzxP);
         NxyP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::N2xy0,dp::N2xy,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NxyP);
         NyyP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::N2yy0,dp::N2yy,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NyyP);
         NzyP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::N2zy0,dp::N2zy,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NzyP);
         NxzP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::N2xz0,dp::N2xz,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NxzP);
         NyzP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::N2yz0,dp::N2yz,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NyzP);
         NzzP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::N2zz0,dp::N2zz,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NzzP);

      }
     #endif
   }

   void dipole::internal::update_field_fft(){

   #ifdef FFT
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


   #endif
   }
}
