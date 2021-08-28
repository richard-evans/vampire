//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins, Andrea Meo and Richard F L Evans 2017.
//       All rights reserved.
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

#ifdef FFT
#include <fftw3.h>
#endif

// dipole module headers
#include "internal.hpp"

// shorthand form for dipole module
namespace dp = dipole::internal;

namespace dipole{

namespace internal{


   //------------------------------------------------------------------------------
   // Calculation of the dipole field using a Fast Fourier Transform
   //------------------------------------------------------------------------------
   // The demagnetizing field of a body can be expressed as
   //
   //                                   H = - N . M
   //
   // where H is the dipole field, N is the demagnetizing tensor, and M is the
   // magnetization. Applying discrete fourier transforms (DFT) gives
   //
   //                            DFT(H) = -DFT(N) . DFT(M)
   //
   // By rearrangement we recover the demagnetizing field H:
   //
   //                          H = iDFT[ -DFT(N) . DFT(M) ]
   //
   // Where iDFT is an inverse discrtete fourier transform. If the interaction
   // tensor N can be assumed to be translationally invariant then we can
   // accelerate the calculation of the field with a Fast Fourier Transform (FFT):
   //
   //                          H = iFFT[ -FFT(N) . FFT(M) ]
   //
   // where FFT(N) is fixed and scales with n log n numercial operations where n
   // is the number of cells.
   //
   //------------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Function to initialise dipole field calculation using FFT solver
   //-----------------------------------------------------------------------------

void initialize_fft_solver(){

   #ifdef FFT
         // std::cout << "fft" << std::endl;
         // std::cin.get();
         // calculate matrix prefactors
         zlog << zTs() << "Precalculating rij matrix for dipole calculation... " << std::endl;


   // determine number of cells in x and y (global)
   cells::num_macro_cells_fft[0] = static_cast<unsigned int>(ceil((cs::system_dimensions[0]+0.01)/cells::macro_cell_size_x));
   cells::num_macro_cells_fft[1] = static_cast<unsigned int>(ceil((cs::system_dimensions[1]+0.01)/cells::macro_cell_size_y));
   cells::num_macro_cells_fft[2] = static_cast<unsigned int>(ceil((cs::system_dimensions[2]+0.01)/cells::macro_cell_size_z));


   // determine number of cells in each direction (with small shift to prevent the fence post problem)
   dp::num_macro_cells_x = cells::num_macro_cells_fft[0];
   dp::num_macro_cells_y = cells::num_macro_cells_fft[1];
   dp::num_macro_cells_z = cells::num_macro_cells_fft[2];


   dipole::cells_field_array_x.resize(cells_num_cells,0.0);
   dipole::cells_field_array_y.resize(cells_num_cells,0.0);
   dipole::cells_field_array_z.resize(cells_num_cells,0.0);

   // resize mu_0*Hd-field cells array
   dipole::cells_mu0Hd_field_array_x.resize(cells_num_cells,0.0);
   dipole::cells_mu0Hd_field_array_y.resize(cells_num_cells,0.0);
   dipole::cells_mu0Hd_field_array_z.resize(cells_num_cells,0.0);

   //calcualtes 8 times number of cells
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



  for(unsigned int k = 0 ; k < 2*dp::num_macro_cells_z ; k++){
     for(unsigned int j = 0 ; j < 2*dp::num_macro_cells_y; j++){
        for(unsigned int i = 0 ; i < 2*dp::num_macro_cells_x ; i++){
          int id = k * 2*dp::num_macro_cells_x*2*dp::num_macro_cells_y + j * 2*dp::num_macro_cells_x + i;
      //    std::cout << id << '\t' << id2 << std::endl;
           dp::N2xx0[id][0] = 0.0;
           dp::N2xx0[id][1] = 0.0;
           dp::N2xx[id][0]  = 0.0;
           dp::N2xx[id][1]  = 0.0;
           dp::N2xy0[id][0] = 0.0;
           dp::N2xy0[id][1] = 0.0;
           dp::N2xy[id][0]  = 0.0;
           dp::N2xy[id][1]  = 0.0;
           dp::N2xz0[id][0] = 0.0;
           dp::N2xz0[id][1] = 0.0;

           dp::N2yx0[id][0] = 0.0;
           dp::N2yx0[id][1] = 0.0;
           dp::N2yx[id][0]  = 0.0;
           dp::N2yx[id][1]  = 0.0;
           dp::N2yy0[id][0] = 0.0;
           dp::N2yy0[id][1] = 0.0;
           dp::N2yy[id][0]  = 0.0;
           dp::N2yy[id][1]  = 0.0;
           dp::N2yz0[id][0] = 0.0;
           dp::N2yz0[id][1] = 0.0;

           dp::N2zx0[id][0] = 0.0;
           dp::N2zx0[id][1] = 0.0;
           dp::N2zx[id][0]  = 0.0;
           dp::N2zx[id][1]  = 0.0;
           dp::N2zy0[id][0] = 0.0;
           dp::N2zy0[id][1] = 0.0;
           dp::N2zy[id][0]  = 0.0;
           dp::N2zy[id][1]  = 0.0;
           dp::N2zz0[id][0] = 0.0;
           dp::N2zz0[id][1] = 0.0;
         }
       }
     }

     for(unsigned int k = 0 ; k < 2*dp::num_macro_cells_z ; k++){
        for(unsigned int j = 0 ; j < 2*dp::num_macro_cells_y; j++){
           for(unsigned int i = 0 ; i < 2*dp::num_macro_cells_x ; i++){
             int id = k * 2*dp::num_macro_cells_x*2*dp::num_macro_cells_y + j * 2*dp::num_macro_cells_x + i;

             int ii,jj,kk;
             if (i >= dp::num_macro_cells_x) ii = i - 2*dp::num_macro_cells_x;
             else ii = i;
             if (j >= dp::num_macro_cells_y) jj = j - 2*dp::num_macro_cells_y;
             else jj = j;
             if (k >= dp::num_macro_cells_z) kk = k - 2*dp::num_macro_cells_z;
             else kk = k;

             // std::cout << "enter" << std::endl;
            const double rx = double(ii) * cells::macro_cell_size_x; // Angstroms
            const double ry = double(jj) * cells::macro_cell_size_y;
            const double rz = double(kk) * cells::macro_cell_size_z;

            // calculate inverse distance
            const double rij = sqrt(rx*rx + ry*ry + rz*rz);
            if (rij > 0.1){
            const double irij = 1.0/rij;
            // if ( ii == i && jj ==j && kk==k )  std::cout << rx << '\t' << ry << '\t' << rz << '\t' << irij << std::endl;
            // calculate unit vector to cells
            const double ex = rx * irij;
            const double ey = ry * irij;
            const double ez = rz * irij;

          // calculate cube of interaction range
          const double irij3 = irij * irij * irij; // Angstroms

          dp::N2xx0[id][0] = (3.0*ex*ex - 1.0)*irij3;
          dp::N2xy0[id][0] = (3.0*ex*ey      )*irij3;
          dp::N2xz0[id][0] = (3.0*ex*ez      )*irij3;

          dp::N2yx0[id][0] = (3.0*ey*ex      )*irij3;
          dp::N2yy0[id][0] = (3.0*ey*ey - 1.0)*irij3;
          dp::N2yz0[id][0] = (3.0*ey*ez      )*irij3;

          dp::N2zx0[id][0] = (3.0*ez*ex      )*irij3;
          dp::N2zy0[id][0] = (3.0*ez*ey      )*irij3;
          dp::N2zz0[id][0] = (3.0*ez*ez - 1.0)*irij3;
        //   if (dp::N2zz0[id][0] > 0)
      //  std::cout << ii << '\t' << jj << '\t' << kk << "\t" << dp::N2xx0[id][0] << '\t' << dp::N2xy0[id][0] << '\t' << dp::N2xz0[id][0] << '\t' << dp::N2yy0[id][0] << '\t' << dp::N2yz0[id][0] << '\t' << dp::N2zz0[id][0] << std::endl;
       }
    }
  }
}

  //--------------------------------------------------------------------------------------
  // Calculate FFT of dipole tensor FFT(N) -> result N2xx etc
  //--------------------------------------------------------------------------------------

  // Set temporary plans to calculate FFT of the dipole tensor FFT(N)
  fftw_plan NxxP, NxyP, NxzP, NyxP, NyyP, NyzP, NzxP, NzyP, NzzP;

  // temporary constants for FFT size for readability

  const unsigned int nx = 2*dp::num_macro_cells_x;
  const unsigned int ny = 2*dp::num_macro_cells_y;
  const unsigned int nz = 2*dp::num_macro_cells_z;

  NxxP = fftw_plan_dft_3d(nz, ny, nx, dp::N2xx0, dp::N2xx, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(NxxP);

  NyxP = fftw_plan_dft_3d(nz, ny, nx, dp::N2yx0, dp::N2yx, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(NyxP);

  NzxP = fftw_plan_dft_3d(nz, ny, nx, dp::N2zx0, dp::N2zx, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(NzxP);

  NxyP = fftw_plan_dft_3d(nz, ny, nx, dp::N2xy0, dp::N2xy, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(NxyP);

  NyyP = fftw_plan_dft_3d(nz, ny, nx, dp::N2yy0, dp::N2yy, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(NyyP);

  NzyP = fftw_plan_dft_3d(nz, ny, nx, dp::N2zy0, dp::N2zy, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(NzyP);

  NxzP = fftw_plan_dft_3d(nz, ny, nx, dp::N2xz0, dp::N2xz, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(NxzP);

  NyzP = fftw_plan_dft_3d(nz, ny, nx, dp::N2yz0, dp::N2yz, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(NyzP);

  NzzP = fftw_plan_dft_3d(nz, ny, nx, dp::N2zz0, dp::N2zz, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(NzzP);

  // free memory from FFTW plans
  fftw_destroy_plan(NxxP);
  fftw_destroy_plan(NxyP);
  fftw_destroy_plan(NxzP);
  fftw_destroy_plan(NyxP);
  fftw_destroy_plan(NyyP);
  fftw_destroy_plan(NyzP);
  fftw_destroy_plan(NzxP);
  fftw_destroy_plan(NzyP);
  fftw_destroy_plan(NzzP);

  dp::MxP = fftw_plan_dft_3d(2*dp::num_macro_cells_z,2*dp::num_macro_cells_y,2*dp::num_macro_cells_x,dp::Mx_in,dp::Mx_out,FFTW_FORWARD,FFTW_ESTIMATE);

  dp::MyP = fftw_plan_dft_3d(2*dp::num_macro_cells_z,2*dp::num_macro_cells_y,2*dp::num_macro_cells_x,dp::My_in,dp::My_out,FFTW_FORWARD,FFTW_ESTIMATE);

  dp::MzP = fftw_plan_dft_3d(2*dp::num_macro_cells_z,2*dp::num_macro_cells_y,2*dp::num_macro_cells_x,dp::Mz_in,dp::Mz_out,FFTW_FORWARD,FFTW_ESTIMATE);

  dp::HxP = fftw_plan_dft_3d(2*dp::num_macro_cells_z,2*dp::num_macro_cells_y,2*dp::num_macro_cells_x,dp::Hx_in,dp::Hx_out,FFTW_BACKWARD,FFTW_ESTIMATE);

  dp::HyP = fftw_plan_dft_3d(2*dp::num_macro_cells_z,2*dp::num_macro_cells_y,2*dp::num_macro_cells_x,dp::Hy_in,dp::Hy_out,FFTW_BACKWARD,FFTW_ESTIMATE);

  dp::HzP = fftw_plan_dft_3d(2*dp::num_macro_cells_z,2*dp::num_macro_cells_y,2*dp::num_macro_cells_x,dp::Hz_in,dp::Hz_out,FFTW_BACKWARD,FFTW_ESTIMATE);

  zlog << zTs() << "dipole field calulation with FFT has been initalised " << std::endl;


  //--------------------------------------------------------------------------------------
  // Calulate FFT of dipole tensor FFT(N) -> result N2xx etc
  //--------------------------------------------------------------------------------------

   #endif

}

void update_field_fft(){


     //std::cout << "update" << std::endl;
   #ifdef FFT

  for (int id = 0; id < dp::eight_num_cells; id++){

     Mx_in[id][0] = 0.0;
     Mx_in[id][1] = 0.0;
     My_in[id][0] = 0.0;
     My_in[id][1] = 0.0;
     Mz_in[id][0] = 0.0;
     Mz_in[id][1] = 0.0;


     Mx_out[id][0] = 0.0;
     Mx_out[id][1] = 0.0;
     My_out[id][0] = 0.0;
     My_out[id][1] = 0.0;
     Mz_out[id][0] = 0.0;
     Mz_out[id][1] = 0.0;


     Hx_in[id][0] = 0.0;
     Hx_in[id][1] = 0.0;
     Hy_in[id][0] = 0.0;
     Hy_in[id][1] = 0.0;
     Hz_in[id][0] = 0.0;
     Hz_in[id][1] = 0.0;


     Hx_out[id][0] = 0.0;
     Hx_out[id][1] = 0.0;
     Hy_out[id][0] = 0.0;
     Hy_out[id][1] = 0.0;
     Hz_out[id][0] = 0.0;
     Hz_out[id][1] = 0.0;

  }

   const double imuB = 1.0/9.27400915e-24;

   for(unsigned int k = 0 ; k < dp::num_macro_cells_z ; k++){
      for(unsigned int j = 0 ; j < dp::num_macro_cells_y; j++){
         for(unsigned int i = 0 ; i < dp::num_macro_cells_x ; i++){
           int id = k * 2*dp::num_macro_cells_x*2*dp::num_macro_cells_y + j * 2*dp::num_macro_cells_x + i;
           int cell = k * dp::num_macro_cells_x*dp::num_macro_cells_y + j * dp::num_macro_cells_x + i;;
    // get cell index
     //int cell = cells::cell_id_array[lc];
     //int id = cell;
   //  std::cout << id << '\t' << cells::mag_array_x[cell] * imuB << '\t' << cells::mag_array_y[cell] * imuB << '\t' << cells::mag_array_z[cell] * imuB <<std::endl;
          Mx_in[id][0] = cells::mag_array_x[cell] * imuB;
          My_in[id][0] = cells::mag_array_y[cell] * imuB;
          Mz_in[id][0] = cells::mag_array_z[cell] * imuB;
        //  std::cout << cell << '\t' << id << '\t' << i << '\t' << j << '\t' << k << "\t" << cells::mag_array_x[cell] * imuB << "\t" << cells::mag_array_y[cell] * imuB << "\t" << cells::mag_array_z[cell] * imuB  << "\t" << cells::pos_and_mom_array[4*cell+0] << '\t' << cells_pos_and_mom_array[4*cell+1] << '\t' << cells_pos_and_mom_array[4*cell+2] <<  std::endl;
        }
      }
    }


    fftw_execute(MxP);
    fftw_execute(MyP);
    fftw_execute(MzP);

    for(unsigned int k = 0 ; k < 2*dp::num_macro_cells_z ; k++){
       for(unsigned int j = 0 ; j < 2*dp::num_macro_cells_y; j++){
          for(unsigned int i = 0 ; i < 2*dp::num_macro_cells_x ; i++){
            int id = k * 2*dp::num_macro_cells_x*2*dp::num_macro_cells_y + j * 2*dp::num_macro_cells_x + i;


          //  int cell = id;//k * dp::num_macro_cells_x*dp::num_macro_cells_y + j * dp::num_macro_cells_x + i;
            //int id = k * 2*dp::num_macro_cells_x*2*dp::num_macro_cells_y + j * 2*dp::num_macro_cells_x + i;
            //int id2 = dp::idarray[i][j][k];
            Hx_in[id][0]  =  dp::N2xx[id][0] * dp::Mx_out[id][0] + dp::N2xy[id][0] * dp::My_out[id][0] + dp::N2xz[id][0] * dp::Mz_out[id][0]; //summing the real part
            Hx_in[id][0] -= (dp::N2xx[id][1] * dp::Mx_out[id][1] + dp::N2xy[id][1] * dp::My_out[id][1] + dp::N2xz[id][1] * dp::Mz_out[id][1]);

            Hx_in[id][1]  =  dp::N2xx[id][0] * dp::Mx_out[id][1] + dp::N2xy[id][0] * dp::My_out[id][1] + dp::N2xz[id][0] * dp::Mz_out[id][1];
            Hx_in[id][1] += (dp::N2xx[id][1] * dp::Mx_out[id][0] + dp::N2xy[id][1] * dp::My_out[id][0] + dp::N2xz[id][1] * dp::Mz_out[id][0]);

            Hy_in[id][0]  =  dp::N2yx[id][0] * dp::Mx_out[id][0] + dp::N2yy[id][0] * dp::My_out[id][0] + dp::N2yz[id][0] * dp::Mz_out[id][0];
            Hy_in[id][0] -= (dp::N2yx[id][1] * dp::Mx_out[id][1] + dp::N2yy[id][1] * dp::My_out[id][1] + dp::N2yz[id][1] * dp::Mz_out[id][1]);

            Hy_in[id][1]  =  dp::N2yx[id][0] * dp::Mx_out[id][1] + dp::N2yy[id][0] * dp::My_out[id][1] + dp::N2yz[id][0] * dp::Mz_out[id][1];
            Hy_in[id][1] += (dp::N2yx[id][1] * dp::Mx_out[id][0] + dp::N2yy[id][1] * dp::My_out[id][0] + dp::N2yz[id][1] * dp::Mz_out[id][0]);

            Hz_in[id][0]  =  dp::N2zx[id][0] * dp::Mx_out[id][0] + dp::N2zy[id][0] * dp::My_out[id][0] + dp::N2zz[id][0] * dp::Mz_out[id][0]; //summing the real part
            Hz_in[id][0] -= (dp::N2zx[id][1] * dp::Mx_out[id][1] + dp::N2zy[id][1] * dp::My_out[id][1] + dp::N2zz[id][1] * dp::Mz_out[id][1]);

            Hz_in[id][1]  =  dp::N2zx[id][0] * dp::Mx_out[id][1] + dp::N2zy[id][0] * dp::My_out[id][1] + dp::N2zz[id][0] * dp::Mz_out[id][1];
            Hz_in[id][1] += (dp::N2zx[id][1] * dp::Mx_out[id][0] + dp::N2zy[id][1] * dp::My_out[id][0] + dp::N2zz[id][1] * dp::Mz_out[id][0]);
      //      std::cout <<dp::Mx_out[id][0] << '\t' << dp::My_out[id][0] << '\t' << dp::Mz_out[id][0] <<"\t" <<   Hx_in[id][0] << '\t' <<  Hy_in[id][0] << '\t' <<  Hy_in[id][0] << "\t" << cells::pos_and_mom_array[4*cell+0] << '\t' << cells_pos_and_mom_array[4*cell+1] << '\t' << cells_pos_and_mom_array[4*cell+2] << std::endl;
          }
        }
    }


   fftw_execute(HxP);
   fftw_execute(HyP);
   fftw_execute(HzP);


   for(unsigned int k = 0 ; k < dp::num_macro_cells_z ; k++){
      for(unsigned int j = 0 ; j < dp::num_macro_cells_y; j++){
         for(unsigned int i = 0 ; i < dp::num_macro_cells_x ; i++){
           int id = k * 2*dp::num_macro_cells_x*2*dp::num_macro_cells_y + j * 2*dp::num_macro_cells_x + i;
           int cell = k * dp::num_macro_cells_x*dp::num_macro_cells_y + j * dp::num_macro_cells_x + i;;
     //int cell = k * dp::num_macro_cells_x*dp::num_macro_cells_y + j * dp::num_macro_cells_x + i;
    // int id = k * 2*dp::num_macro_cells_x*2*dp::num_macro_cells_y + j * 2*dp::num_macro_cells_x + i;
    // int id2 = dp::idarray[i][j][k];
     dipole::cells_field_array_x[cell] = Hx_out[id][0]/dp::eight_num_cells;
     dipole::cells_field_array_y[cell] = Hy_out[id][0]/dp::eight_num_cells;
     dipole::cells_field_array_z[cell] = Hz_out[id][0]/dp::eight_num_cells;

     dipole::cells_field_array_x[cell] *= 9.27400915e-01;
     dipole::cells_field_array_y[cell] *= 9.27400915e-01;
     dipole::cells_field_array_z[cell] *= 9.27400915e-01;

    //dp_fields <<sim::time<< '\t' <<  dipole::cells_field_array_x[cell] << '\t' << dipole::cells_field_array_y[cell] << '\t' << dipole::cells_field_array_z[cell]  << "\t" << cells::pos_and_mom_array[4*cell+0] << '\t' << cells_pos_and_mom_array[4*cell+1] << '\t' << cells_pos_and_mom_array[4*cell+2] << std::endl;
    cell++;

  }
}

  }
//std::cin.get();
   #endif

   return;

}

   //-----------------------------------------------------------------------------
   // Function to finalize FFT solver and release memory
   //-----------------------------------------------------------------------------
   void finalize_fft_solver(){
   #ifdef FFT
      // // Print informative message to log file and screen
      std::cout << "Deallocating memory for FFT dipole calculation" << std::endl;
      zlog << zTs() << "Deallocating memory for FFT dipole calculation" << std::endl;
      //
      // // Free memory from FFT complex variables
      fftw_free(dp::Mx_out);
      fftw_free(dp::My_out);
      fftw_free(dp::Mz_out);
      fftw_free(dp::Hx_out);
      fftw_free(dp::Hy_out);
      fftw_free(dp::Hz_out);
      fftw_free(dp::Mx_in);
      fftw_free(dp::My_in);
      fftw_free(dp::Mz_in);
      fftw_free(dp::Hx_in);
      fftw_free(dp::Hy_in);
      fftw_free(dp::Hz_in);

      fftw_free(dp::N2xx );
      fftw_free(dp::N2xy );
      fftw_free(dp::N2xz );
      fftw_free(dp::N2yx );
      fftw_free(dp::N2yy );
      fftw_free(dp::N2yz );
      fftw_free(dp::N2zx );
      fftw_free(dp::N2zy );
      fftw_free(dp::N2zz );
      fftw_free(dp::N2xx0);
      fftw_free(dp::N2xy0);
      fftw_free(dp::N2xz0);
      fftw_free(dp::N2yx0);
      fftw_free(dp::N2yy0);
      fftw_free(dp::N2yz0);
      fftw_free(dp::N2zx0);
      fftw_free(dp::N2zy0);
      fftw_free(dp::N2zz0);

      fftw_destroy_plan(dp::MxP);
      fftw_destroy_plan(dp::MyP);
      fftw_destroy_plan(dp::MzP);

      fftw_destroy_plan(dp::HxP);
      fftw_destroy_plan(dp::HyP);
      fftw_destroy_plan(dp::HzP);

   #endif
      return;

   }
   }
   }
