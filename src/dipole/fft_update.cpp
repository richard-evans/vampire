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
   // allocate arrays to store data [nloccell x ncells]
   if(dipole::fft==true) {

      // calculate matrix prefactors
      zlog << zTs() << "Precalculating rij matrix for dipole calculation... " << std::endl;

      // determine number of cells in each direction (with small shift to prevent the fence post problem)
      dp::num_macro_cells_x = static_cast<unsigned int>(ceil((cs::system_dimensions[0]+0.01)/cells::macro_cell_size[0]));
      dp::num_macro_cells_y = static_cast<unsigned int>(ceil((cs::system_dimensions[1]+0.01)/cells::macro_cell_size[1]));
      dp::num_macro_cells_z = static_cast<unsigned int>(ceil((cs::system_dimensions[2]+0.01)/cells::macro_cell_size[2]));

      //calcualtes 8 times number of cells
      dp::eight_num_cells = 8*dp::num_macro_cells_x*dp::num_macro_cells_y*dp::num_macro_cells_z;

      // Alocate memory for complex arrays to store the dipole tensor, magnetization and field arrays (in - N2xx0 and out N2xx etc.)
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

      //--------------------------------------------------------------------------------------
      // Initialise real and imaginary components of the dipole tensor to 0.0
      //--------------------------------------------------------------------------------------
      for(unsigned int i = 0 ; i < 2*dp::num_macro_cells_x ; i++){
         for(unsigned int j = 0 ; j < 2*dp::num_macro_cells_y; j++){
            for(unsigned int k = 0 ; k < 2*dp::num_macro_cells_z ; k++){

               int id = (i*dp::num_macro_cells_y*2 + j)*dp::num_macro_cells_z*2+k;

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

      //--------------------------------------------------------------------------------------
      // Calculate the dipole tensor between cells (assumes translational invariance for FFT)
      //--------------------------------------------------------------------------------------
      double ii,jj,kk;

      // perform a triple loop over all macrocells in x,y,z
      for(int i = 0; i < dp::num_macro_cells_x*2; i++){

         // apply periodic boundary conditions in x
         if (i >= dp::num_macro_cells_x) ii = i - 2*dp::num_macro_cells_x;
         else ii = i;

         for(int j = 0; j < dp::num_macro_cells_y*2; j++){

            // apply periodic boundary conditions in y
            if (j >= dp::num_macro_cells_y) jj = j - 2*dp::num_macro_cells_y;
            else jj = j;

            for(int k = 0; k < dp::num_macro_cells_z*2; k++){

               // apply periodic boundary conditions in z
               if (k >= dp::num_macro_cells_z) kk = k - 2*dp::num_macro_cells_z;
               else kk = k;

               // check that i != j != k to avoid NaN in sqrt()
               if( (ii != jj) && (jj != kk) ){

                  // calculate position vector to neighbouring cells
                  const double rx = double(ii) * cells::macro_cell_size[0]; // Angstroms
                  const double ry = double(jj) * cells::macro_cell_size[1];
                  const double rz = double(kk) * cells::macro_cell_size[2];

                  // calculate inverse distance
                  const double irij = 1.0/sqrt(rx*rx + ry*ry + rz*rz);

                  // calculate unit vector to cells
                  const double ex = rx * irij;
                  const double ey = ry * irij;
                  const double ez = rz * irij;

                  // calculate cube of interaction range
                  const double irij3 = irij * irij * irij; // Angstroms

                  // calculate index in 1D array
                  int id = (i*2*dp::num_macro_cells_y + j)*2*dp::num_macro_cells_z + k;

                  // set dipole tensor components
                  dp::N2xx0[id][0] = (3.0*ex*ex - 1.0)*irij3;
                  dp::N2xy0[id][0] = (3.0*ex*ey      )*irij3;
                  dp::N2xz0[id][0] = (3.0*ex*ez      )*irij3;

                  dp::N2yx0[id][0] = (3.0*ey*ex      )*irij3;
                  dp::N2yy0[id][0] = (3.0*ey*ey - 1.0)*irij3;
                  dp::N2yz0[id][0] = (3.0*ey*ez      )*irij3;

                  dp::N2zx0[id][0] = (3.0*ez*ex      )*irij3;
                  dp::N2zy0[id][0] = (3.0*ez*ey      )*irij3;
                  dp::N2zz0[id][0] = (3.0*ez*ez - 1.0)*irij3;

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

      NxxP = fftw_plan_dft_3d(nx, ny, nz, dp::N2xx0, dp::N2xx, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_execute(NxxP);

      NyxP = fftw_plan_dft_3d(nx, ny, nz, dp::N2yx0, dp::N2yx, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_execute(NyxP);

      NzxP = fftw_plan_dft_3d(nx, ny, nz, dp::N2zx0, dp::N2zx, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_execute(NzxP);

      NxyP = fftw_plan_dft_3d(nx, ny, nz, dp::N2xy0, dp::N2xy, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_execute(NxyP);

      NyyP = fftw_plan_dft_3d(nx, ny, nz, dp::N2yy0, dp::N2yy, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_execute(NyyP);

      NzyP = fftw_plan_dft_3d(nx, ny, nz, dp::N2zy0, dp::N2zy, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_execute(NzyP);

      NxzP = fftw_plan_dft_3d(nx, ny, nz, dp::N2xz0, dp::N2xz, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_execute(NxzP);

      NyzP = fftw_plan_dft_3d(nx, ny, nz, dp::N2yz0, dp::N2yz, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_execute(NyzP);

      NzzP = fftw_plan_dft_3d(nx, ny, nz, dp::N2zz0, dp::N2zz, FFTW_FORWARD, FFTW_ESTIMATE);
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

      zlog << zTs() << "dipole field calulation with FFT has been initalised " << std::endl;


   //--------------------------------------------------------------------------------------
   // Calulate FFT of dipole tensor FFT(N) -> result N2xx etc
   //--------------------------------------------------------------------------------------
   }
        #endif
   return;

}

//-----------------------------------------------------------------------------
// Function to update dipole fields using FFT solver
//-----------------------------------------------------------------------------
void update_field_fft(){

#ifdef FFT

   std::cerr << "doing fft update" << std::endl;

   //---------------------------------------------------------------
   // Initalise all the magnetization and field components to zero
   //---------------------------------------------------------------
   for (int i=0 ; i<2*dp::num_macro_cells_x ; i++){
      for (int j=0 ; j<2*dp::num_macro_cells_y ; j++){
         for (int k=0 ; k<2*dp::num_macro_cells_z ; k++){

            //int id = (i*dp::num_macro_cells_x+j)*dp::num_macro_cells_y+k;
            int id = (i*2*dp::num_macro_cells_y+j)*2*dp::num_macro_cells_z+k;

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
      }
   }

   //-----------------------------------------------
   // Set the magnetization inside the system
   //-----------------------------------------------

   // temporary counter for cell number
   int cell = 0;

   // inverse Bohr magneton to convert cell magnetization to J/T
   const double imuB = 1.0/9.27400915e-24;

   for (int i = 0; i < dp::num_macro_cells_x; i++){
      for (int j = 0; j < dp::num_macro_cells_y; j++){
         for (int k = 0; k < dp::num_macro_cells_z; k++){

            // Calculate cell id
            int id = (i*2*dp::num_macro_cells_y + j)*2*dp::num_macro_cells_z + k;

            Mx_in[id][0] = cells::mag_array_x[cell] * imuB;
            My_in[id][0] = cells::mag_array_y[cell] * imuB;
            Mz_in[id][0] = cells::mag_array_z[cell] * imuB;

            // increment cell counter
            cell ++;

         }
      }
   }

   fftw_plan MxP,MyP,MzP;

   //---------------------------------------------------------------------------
   // Calculate Fourier Transform of the magnetization FFT(M)
   //---------------------------------------------------------------------------

   MxP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::Mx_in,dp::Mx_out,FFTW_FORWARD,FFTW_ESTIMATE);
   fftw_execute(MxP);
   MyP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::My_in,dp::My_out,FFTW_FORWARD,FFTW_ESTIMATE);
   fftw_execute(MyP);
   MzP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::Mz_in,dp::Mz_out,FFTW_FORWARD,FFTW_ESTIMATE);
   fftw_execute(MzP);

   fftw_destroy_plan(MxP);
   fftw_destroy_plan(MyP);
   fftw_destroy_plan(MzP);

   //---------------------------------------------------------------------------
   // Perform the convolution between N and M [ FFT(N) . FFT(M) ]
   //---------------------------------------------------------------------------

   // loop over all cells in x,y,z
   for (int i = 0 ; i < 2*dp::num_macro_cells_x ; i++){
      for (int j = 0 ; j < 2*dp::num_macro_cells_y ; j++){
         for (int k = 0 ; k < 2*dp::num_macro_cells_z ; k++){

            // calculate array index
            int id = ((i*2*dp::num_macro_cells_y + j)*2*dp::num_macro_cells_z) + k;

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

         }
      }
   }

   //------------------------------------------------------------------------------------
   // Perform the backward transform to give the dipole field, H = iFFT( FFT(N).FFT(M) )
   //------------------------------------------------------------------------------------

   fftw_plan HxP,HyP,HzP;

   HxP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::Hx_in,dp::Hx_out,FFTW_BACKWARD,FFTW_ESTIMATE);
   fftw_execute(HxP);
   HyP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::Hy_in,dp::Hy_out,FFTW_BACKWARD,FFTW_ESTIMATE);
   fftw_execute(HyP);
   HzP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,dp::Hz_in,dp::Hz_out,FFTW_BACKWARD,FFTW_ESTIMATE);
   fftw_execute(HzP);

   fftw_destroy_plan(HxP);
   fftw_destroy_plan(HyP);
   fftw_destroy_plan(HzP);

   //-------------------------------------------------------------------------------------
   // loop over all local cells to initialise field with self term
   //-------------------------------------------------------------------------------------
   for(int lc = 0; lc < dipole::internal::cells_num_local_cells; lc++){

      // get cell index
      int i = cells::cell_id_array[lc];

      // check cell contains atoms
      if(dipole::internal::cells_num_atoms_in_cell[i]>0){

         // set constant for self dipole field
         const double eightPI_three_cell_volume = 8.0*M_PI/(3.0*dipole::internal::cells_volume_array[i]);
         const double self_demag = imuB * eightPI_three_cell_volume;

         // Add self-demagnetisation as mu_0/4_PI * 8PI/3V
         dipole::cells_field_array_x[i] = self_demag * cells::mag_array_x[i];
         dipole::cells_field_array_y[i] = self_demag * cells::mag_array_y[i];
         dipole::cells_field_array_z[i] = self_demag * cells::mag_array_z[i];

      }
   }

   //-------------------------------------------------------------------------------------
   // Add calculated field to field array
   //-------------------------------------------------------------------------------------

   // reset cell index counter
   cell = 0;

   // loop over all cells in x,y and z
   for( int i = 0 ; i < dp::num_macro_cells_x; i++){
      for( int j = 0 ; j < dp::num_macro_cells_y; j++){
         for( int k = 0 ; k < dp::num_macro_cells_z; k++){

            // calculate cell id
            int id = (i*2*dp::num_macro_cells_y + j)*2*dp::num_macro_cells_z + k;

            dipole::cells_field_array_x[cell] += Hx_out[id][0]/dp::eight_num_cells;
            dipole::cells_field_array_y[cell] += Hy_out[id][0]/dp::eight_num_cells;
            dipole::cells_field_array_z[cell] += Hz_out[id][0]/dp::eight_num_cells;

            dipole::cells_field_array_x[cell] *= 9.27400915e-01;
            dipole::cells_field_array_y[cell] *= 9.27400915e-01;
            dipole::cells_field_array_z[cell] *= 9.27400915e-01;

            // increment cell counter
            cell++;

         }
      }
   }

#endif

   std::cerr << "done" << std::endl;

   return;
}

//-----------------------------------------------------------------------------
// Function to finalize FFT solver and release memory
//-----------------------------------------------------------------------------
void finalize_fft_solver(){
#ifdef FFT
   // Print informative message to log file and screen
   std::cout << "Deallocating memory for FFT dipole calculation" << std::endl;
   zlog << zTs() << "Deallocating memory for FFT dipole calculation" << std::endl;

   // Free memory from FFT complex variables
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

#endif
std::cout << "A" <<std::endl;
   return;

}

} // end of internal namespace

} // end of dipole namespace
