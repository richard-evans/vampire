//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//
// Vampire headers
#include "environment.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "sim.hpp"
#include <math.h>
#include "cells.hpp"
#ifdef FFT
#include <fftw3.h>
#endif
namespace environment{

   namespace internal{

      int initialise_demag_fields(){

         std::cout << "Initialising dipole fields in environment module..." << std::endl;

         //only complie this section of FFT is enabled else don't
         #ifdef FFT
         //save eight times number cells and 8 pi/3V to use later.
         eight_num_cells = 8*num_cells_x*num_cells_y*num_cells_z;
         eightPI_three_cell_volume = 8.0*M_PI/(3.0*cell_volume);

         std::cout << num_cells_x << '\t' << num_cells_y << '\t' << num_cells_z <<std::endl;

         //Resize arrays for fft
         N2xx =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2xy =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2xz =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2yx =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2yy =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2yz =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2zx =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2zy =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2zz =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);

         N2xx0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2xy0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2xz0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2yx0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2yy0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2yz0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2zx0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2zy0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         N2zz0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);

         Mx_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         My_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         Mz_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);

         Mx_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         My_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         Mz_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);

         Hx_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         Hy_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         Hz_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);

         Hx_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         Hy_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         Hz_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);



         //initialises all the demag tensor components to 0
         for (int id = 0; id < eight_num_cells; id++){


                  N2xx0[id][0]=0;
                  N2xx0[id][1]=0;
                  N2xx[id][0] =0;
                  N2xx[id][1] =0;
                  N2xy0[id][0]=0;
                  N2xy0[id][1]=0;
                  N2xy[id][0] =0;
                  N2xy[id][1] =0;
                  N2xz0[id][0]=0;
                  N2xz0[id][1]=0;

                  N2yx0[id][0]=0;
                  N2yx0[id][1]=0;
                  N2yx[id][0] =0;
                  N2yx[id][1] =0;
                  N2yy0[id][0]=0;
                  N2yy0[id][1]=0;
                  N2yy[id][0] =0;
                  N2yy[id][1] =0;
                  N2yz0[id][0]=0;
                  N2yz0[id][1]=0;

                  N2zx0[id][0]=0;
                  N2zx0[id][1]=0;
                  N2zx[id][0] =0;
                  N2zx[id][1] =0;
                  N2zy0[id][0]=0;
                  N2zy0[id][1]=0;
                  N2zy[id][0] =0;
                  N2zy[id][1] =0;
                  N2zz0[id][0]=0;
                  N2zz0[id][1]=0;
               }


         //initalises all the non-zero demag tensor components as with the normal demag field calcualtion
         for(unsigned int k = 0 ; k < 2*num_cells_z ; k++){
            for(unsigned int j = 0 ; j < 2*num_cells_y; j++){
               for(unsigned int i = 0 ; i < 2*num_cells_x ; i++){
                  int id = k * 2*num_cells_x*2*num_cells_y + j * 2*num_cells_x + i;

                  int ii,jj,kk;
                  if (i >= num_cells_x) ii = i - 2*num_cells_x;
                  else ii = i;
                  if (j >= num_cells_y) jj = j - 2*num_cells_y;
                  else jj = j;
                  if (k >= num_cells_z) kk = k - 2*num_cells_z;
                  else kk = k;

                   const double rx = ii*cell_size[0]; // Angstroms
                   const double ry = jj*cell_size[1];
                   const double rz = kk*cell_size[2];

                   const double rij = sqrt(rx*rx + ry*ry + rz*rz);
                   if (rij > 0.1){
                   const double irij = 1.0/rij;

                   const double ex = rx*irij;
                   const double ey = ry*irij;
                   const double ez = rz*irij;

                   const double rij3 = irij*irij*irij; // Angstroms
                   int id = (2*i*num_cells_y+j)*2*num_cells_z+k;

                   N2xx0[id][0] = (3.0*ex*ex - 1.0)*rij3;
                   N2xy0[id][0] = (3.0*ex*ey      )*rij3;
                   N2xz0[id][0] = (3.0*ex*ez      )*rij3;

                   N2yx0[id][0] = (3.0*ey*ex      )*rij3;
                   N2yy0[id][0] = (3.0*ey*ey - 1.0)*rij3;
                   N2yz0[id][0] = (3.0*ey*ez      )*rij3;

                   N2zx0[id][0] = (3.0*ez*ex      )*rij3;
                   N2zy0[id][0] = (3.0*ez*ey      )*rij3;
                   N2zz0[id][0] = (3.0*ez*ez - 1.0)*rij3;
              //     std::cout << ii << '\t' << jj << '\t' << kk << "\t" << N2xx0[id][0] << '\t' << N2xy0[id][0] << '\t' << N2xz0[id][0] << '\t' << N2yy0[id][0] << '\t' << N2yz0[id][0] << '\t' << N2zz0[id][0] << std::endl;

                  }
               }
             }
           }

         // fft calculations
         fftw_plan NxxP,NxyP,NxzP,NyxP,NyyP,NyzP,NzxP,NzyP,NzzP;

         //deterines the forward transform for the demag field arrays from N2xx0 to N2xx arrays
         NxxP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2xx0,N2xx,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NxxP);
         NyxP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2yx0,N2yx,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NyxP);
         NzxP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2zx0,N2zx,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NzxP);
         NxyP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2xy0,N2xy,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NxyP);
         NyyP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2yy0,N2yy,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NyyP);
         NzyP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2zy0,N2zy,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NzyP);
         NxzP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2xz0,N2xz,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NxzP);
         NyzP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2yz0,N2yz,FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NyzP);
         NzzP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2zz0,N2zz,FFTW_FORWARD,FFTW_ESTIMATE);
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


         MxP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x,Mx_in,Mx_out,FFTW_FORWARD,FFTW_ESTIMATE);
         MyP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x,My_in,My_out,FFTW_FORWARD,FFTW_ESTIMATE);
         MzP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x,Mz_in,Mz_out,FFTW_FORWARD,FFTW_ESTIMATE);

         HxP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x,Hx_in,Hx_out,FFTW_BACKWARD,FFTW_ESTIMATE);
         HyP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x,Hy_in,Hy_out,FFTW_BACKWARD,FFTW_ESTIMATE);
         HzP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x,Hz_in,Hz_out,FFTW_BACKWARD,FFTW_ESTIMATE);

         std::cout << "End of dipole fields initialisation..." << std::endl;

         #endif
         return 0;

      }

      int calculate_demag_fields(){

         #ifdef FFT
         //initalise all components of M and H arrays to 0
         for (int id = 0; id < eight_num_cells; id++){

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


        // std::ofstream pfile;
      //   pfile.open("m.txt");
         //initialised the in components for the FT to the magnetisation of each cell
         int cell = 0;
         const double imuB = 1.0/9.27400915e-24;

         for(unsigned int k = 0 ; k < num_cells_z ; k++){
            for(unsigned int j = 0 ; j < num_cells_y; j++){
               for(unsigned int i = 0 ; i < num_cells_x ; i++){
                  int id = k * 2*num_cells_x*2*num_cells_y + j * 2*num_cells_x + i;

                  // get cell index
                  int cell = k * num_cells_x*num_cells_y + j * num_cells_x + i;

                  Mx_in[id][0] = x_mag_array[cell] * imuB;
                  My_in[id][0] = y_mag_array[cell] * imuB;
                  Mz_in[id][0] = z_mag_array[cell] * imuB;
              //    pfile <<i << '\t' << j << '\t' << k << "\t"<< cell << '\t' <<   Mx_in[id][0] << "\t" << My_in[id][0] << "\t" << Mz_in[id][0] << "\t" << std::endl;
                }
              }
            }


         fftw_execute(MxP);
         fftw_execute(MyP);
         fftw_execute(MzP);


         cell = 0;

         // performs the converlusion between Nk and Mk
         for(unsigned int k = 0 ; k < 2*num_cells_z ; k++){
            for(unsigned int j = 0 ; j < 2*num_cells_y; j++){
               for(unsigned int i = 0 ; i < 2*num_cells_x ; i++){
                 int id = k * 2*num_cells_x*2*num_cells_y + j * 2*num_cells_x + i;

                  Hx_in[id][0] = N2xx[id][0]*Mx_out[id][0] + N2xy[id][0]*My_out[id][0] + N2xz[id][0]*Mz_out[id][0]; //summing the real part
                  Hx_in[id][0] -= (N2xx[id][1]*Mx_out[id][1] + N2xy[id][1]*My_out[id][1] + N2xz[id][1]*Mz_out[id][1]);

                  Hx_in[id][1] = N2xx[id][0]*Mx_out[id][1] + N2xy[id][0]*My_out[id][1] + N2xz[id][0]*Mz_out[id][1];
                  Hx_in[id][1] += (N2xx[id][1]*Mx_out[id][0] + N2xy[id][1]*My_out[id][0] + N2xz[id][1]*Mz_out[id][0]);

                  Hy_in[id][0] = N2yx[id][0]*Mx_out[id][0] + N2yy[id][0]*My_out[id][0] + N2yz[id][0]*Mz_out[id][0];
                  Hy_in[id][0] -= (N2yx[id][1]*Mx_out[id][1] + N2yy[id][1]*My_out[id][1] + N2yz[id][1]*Mz_out[id][1]);

                  Hy_in[id][1] = N2yx[id][0]*Mx_out[id][1] + N2yy[id][0]*My_out[id][1] + N2yz[id][0]*Mz_out[id][1];
                  Hy_in[id][1] += (N2yx[id][1]*Mx_out[id][0] + N2yy[id][1]*My_out[id][0] + N2yz[id][1]*Mz_out[id][0]);

                  Hz_in[id][0] = N2zx[id][0]*Mx_out[id][0] + N2zy[id][0]*My_out[id][0] + N2zz[id][0]*Mz_out[id][0]; //summing the real part
                  Hz_in[id][0] -= (N2zx[id][1]*Mx_out[id][1] + N2zy[id][1]*My_out[id][1] + N2zz[id][1]*Mz_out[id][1]);

                  Hz_in[id][1] = N2zx[id][0]*Mx_out[id][1] + N2zy[id][0]*My_out[id][1] + N2zz[id][0]*Mz_out[id][1];
                  Hz_in[id][1] += (N2zx[id][1]*Mx_out[id][0] + N2zy[id][1]*My_out[id][0] + N2zz[id][1]*Mz_out[id][0]);
                  cell++;
               }
            }
         }
         fftw_execute(HxP);
         fftw_execute(HyP);
         fftw_execute(HzP);

         double constant = 9.27400915e-01/eight_num_cells;


         // for (int i = 0; i< num_cells; i++){
         //
         //    // Add self-demagnetisation as mu_0/4_PI * 8PI/3V
         //    dipole_field_x[i]=0.0;//eightPI_three_cell_volume*(x_mag_array[i]/9.27400915e-24);
         //    dipole_field_y[i]=0.0;//eightPI_three_cell_volume*(y_mag_array[i]/9.27400915e-24);
         //    dipole_field_z[i]=0.0;//eightPI_three_cell_volume*(z_mag_array[i]/9.27400915e-24);
         // }

      //   std::ofstream ofile;
      //   ofile.open("field.txt");

         //sums the dipole field N.m + self demag/eightnumcells
         for(unsigned int k = 0 ; k < num_cells_z ; k++){
            for(unsigned int j = 0 ; j < num_cells_y; j++){
               for(unsigned int i = 0 ; i < num_cells_x ; i++){
                  int id = k * 2*num_cells_x*2*num_cells_y + j * 2*num_cells_x + i;


                // get cell index
                int cell = k * num_cells_x*num_cells_y + j * num_cells_x + i;

               //   std:: cout << dipole_field_x[cell] << '\t' << dipole_field_y[cell] << '\t' << dipole_field_z[cell] << '\t' << Hx_out[id][0] << '\t' << Hy_out[id][0] << '\t' << Hz_out[id][0] << std::endl;
                 dipole_field_x[cell] = Hx_out[id][0]*constant;
                 dipole_field_y[cell] = Hy_out[id][0]*constant;
                 dipole_field_z[cell] = Hz_out[id][0]*constant;
              //   ofile <<i << '\t' << j << '\t' << k << "\t" <<  cell << '\t' << dipole_field_x[cell] << '\t' << dipole_field_y[cell] << '\t' << dipole_field_z[cell] << '\t' << std::endl;
               }
             }
           }
//std::cin.get();
         //saves the dipole field for each cell to the environment cell for use in the environment module
         for(int lc=0; lc<cells::num_local_cells; lc++){
            int cell = cells::cell_id_array[lc];
            int env_cell = list_env_cell_atomistic_cell[cell];
            //std::cout << cell << '\t' << env_cell <<std::endl;
            environment_field_x[cell] = sbias_field_x[env_cell];
            environment_field_y[cell] = dipole_field_y[env_cell] + bias_field_y[env_cell];
            environment_field_z[cell] = dipole_field_z[env_cell] + bias_field_z[env_cell];

          //  ofile << cells::pos_and_mom_array[4*cell+0] << '\t' << cells::pos_and_mom_array[4*cell+1] << '\t' << cells::pos_and_mom_array[4*cell+2] << '\t' << environment_field_x[cell] << '\t' << environment_field_y[cell] << '\t' <<environment_field_z[cell] << '\t' << std::endl;

         }
   //      //end the FFT only compilation
         #endif
         return 0;

      }
   }
}
