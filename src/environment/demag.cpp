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
#include "atoms.hpp"
#include <math.h>
#include "cells.hpp"
#include "vio.hpp"
#ifdef FFT
#include <fftw3.h>
#endif
namespace environment{

   namespace internal{

      int initialise_demag_fields(){


        std::cout << "initialise demag env " <<std::endl;

        int num_interactions = num_cells*num_cells;
        rij_tensor_xx.resize(num_interactions,0.0);
        rij_tensor_xy.resize(num_interactions,0.0);
        rij_tensor_xz.resize(num_interactions,0.0);
        rij_tensor_yy.resize(num_interactions,0.0);
        rij_tensor_zz.resize(num_interactions,0.0);
        rij_tensor_yz.resize(num_interactions,0.0);
        int interaction_no = 0;
        for(int celli=0; celli<num_cells; celli++){

          for(int cellj=0; cellj<num_cells; cellj++){

              if (celli != cellj){

                double rx2 = cell_coords_array_x[cellj] - cell_coords_array_x[celli];
                double ry2 = cell_coords_array_y[cellj] - cell_coords_array_y[celli];
                double rz2 = cell_coords_array_z[cellj] - cell_coords_array_z[celli];
                double rij = sqrt(rx2*rx2+ry2*ry2+rz2*rz2);
                double rij_1 = 1.0/rij;//Reciprocal of the distance

                  // define unitarian distance vectors

                  const double ex = rx2*rij_1;
                  const double ey = ry2*rij_1;
                  const double ez = rz2*rij_1;

                  const double rij3 = (rij_1*rij_1*rij_1); // Angstroms

                  // calculate dipolar matrix for 6 entries because of symmetry
                  rij_tensor_xx[interaction_no] = ((3.0*ex*ex - 1.0)*rij3);
                  rij_tensor_xy[interaction_no] = ((3.0*ex*ey      )*rij3);
                  rij_tensor_xz[interaction_no] = ((3.0*ex*ez      )*rij3);
                  //
                  rij_tensor_yy[interaction_no] = ((3.0*ey*ey - 1.0)*rij3);
                  rij_tensor_yz[interaction_no] = ((3.0*ey*ez      )*rij3);
                  rij_tensor_zz[interaction_no] = ((3.0*ez*ez - 1.0)*rij3);
               //   std::cout << interaction_no << '\t' << cell_i << '\t' << cell_j << rij_tensor_xx[interaction_no] <<  std::endl;
                  interaction_no++;
              }
              else if (celli == cellj){

                rij_tensor_xx[interaction_no] = 0.0;
                rij_tensor_xy[interaction_no] = 0.0;
                rij_tensor_xz[interaction_no] = 0.0;
                //
                rij_tensor_yy[interaction_no] = 0.0;
                rij_tensor_yz[interaction_no] = 0.0;
                rij_tensor_zz[interaction_no] = 0.0;
                interaction_no++;

              }

          }
        }



         // std::cout << "Initialising dipole fields in environment module..." << std::endl;
         //
         // //only complie this section of FFT is enabled else don't
         // #ifdef FFT
         // //save eight times number cells and 8 pi/3V to use later.
         // eight_num_cells = 8*num_cells_x*num_cells_y*num_cells_z;
         // eightPI_three_cell_volume = 8.0*M_PI/(3.0*cell_volume);
         //
         // std::cout << num_cells_x << '\t' << num_cells_y << '\t' << num_cells_z <<std::endl;
         //
         // //Resize arrays for fft
         // N2xx =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2xy =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2xz =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2yx =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2yy =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2yz =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2zx =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2zy =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2zz =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         //
         // N2xx0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2xy0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2xz0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2yx0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2yy0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2yz0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2zx0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2zy0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // N2zz0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         //
         // Mx_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // My_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // Mz_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         //
         // Mx_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // My_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // Mz_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         //
         // Hx_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // Hy_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // Hz_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         //
         // Hx_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // Hy_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         // Hz_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * eight_num_cells);
         //
         //
         //
         // //initialises all the demag tensor components to 0
         // for (int id = 0; id < eight_num_cells; id++){
         //
         //
         //          N2xx0[id][0]=0;
         //          N2xx0[id][1]=0;
         //          N2xx[id][0] =0;
         //          N2xx[id][1] =0;
         //          N2xy0[id][0]=0;
         //          N2xy0[id][1]=0;
         //          N2xy[id][0] =0;
         //          N2xy[id][1] =0;
         //          N2xz0[id][0]=0;
         //          N2xz0[id][1]=0;
         //
         //          N2yx0[id][0]=0;
         //          N2yx0[id][1]=0;
         //          N2yx[id][0] =0;
         //          N2yx[id][1] =0;
         //          N2yy0[id][0]=0;
         //          N2yy0[id][1]=0;
         //          N2yy[id][0] =0;
         //          N2yy[id][1] =0;
         //          N2yz0[id][0]=0;
         //          N2yz0[id][1]=0;
         //
         //          N2zx0[id][0]=0;
         //          N2zx0[id][1]=0;
         //          N2zx[id][0] =0;
         //          N2zx[id][1] =0;
         //          N2zy0[id][0]=0;
         //          N2zy0[id][1]=0;
         //          N2zy[id][0] =0;
         //          N2zy[id][1] =0;
         //          N2zz0[id][0]=0;
         //          N2zz0[id][1]=0;
         //       }
         //
         //
         // //initalises all the non-zero demag tensor components as with the normal demag field calcualtion
         // for(unsigned int k = 0 ; k < 2*num_cells_z ; k++){
         //    for(unsigned int j = 0 ; j < 2*num_cells_y; j++){
         //       for(unsigned int i = 0 ; i < 2*num_cells_x ; i++){
         //          int id = k * 2*num_cells_x*2*num_cells_y + j * 2*num_cells_x + i;
         //
         //          int ii,jj,kk;
         //          if (i >= num_cells_x) ii = i - 2*num_cells_x;
         //          else ii = i;
         //          if (j >= num_cells_y) jj = j - 2*num_cells_y;
         //          else jj = j;
         //          if (k >= num_cells_z) kk = k - 2*num_cells_z;
         //          else kk = k;
         //
         //           const double rx = ii*cell_size_x[id]; // Angstroms
         //           const double ry = jj*cell_size_y[id];
         //           const double rz = kk*cell_size_z[id];
         //
         //           const double rij = sqrt(rx*rx + ry*ry + rz*rz);
         //           if (rij > 0.1){
         //           const double irij = 1.0/rij;
         //
         //           const double ex = rx*irij;
         //           const double ey = ry*irij;
         //           const double ez = rz*irij;
         //
         //           const double rij3 = irij*irij*irij; // Angstroms
         //           int id = (2*i*num_cells_y+j)*2*num_cells_z+k;
         //
         //           N2xx0[id][0] = (3.0*ex*ex - 1.0)*rij3;
         //           N2xy0[id][0] = (3.0*ex*ey      )*rij3;
         //           N2xz0[id][0] = (3.0*ex*ez      )*rij3;
         //
         //           N2yx0[id][0] = (3.0*ey*ex      )*rij3;
         //           N2yy0[id][0] = (3.0*ey*ey - 1.0)*rij3;
         //           N2yz0[id][0] = (3.0*ey*ez      )*rij3;
         //
         //           N2zx0[id][0] = (3.0*ez*ex      )*rij3;
         //           N2zy0[id][0] = (3.0*ez*ey      )*rij3;
         //           N2zz0[id][0] = (3.0*ez*ez - 1.0)*rij3;
         //      //     std::cout << ii << '\t' << jj << '\t' << kk << "\t" << N2xx0[id][0] << '\t' << N2xy0[id][0] << '\t' << N2xz0[id][0] << '\t' << N2yy0[id][0] << '\t' << N2yz0[id][0] << '\t' << N2zz0[id][0] << std::endl;
         //
         //          }
         //       }
         //     }
         //   }
         //
         // // fft calculations
         // fftw_plan NxxP,NxyP,NxzP,NyxP,NyyP,NyzP,NzxP,NzyP,NzzP;
         //
         // //deterines the forward transform for the demag field arrays from N2xx0 to N2xx arrays
         // NxxP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2xx0,N2xx,FFTW_FORWARD,FFTW_ESTIMATE);
         // fftw_execute(NxxP);
         // NyxP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2yx0,N2yx,FFTW_FORWARD,FFTW_ESTIMATE);
         // fftw_execute(NyxP);
         // NzxP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2zx0,N2zx,FFTW_FORWARD,FFTW_ESTIMATE);
         // fftw_execute(NzxP);
         // NxyP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2xy0,N2xy,FFTW_FORWARD,FFTW_ESTIMATE);
         // fftw_execute(NxyP);
         // NyyP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2yy0,N2yy,FFTW_FORWARD,FFTW_ESTIMATE);
         // fftw_execute(NyyP);
         // NzyP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2zy0,N2zy,FFTW_FORWARD,FFTW_ESTIMATE);
         // fftw_execute(NzyP);
         // NxzP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2xz0,N2xz,FFTW_FORWARD,FFTW_ESTIMATE);
         // fftw_execute(NxzP);
         // NyzP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2yz0,N2yz,FFTW_FORWARD,FFTW_ESTIMATE);
         // fftw_execute(NyzP);
         // NzzP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x, N2zz0,N2zz,FFTW_FORWARD,FFTW_ESTIMATE);
         // fftw_execute(NzzP);
         //
         // // free memory from FFTW plans
         // fftw_destroy_plan(NxxP);
         // fftw_destroy_plan(NxyP);
         // fftw_destroy_plan(NxzP);
         // fftw_destroy_plan(NyxP);
         // fftw_destroy_plan(NyyP);
         // fftw_destroy_plan(NyzP);
         // fftw_destroy_plan(NzxP);
         // fftw_destroy_plan(NzyP);
         // fftw_destroy_plan(NzzP);
         //
         //
         // MxP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x,Mx_in,Mx_out,FFTW_FORWARD,FFTW_ESTIMATE);
         // MyP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x,My_in,My_out,FFTW_FORWARD,FFTW_ESTIMATE);
         // MzP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x,Mz_in,Mz_out,FFTW_FORWARD,FFTW_ESTIMATE);
         //
         // HxP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x,Hx_in,Hx_out,FFTW_BACKWARD,FFTW_ESTIMATE);
         // HyP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x,Hy_in,Hy_out,FFTW_BACKWARD,FFTW_ESTIMATE);
         // HzP = fftw_plan_dft_3d(2*num_cells_z,2*num_cells_y,2*num_cells_x,Hz_in,Hz_out,FFTW_BACKWARD,FFTW_ESTIMATE);
         //
         // std::cout << "End of dipole fields initialisation..." << std::endl;
         //
         // #endif
         return 0;

      }

      int calculate_demag_fields(){

         // For MPI version, only add local atoms
          #ifdef MPICF
             int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
          #else
             int num_local_atoms = atoms::num_atoms;
          #endif

        int interaction_no = 0;
           const double imuB = 1.0/9.27400915e-24;

      //  std::cout << "calcualte demag env " <<std::endl;

      for(int cell_i = 0; cell_i<num_cells;cell_i++){

        //const double V = cell_volume[cell_i];
        //const double eightPI_three_cell_volume = 8.0*M_PI/(3.0*V);
        //double self_demag = eightPI_three_cell_volume;
        // Normalise cell magnetisation by the Bohr magneton

        //std::cout << cell_i << '\t' << mx_i << '\t' << my_i << '\t' << mz_i << std::endl;
        //Add self-demagnetisation as mu_0/4_PI * 8PI*m_cell/3V
        dipole_field_x[cell_i] = 0.0;//self_demag * mx_i*0.0; //*0.0
        dipole_field_y[cell_i] = 0.0;//self_demag * my_i*0.0; //*0.0
        dipole_field_z[cell_i] = 0.0;//self_demag * mz_i*0.0; //*0.0


        for(int cell_j = 0; cell_j<num_cells;cell_j++){
          int shield = shield_number[cell_j];
          //dont loop over micromagnetic cells to stop double couting - only environment cell_size
          //mm -> mm calculated in the main micromagnetic module
          if (shield != num_shields){
             int j = interaction_no;

             const double mx = x_mag_array[cell_j]*imuB;
             const double my = y_mag_array[cell_j]*imuB;
             const double mz = z_mag_array[cell_j]*imuB;
           //std::cout<< cell_i << '\t' << mx_i << '\t' << my_i << '\t' << mz_i << "\t" <<  cell_j << '\t' << mx << '\t' << my << '\t' << mz <<std::endl;
             dipole_field_x[cell_i]      +=(mx*rij_tensor_xx[j] + my*rij_tensor_xy[j] + mz*rij_tensor_xz[j]);
             dipole_field_y[cell_i]      +=(mx*rij_tensor_xy[j] + my*rij_tensor_yy[j] + mz*rij_tensor_yz[j]);
             dipole_field_z[cell_i]      +=(mx*rij_tensor_xz[j] + my*rij_tensor_yz[j] + mz*rij_tensor_zz[j]);
           // Demag field
            //std::cout << rij_tensor_xx[j] << '\t' << rij_tensor_xy[j] << '\t' << rij_tensor_xz[j] << '\t' << mx << std::endl;
            }
            interaction_no ++;
         }
        //std::cout <<"D\t" << cell_i << '\t' << dipole::cells_field_array_x[cell_i] <<'\t' << dipole::cells_field_array_y[cell_i] <<'\t' << dipole::cells_field_array_z[cell_i] <<std::endl;

        dipole_field_x[cell_i] = dipole_field_x[cell_i]*9.27400915e-01;
        dipole_field_y[cell_i] = dipole_field_y[cell_i]*9.27400915e-01;
        dipole_field_z[cell_i] = dipole_field_z[cell_i]*9.27400915e-01;
    //   std::cout << sim::time << cell_i << '\t' << dipole_field_x[cell_i] << '\t' << dipole_field_y[cell_i] << '\t' << dipole_field_z[cell_i] << '\t' << std::endl;
      }
    // std::cout << interaction_no << '\t' << num_cells*num_cells << std::endl;
    // std::cin.get();
    // saves the dipole field for each cell to the environment cell for use in the environment module
      for(int lc=0; lc<cells::num_local_cells; lc++){
         int cell = cells::cell_id_array[lc];
         int env_cell = list_env_cell_atomistic_cell[cell];
         //std::cout << cell << '\t' << env_cell <<std::endl;
         environment_field_x[cell] = dipole_field_x[env_cell];// + bias_field_x[env_cell];
         environment_field_y[cell] = dipole_field_y[env_cell];// + bias_field_y[env_cell];
         environment_field_z[cell] = dipole_field_z[env_cell];//+ bias_field_z[env_cell];

      //  std::cout << cells::pos_and_mom_array[4*cell+0] << '\t' << cells::pos_and_mom_array[4*cell+1] << '\t' << cells::pos_and_mom_array[4*cell+2] << '\t' << environment_field_x[cell] << '\t' << environment_field_y[cell] << '\t' <<environment_field_z[cell] << '\t' << std::endl;

      }
      //  std::cin.get();
//          #ifdef FFT
//          //initalise all components of M and H arrays to 0
//          for (int id = 0; id < eight_num_cells; id++){
//
//                   Mx_in[id][0]=0;
//                   Mx_in[id][1]=0;
//                   My_in[id][0]=0;
//                   My_in[id][1]=0;
//                   Mz_in[id][0]=0;
//                   Mz_in[id][1]=0;
//
//
//                   Mx_out[id][0]=0;
//                   Mx_out[id][1]=0;
//                   My_out[id][0]=0;
//                   My_out[id][1]=0;
//                   Mz_out[id][0]=0;
//                   Mz_out[id][1]=0;
//
//
//                   Hx_in[id][0]=0;
//                   Hx_in[id][1]=0;
//                   Hy_in[id][0]=0;
//                   Hy_in[id][1]=0;
//                   Hz_in[id][0]=0;
//                   Hz_in[id][1]=0;
//
//
//                   Hx_out[id][0]=0;
//                   Hx_out[id][1]=0;
//                   Hy_out[id][0]=0;
//                   Hy_out[id][1]=0;
//                   Hz_out[id][0]=0;
//                   Hz_out[id][1]=0;
//                }
//
//
//         // std::ofstream pfile;
//       //   pfile.open("m.txt");
//          //initialised the in components for the FT to the magnetisation of each cell
//          int cell = 0;
//          const double imuB = 1.0/9.27400915e-24;
//
//          for(unsigned int k = 0 ; k < num_cells_z ; k++){
//             for(unsigned int j = 0 ; j < num_cells_y; j++){
//                for(unsigned int i = 0 ; i < num_cells_x ; i++){
//                   int id = k * 2*num_cells_x*2*num_cells_y + j * 2*num_cells_x + i;
//
//                   // get cell index
//                   int cell = k * num_cells_x*num_cells_y + j * num_cells_x + i;
//
//                   Mx_in[id][0] = x_mag_array[cell] * imuB;
//                   My_in[id][0] = y_mag_array[cell] * imuB;
//                   Mz_in[id][0] = z_mag_array[cell] * imuB;
//               //    pfile <<i << '\t' << j << '\t' << k << "\t"<< cell << '\t' <<   Mx_in[id][0] << "\t" << My_in[id][0] << "\t" << Mz_in[id][0] << "\t" << std::endl;
//                 }
//               }
//             }
//
//
//          fftw_execute(MxP);
//          fftw_execute(MyP);
//          fftw_execute(MzP);
//
//
//          cell = 0;
//
//          // performs the converlusion between Nk and Mk
//          for(unsigned int k = 0 ; k < 2*num_cells_z ; k++){
//             for(unsigned int j = 0 ; j < 2*num_cells_y; j++){
//                for(unsigned int i = 0 ; i < 2*num_cells_x ; i++){
//                  int id = k * 2*num_cells_x*2*num_cells_y + j * 2*num_cells_x + i;
//
//                   Hx_in[id][0] = N2xx[id][0]*Mx_out[id][0] + N2xy[id][0]*My_out[id][0] + N2xz[id][0]*Mz_out[id][0]; //summing the real part
//                   Hx_in[id][0] -= (N2xx[id][1]*Mx_out[id][1] + N2xy[id][1]*My_out[id][1] + N2xz[id][1]*Mz_out[id][1]);
//
//                   Hx_in[id][1] = N2xx[id][0]*Mx_out[id][1] + N2xy[id][0]*My_out[id][1] + N2xz[id][0]*Mz_out[id][1];
//                   Hx_in[id][1] += (N2xx[id][1]*Mx_out[id][0] + N2xy[id][1]*My_out[id][0] + N2xz[id][1]*Mz_out[id][0]);
//
//                   Hy_in[id][0] = N2yx[id][0]*Mx_out[id][0] + N2yy[id][0]*My_out[id][0] + N2yz[id][0]*Mz_out[id][0];
//                   Hy_in[id][0] -= (N2yx[id][1]*Mx_out[id][1] + N2yy[id][1]*My_out[id][1] + N2yz[id][1]*Mz_out[id][1]);
//
//                   Hy_in[id][1] = N2yx[id][0]*Mx_out[id][1] + N2yy[id][0]*My_out[id][1] + N2yz[id][0]*Mz_out[id][1];
//                   Hy_in[id][1] += (N2yx[id][1]*Mx_out[id][0] + N2yy[id][1]*My_out[id][0] + N2yz[id][1]*Mz_out[id][0]);
//
//                   Hz_in[id][0] = N2zx[id][0]*Mx_out[id][0] + N2zy[id][0]*My_out[id][0] + N2zz[id][0]*Mz_out[id][0]; //summing the real part
//                   Hz_in[id][0] -= (N2zx[id][1]*Mx_out[id][1] + N2zy[id][1]*My_out[id][1] + N2zz[id][1]*Mz_out[id][1]);
//
//                   Hz_in[id][1] = N2zx[id][0]*Mx_out[id][1] + N2zy[id][0]*My_out[id][1] + N2zz[id][0]*Mz_out[id][1];
//                   Hz_in[id][1] += (N2zx[id][1]*Mx_out[id][0] + N2zy[id][1]*My_out[id][0] + N2zz[id][1]*Mz_out[id][0]);
//                   cell++;
//                }
//             }
//          }
//          fftw_execute(HxP);
//          fftw_execute(HyP);
//          fftw_execute(HzP);
//
//          double constant = 9.27400915e-01/eight_num_cells;
//
//
//          // for (int i = 0; i< num_cells; i++){
//          //
//          //    // Add self-demagnetisation as mu_0/4_PI * 8PI/3V
//          //    dipole_field_x[i]=0.0;//eightPI_three_cell_volume*(x_mag_array[i]/9.27400915e-24);
//          //    dipole_field_y[i]=0.0;//eightPI_three_cell_volume*(y_mag_array[i]/9.27400915e-24);
//          //    dipole_field_z[i]=0.0;//eightPI_three_cell_volume*(z_mag_array[i]/9.27400915e-24);
//          // }
//
//       //   std::ofstream ofile;
//       //   ofile.open("field.txt");
//
//          //sums the dipole field N.m + self demag/eightnumcells
//          for(unsigned int k = 0 ; k < num_cells_z ; k++){
//             for(unsigned int j = 0 ; j < num_cells_y; j++){
//                for(unsigned int i = 0 ; i < num_cells_x ; i++){
//                   int id = k * 2*num_cells_x*2*num_cells_y + j * 2*num_cells_x + i;
//
//
//                 // get cell index
//                 int cell = k * num_cells_x*num_cells_y + j * num_cells_x + i;
//
//                //   std:: cout << dipole_field_x[cell] << '\t' << dipole_field_y[cell] << '\t' << dipole_field_z[cell] << '\t' << Hx_out[id][0] << '\t' << Hy_out[id][0] << '\t' << Hz_out[id][0] << std::endl;
//                  dipole_field_x[cell] = Hx_out[id][0]*constant;
//                  dipole_field_y[cell] = Hy_out[id][0]*constant;
//                  dipole_field_z[cell] = Hz_out[id][0]*constant;
//               //   ofile <<i << '\t' << j << '\t' << k << "\t" <<  cell << '\t' << dipole_field_x[cell] << '\t' << dipole_field_y[cell] << '\t' << dipole_field_z[cell] << '\t' << std::endl;
//                }
//              }
//            }
// //std::cin.get();
  //  std::ofstream ofile;
  // ofile.open("field.txt");
//          //saves the dipole field for each cell to the environment cell for use in the environment module
         for(int lc=0; lc<cells::num_local_cells; lc++){
            int cell = cells::cell_id_array[lc];
            int env_cell = list_env_cell_atomistic_cell[cell];
            //std::cout << cell << '\t' << env_cell <<std::endl;
            environment_field_x[cell] = dipole_field_x[env_cell];// + bias_field_x[env_cell];
            environment_field_y[cell] = dipole_field_y[env_cell];// + bias_field_y[env_cell];
            environment_field_z[cell] = dipole_field_z[env_cell];// + bias_field_z[env_cell];

         //   std::cout << cells::pos_and_mom_array[4*cell+0] << '\t' << cells::pos_and_mom_array[4*cell+1] << '\t' << cells::pos_and_mom_array[4*cell+2] << '\t' << environment_field_x[cell] << '\t' << environment_field_y[cell] << '\t' <<environment_field_z[cell] << '\t' << std::endl;

         }
         for (int atom =0; atom < num_local_atoms; atom++){
            int cell = cells::atom_cell_id_array[atom];
            atomistic_environment_field_x[atom] = environment_field_x[cell];
            atomistic_environment_field_y[atom] = environment_field_y[cell];
            atomistic_environment_field_z[atom] = environment_field_z[cell];
         }

         return 0;

      }
   }
}
