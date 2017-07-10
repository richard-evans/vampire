// Vampire headers
#include "environment.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "sim.hpp"
#include <math.h>
#include "cells.hpp"
#include "array3d.h"
namespace environment{

     namespace internal{

       int initialise_demag_fields(){

       eight_num_cells = 8*num_cells_x*num_cells_y*num_cells_z;
       eightPI_three_cell_volume = 8.0*M_PI/(3.0*cell_volume);

       Nxx0.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
       Nxy0.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
       Nxz0.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);

       Nyx0.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
       Nyy0.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
       Nyz0.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);

       Nzx0.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
       Nzy0.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
       Nzz0.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);

       Nxx.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
       Nxy.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
       Nxz.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);

       Nyx.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
       Nyy.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
       Nyz.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);

       Nzx.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
       Nzy.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
       Nzz.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);


       Nxx.IFill(0.0);
       Nxy.IFill(0.0);
       Nxz.IFill(0.0);
       Nyx.IFill(0.0);
       Nyy.IFill(0.0);
       Nyz.IFill(0.0);
       Nzx.IFill(0.0);
       Nzy.IFill(0.0);
       Nzz.IFill(0.0);

       Nxx0.IFill(0.0);
       Nxy0.IFill(0.0);
       Nxz0.IFill(0.0);
       Nyx0.IFill(0.0);
       Nyy0.IFill(0.0);
       Nyz0.IFill(0.0);
       Nzx0.IFill(0.0);
       Nzy0.IFill(0.0);
       Nzz0.IFill(0.0);

      double ii,jj,kk;
         for(int i=0;i<num_cells_x*2;i++){
            if (i >= num_cells_x) ii = i - 2*num_cells_x;
            else ii = i;
            for(int j=0;j<num_cells_y*2;j++){
               if (j >= num_cells_y) jj = j - 2*num_cells_y;
               else jj = j;
               for(int k=0;k<num_cells_z*2;k++){
                  if (k>= num_cells_z) kk = k - 2*num_cells_z;
                  else kk = k;
                  if((ii!=jj) && (jj != kk)){

                     const double rx = ii*cell_size[0]; // Angstroms
                     const double ry = jj*cell_size[1];
                     const double rz = kk*cell_size[2];
                  //   std::cout << "r" << rx << '\t' << ry << '\t' << rz << "\t" << cells::macro_cell_size << '\t' << rij << '\t' << rij3 << '\t' << ex << '\t' << ey << '\t' << ez <<std::endl;
                     const double rij = 1.0/pow(rx*rx+ry*ry+rz*rz,0.5);

                     const double ex = rx*rij;
                     const double ey = ry*rij;
                     const double ez = rz*rij;

                     const double rij3 = rij*rij*rij; // Angstroms
                  // std::cout << "r" << rx << '\t' << ry << '\t' << rz << "\t" << cell_size[0] << '\t' << rij << '\t' << rij3 << '\t' << ex << '\t' << ey << '\t' << ez <<std::endl;

                     Nxx0(i,j,k)[0] = (3.0*ex*ex - 1.0)*rij3;
                     Nxy0(i,j,k)[0] = (3.0*ex*ey      )*rij3;
                     Nxz0(i,j,k)[0] = (3.0*ex*ez      )*rij3;

                     Nyx0(i,j,k)[0] = (3.0*ey*ex - 1.0)*rij3;
                     Nyy0(i,j,k)[0] = (3.0*ey*ey      )*rij3;
                     Nyz0(i,j,k)[0] = (3.0*ey*ez      )*rij3;

                     Nzx0(i,j,k)[0] = (3.0*ez*ex - 1.0)*rij3;
                     Nzy0(i,j,k)[0] = (3.0*ez*ey      )*rij3;
                     Nzz0(i,j,k)[0] = (3.0*ez*ez      )*rij3;


            //   std::cout << 	i << '\t' << j << "\t" << k << '\t' << Nxx0(i,j,k)[0] << '\t' << Nxy0(i,j,k)[0] << '\t' << Nxz0(i,j,k)[0] << '\t' << Nyy(i,j,k)[0] << '\t' << Nyz0(i,j,k)[0] << '\t' << Nzz0(i,j,k)[0] <<std::endl;

                  }
               }
            }
         }



         // fft calculations
         fftw_plan NxxP,NxyP,NxzP,NyxP,NyyP,NyzP,NzxP,NzyP,NzzP;
         int i = 3;
         int j = 2;
         int k = 3;
          std::cout << "INIT ENV" <<std::endl;
          std::cout << Nxx0(i,j,k)[0] << '\t' << Nxy0(i,j,k)[0] << '\t' << Nxz0(i,j,k)[0] << '\t' << Nyy0(i,j,k)[0] << '\t' << Nyz0(i,j,k)[0] << '\t' << Nzz0(i,j,k)[0] <<std::endl;
          std::cout << Nxx(i,j,k)[0] << '\t' << Nxy(i,j,k)[0] << '\t' << Nxz(i,j,k)[0] << '\t' << Nyy(i,j,k)[0] << '\t' << Nyz(i,j,k)[0] << '\t' << Nzz(i,j,k)[0] <<std::endl;
          std::cout << num_cells_x << '\t' << num_cells_y << '\t' << num_cells_z << std::endl;
          std::cout << *Nxx0.ptr() <<std::endl;
         //deterines the forward transform for the N arrays
         NxxP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Nxx0.ptr(),Nxx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NxxP);
         NyxP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Nyx0.ptr(),Nyx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NyxP);
         NzxP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Nzx0.ptr(),Nzx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NzxP);
         NxyP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Nxy0.ptr(),Nxy.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NxyP);
         NyyP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Nyy0.ptr(),Nyy.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NyyP);
         NzyP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Nzy0.ptr(),Nzy.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NzyP);
         NxzP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Nxz0.ptr(),Nxz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NxzP);
         NyzP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Nyz0.ptr(),Nyz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NyzP);
         NzzP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Nzz0.ptr(),Nzz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(NzzP);

          std::cout << "INIT ENV" <<std::endl;
          std::cout << Nxx0.getarrayelement(i,j,k) << Nxx0(i,j,k)[0] << '\t' << Nxy0(i,j,k)[0] << '\t' << Nxz0(i,j,k)[0] << '\t' << Nyy0(i,j,k)[0] << '\t' << Nyz0(i,j,k)[0] << '\t' << Nzz0(i,j,k)[0] <<std::endl;
          std::cout << Nxx(i,j,k)[0] << '\t' << Nxy(i,j,k)[0] << '\t' << Nxz(i,j,k)[0] << '\t' << Nyy(i,j,k)[0] << '\t' << Nyz(i,j,k)[0] << '\t' << Nzz(i,j,k)[0] <<std::endl;
          std::cin.get();
         return 0;

       }

        int calculate_demag_fields(){


         Array3D<fftw_complex> Mx; //3D Array for magneetisation
         Array3D<fftw_complex> My;
         Array3D<fftw_complex> Mz;

         Array3D<fftw_complex> Hx; //3D Array for dipolar field
         Array3D<fftw_complex> Hy;
         Array3D<fftw_complex> Hz;

         Array3D<fftw_complex> Mx2; //3D Array for magneetisation
         Array3D<fftw_complex> My2;
         Array3D<fftw_complex> Mz2;

         Array3D<fftw_complex> Hx2; //3D Array for dipolar field
         Array3D<fftw_complex> Hy2;
         Array3D<fftw_complex> Hz2;

         Hx.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
         Hy.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
         Hz.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);

         Mx.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
         My.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
         Mz.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);

         Hx2.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
         Hy2.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
         Hz2.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);

         Mx2.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
         My2.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);
         Mz2.resize(2*num_cells_x,2*num_cells_y,2*num_cells_z);


         Mx.IFill(0.0);
         My.IFill(0.0);
         Mz.IFill(0.0);

         Mx2.IFill(0.0);
         My2.IFill(0.0);
         Mz2.IFill(0.0);

         Hx.IFill(0.0);
         Hy.IFill(0.0);
         Hz.IFill(0.0);

         Hx2.IFill(0.0);
         Hy2.IFill(0.0);
         Hz2.IFill(0.0);


         int cell = 0;
         for (int i=0 ; i<num_cells_x; i++){
         	for (int j=0 ; j<num_cells_y; j++){
         		for (int k=0 ; k<num_cells_z; k++){

         			Mx(i,j,k)[0] = x_mag_array[cell]/9.27400915e-24;
         			My(i,j,k)[0] = y_mag_array[cell]/9.27400915e-24;
         			Mz(i,j,k)[0] = z_mag_array[cell]/9.27400915e-24;
               //   std::cout << x_mag_array[cell] << '\t' << y_mag_array[cell] << '\t' << z_mag_array[cell] << '\t' << Mx(i,j,k)[0] << '\t' << My(i,j,k)[0] << '\t' << Mz(i,j,k)[0] <<std::endl;
            	   cell ++;

         		}
         	 }
         }

         fftw_plan MxP,MyP,MzP;

         //std::cout << 'g' <<std::endl;
         MxP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Mx.ptr(),Mx2.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(MxP);
         MyP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,My.ptr(),My2.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(MyP);
         MzP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Mz.ptr(),Mz2.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(MzP);
         //std::cout << 'h' <<std::endl;

         cell = 0;

         // performs the converlusion between Nk and Mk
         for (int i=0 ; i<2*num_cells_x ; i++){
           for (int j=0 ; j<2*num_cells_y ; j++){
               for (int k=0 ; k<2*num_cells_z ; k++){
                        //  if (cell == 0) std::cout << "ANx" << Nxx(i,j,k)[0]  << "\t" <<   Nxy(i,j,k)[0]  << "\t" <<   Nxz(i,j,k)[0]  << "\t" << Nxx(i,j,k)[1]  << "\t" <<   Nxy(i,j,k)[1]  << "\t" <<   Nxz(i,j,k)[1]  << "\t" <<  std::endl;
                        //        if (cell == 0) std::cout << "ANy" << Nyx(i,j,k)[0]  << "\t" <<   Nyy(i,j,k)[0]  << "\t" <<   Nyz(i,j,k)[0]  << "\t" << Nyx(i,j,k)[1]  << "\t" <<   Nyy(i,j,k)[1]  << "\t" <<   Nyz(i,j,k)[1]  << "\t" <<  std::endl;
                        //        if (cell == 0) std::cout << "ANz" << Nzx(i,j,k)[0]  << "\t" <<   Nzy(i,j,k)[0]  << "\t" <<   Nzz(i,j,k)[0]  << "\t" << Nzx(i,j,k)[1]  << "\t" <<   Nzy(i,j,k)[1]  << "\t" <<   Nzz(i,j,k)[1]  << "\t" <<  std::endl;
                         //
                        //        if (cell == 0) std::cout << "AM" << Mx2(i,j,k)[0]  << "\t" <<   My2(i,j,k)[0]  << "\t" <<   Mz2(i,j,k)[0]  << "\t" << Mx2(i,j,k)[1]  << "\t" <<   My2(i,j,k)[1]  << "\t" <<   Mz2(i,j,k)[1]  << "\t" <<  std::endl;


                  Hx(i,j,k)[0] = Nxx(i,j,k)[0]*Mx2(i,j,k)[0] + Nxy(i,j,k)[0]*My2(i,j,k)[0] + Nxz(i,j,k)[0]*Mz2(i,j,k)[0]; //summing the real part
                  Hx(i,j,k)[0] -= (Nxx(i,j,k)[1]*Mx2(i,j,k)[1] + Nxy(i,j,k)[1]*My2(i,j,k)[1] + Nxz(i,j,k)[1]*Mz2(i,j,k)[1]);

                  Hx(i,j,k)[1] = Nxx(i,j,k)[0]*Mx2(i,j,k)[1] + Nxy(i,j,k)[0]*My2(i,j,k)[1] + Nxz(i,j,k)[0]*Mz2(i,j,k)[1];
                  Hx(i,j,k)[1] += (Nxx(i,j,k)[1]*Mx2(i,j,k)[0] + Nxy(i,j,k)[1]*My2(i,j,k)[0] + Nxz(i,j,k)[1]*Mz2(i,j,k)[0]);

                  Hy(i,j,k)[0] = Nyx(i,j,k)[0]*Mx2(i,j,k)[0] + Nyy(i,j,k)[0]*My2(i,j,k)[0] + Nyz(i,j,k)[0]*Mz2(i,j,k)[0];
                  Hy(i,j,k)[0] -= (Nyx(i,j,k)[1]*Mx2(i,j,k)[1] + Nyy(i,j,k)[1]*My2(i,j,k)[1] + Nyz(i,j,k)[1]*Mz2(i,j,k)[1]);

                  Hy(i,j,k)[1] = Nyx(i,j,k)[0]*Mx2(i,j,k)[1] + Nyy(i,j,k)[0]*My2(i,j,k)[1] + Nyz(i,j,k)[0]*Mz2(i,j,k)[1];
                  Hy(i,j,k)[1] += (Nyx(i,j,k)[1]*Mx2(i,j,k)[0] + Nyy(i,j,k)[1]*My2(i,j,k)[0] + Nyz(i,j,k)[1]*Mz2(i,j,k)[0]);

                  Hz(i,j,k)[0] = Nzx(i,j,k)[0]*Mx2(i,j,k)[0] + Nzy(i,j,k)[0]*My2(i,j,k)[0] + Nzz(i,j,k)[0]*Mz2(i,j,k)[0]; //summing the real part
                  Hz(i,j,k)[0] -= (Nzx(i,j,k)[1]*Mx2(i,j,k)[1] + Nzy(i,j,k)[1]*My2(i,j,k)[1] + Nzz(i,j,k)[1]*Mz2(i,j,k)[1]);

                  Hz(i,j,k)[1] = Nzx(i,j,k)[0]*Mx2(i,j,k)[1] + Nzy(i,j,k)[0]*My2(i,j,k)[1] + Nzz(i,j,k)[0]*Mz2(i,j,k)[1];
                  Hz(i,j,k)[1] += (Nzx(i,j,k)[1]*Mx2(i,j,k)[0] + Nzy(i,j,k)[1]*My2(i,j,k)[0] + Nzz(i,j,k)[1]*Mz2(i,j,k)[0]);
                  cell++;
               }
           }
         }
         //std::cin.get();
        // performs the backward transform to give the dipole field, Hx, Hy, Hz
        fftw_plan HxP,HyP,HzP;

        HxP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Hx.ptr(),Hx2.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
        fftw_execute(HxP);
        HyP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Hy.ptr(),Hy2.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
        fftw_execute(HyP);
        HzP = fftw_plan_dft_3d(2*num_cells_x,2*num_cells_y,2*num_cells_z,Hz.ptr(),Hz2.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
        fftw_execute(HzP);


	      for (int i = 0; i< num_cells; i++){

     //         	const double self_demag = demag::prefactor*eightPI_three_cell_volume;


           // Add self-demagnetisation as mu_0/4_PI * 8PI/3V
           dipole_field_x[i]=eightPI_three_cell_volume*(x_mag_array[i]/9.27400915e-24);
           dipole_field_y[i]=eightPI_three_cell_volume*(y_mag_array[i]/9.27400915e-24);
           dipole_field_z[i]=eightPI_three_cell_volume*(z_mag_array[i]/9.27400915e-24);

        }

        	cell = 0;
        	for (int i=0 ; i<num_cells_x ; i++){
        		for (int j=0 ; j<num_cells_y ; j++){
        			for (int k=0 ; k<num_cells_z ; k++){
        			  dipole_field_x[cell] += Hx2(i,j,k)[0]/eight_num_cells;
        			  dipole_field_y[cell] += Hy2(i,j,k)[0]/eight_num_cells;
        			  dipole_field_z[cell] += Hz2(i,j,k)[0]/eight_num_cells;
                 dipole_field_x[cell] *= 9.27400915e-01;
                 dipole_field_y[cell] *= 9.27400915e-01;
                 dipole_field_z[cell] *= 9.27400915e-01;
                // std::cout << "SELF" << eightPI_three_cell_volume*(x_mag_array[i]/9.27400915e-24) << '\t' << eightPI_three_cell_volume*(y_mag_array[i]/9.27400915e-24) << '\t' << eightPI_three_cell_volume*(z_mag_array[i]/9.27400915e-24) <<std::endl;
               //  std::cout << "OTHER" << Hx2(i,j,k)[0]/eight_num_cells << '\t' << Hy2(i,j,k)[0]/eight_num_cells << '\t' << Hz2(i,j,k)[0]/eight_num_cells <<std::endl;
               	cell++;
        			}
        		}
        	}
         for (int cell = 0; cell < cells::num_cells; cell++){

            int env_cell = list_env_cell_atomistic_cell[cell];
            environment_field_x[cell] = dipole_field_x[env_cell];
            environment_field_y[cell] = dipole_field_y[env_cell];
            environment_field_z[cell] = dipole_field_z[env_cell];
         //	std::cout << environment_field_x[cell] << "\t" << environment_field_z[cell] <<std::endl;
         }
   //      std::cin.get();
        	Hx.clear();
        	Hy.clear();
        	Hz.clear();
        	Mx.clear();
        	My.clear();
        	Mz.clear();

         return 0;

        }
     }
}
