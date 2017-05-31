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


	void dipole::internal::update_field_fft(){

      if(err::check==true){
   		terminaltextcolor(RED);
   		std::cerr << "demag::fft_update has been called " << vmpi::my_rank << std::endl;
   		terminaltextcolor(WHITE);
   	}

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

   	Hx.resize(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z);
      Hy.resize(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z);
      Hz.resize(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z);

      Mx.resize(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z);
      My.resize(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z);
      Mz.resize(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z);

   	Hx2.resize(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z);
   	Hy2.resize(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z);
   	Hz2.resize(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z);

   	Mx2.resize(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z);
   	My2.resize(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z);
   	Mz2.resize(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z);


         Mx.IFill(0.0);
         My.IFill(0.0);
         Mz.IFill(0.0);

   		Mx2.IFill(0.0);
         My2.IFill(0.0);
         Mz2.IFill(0.0);

   		Hx2.IFill(0.0);
         Hy2.IFill(0.0);
         Hz2.IFill(0.0);
   	//		std::cout << "a" <<std::endl;
   		int cell = 0;
   		for (int i=0 ; i<dp::num_macro_cells_x; i++){
   			for (int j=0 ; j<dp::num_macro_cells_y; j++){
   				for (int k=0 ; k<dp::num_macro_cells_z; k++){

   					Mx(i,j,k)[0] = cells::mag_array_x[cell]/9.27400915e-24;
   					My(i,j,k)[0] = cells::mag_array_y[cell]/9.27400915e-24;
   					Mz(i,j,k)[0] = cells::mag_array_z[cell]/9.27400915e-24;
               //   if (cell == 0) std::cout << Mx(i,j,k)[0]  << "\t" <<   My(i,j,k)[0]  << "\t" <<   Mz(i,j,k)[0]  << "\t" <<  std::endl;
 					   cell ++;

   				}
   			 }
   		}



         Hx.IFill(0.0);
         Hy.IFill(0.0);
         Hz.IFill(0.0);

   	   fftw_plan MxP,MyP,MzP;

   		//std::cout << 'g' <<std::endl;
         MxP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,Mx.ptr(),Mx2.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(MxP);
         MyP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,My.ptr(),My2.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(MyP);
         MzP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,Mz.ptr(),Mz2.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
         fftw_execute(MzP);
   		//std::cout << 'h' <<std::endl;

cell = 0;
         // performs the converlusion between Nk and Mk
      for (int i=0 ; i<2*dp::num_macro_cells_x ; i++){
         for (int j=0 ; j<2*dp::num_macro_cells_y ; j++){
             for (int k=0 ; k<2*dp::num_macro_cells_z ; k++){
         //       if (cell == 0) std::cout << "Nx" << Nxx(i,j,k)[0]  << "\t" <<   Nxy(i,j,k)[0]  << "\t" <<   Nxz(i,j,k)[0]  << "\t" << Nxx(i,j,k)[1]  << "\t" <<   Nxy(i,j,k)[1]  << "\t" <<   Nxz(i,j,k)[1]  << "\t" <<  std::endl;
         //       if (cell == 0) std::cout << "Ny" << Nyx(i,j,k)[0]  << "\t" <<   Nyy(i,j,k)[0]  << "\t" <<   Nyz(i,j,k)[0]  << "\t" << Nyx(i,j,k)[1]  << "\t" <<   Nyy(i,j,k)[1]  << "\t" <<   Nyz(i,j,k)[1]  << "\t" <<  std::endl;
         //       if (cell == 0) std::cout << "Nz" << Nzx(i,j,k)[0]  << "\t" <<   Nzy(i,j,k)[0]  << "\t" <<   Nzz(i,j,k)[0]  << "\t" << Nzx(i,j,k)[1]  << "\t" <<   Nzy(i,j,k)[1]  << "\t" <<   Nzz(i,j,k)[1]  << "\t" <<  std::endl;

            //    if (cell == 0) std::cout << "M" << Mx2(i,j,k)[0]  << "\t" <<   My2(i,j,k)[0]  << "\t" <<   Mz2(i,j,k)[0]  << "\t" << Mx2(i,j,k)[1]  << "\t" <<   My2(i,j,k)[1]  << "\t" <<   Mz2(i,j,k)[1]  << "\t" <<  std::endl;

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
   			//		 		 		std::cout << 	i << '\t' << j << "\t" << k << '\t' << Nxx(i,j,k)[0] << '\t' << Nxy(i,j,k)[0] << '\t' << Nxz(i,j,k)[0] << '\t' << Nyy(i,j,k)[0] << '\t' << Nyz(i,j,k)[0] << '\t' << Nzz(i,j,k)[0] <<std::endl;
    				//	 std::cout  << i << '\t' << j << '\t' << k << '\t' << Mx(i,j,k)[0] << '\t' << Mx(i,j,k)[1]<< "\t" << My(i,j,k)[0] << '\t' << My(i,j,k)[1]<< "\t" << Mz(i,j,k)[0] << '\t' << Mz(i,j,k)[1] << '\t' <<  Hx(i,j,k)[0] << '\t' << Hx(i,j,k)[1] << '\t' << Hy(i,j,k)[0] << '\t' << Hy(i,j,k)[1] << '\t'  << Hz(i,j,k)[0] << '\t' << Hz(i,j,k)[1] << '\t' << std::endl;

            //   if (cell == 0) std::cout << "H" << Hx(i,j,k)[0]  << "\t" <<   Hy(i,j,k)[0]  << "\t" <<   Hz(i,j,k)[0]  << "\t" << Hx(i,j,k)[1]  << "\t" <<   Hy(i,j,k)[1]  << "\t" <<   Hz(i,j,k)[1]  << "\t" <<  std::endl;

               cell++;
             }
         }
      }

//std::cin.get();

      // performs the backward transform to give the dipole field, Hx, Hy, Hz
      fftw_plan HxP,HyP,HzP;

      HxP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,Hx.ptr(),Hx2.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
      fftw_execute(HxP);
      HyP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,Hy.ptr(),Hy2.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
      fftw_execute(HyP);
      HzP = fftw_plan_dft_3d(2*dp::num_macro_cells_x,2*dp::num_macro_cells_y,2*dp::num_macro_cells_z,Hz.ptr(),Hz2.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
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
   				dipole::cells_field_array_x[cell] += Hx2(i,j,k)[0]/dp::eight_num_cells;
   				dipole::cells_field_array_y[cell] += Hy2(i,j,k)[0]/dp::eight_num_cells;
   				dipole::cells_field_array_z[cell] += Hz2(i,j,k)[0]/dp::eight_num_cells;
               dipole::cells_field_array_x[cell] *= 9.27400915e-01;
               dipole::cells_field_array_y[cell] *= 9.27400915e-01;
               dipole::cells_field_array_z[cell] *= 9.27400915e-01;
            //   if (cell ==0) std::cout << "cell field" << '\t' << sim::temperature << '\t'<< dipole::cells_field_array_x[cell] << '\t' << dipole::cells_field_array_y[cell] << '\t' << dipole::cells_field_array_z[cell] << '\t' << std::endl;

   		//	if (cell == 0)	std::cout << "fft" <<   '\t' << cell << '\t' << cells::num_cells << '\t' << dipole::cells_field_array_x[cell] <<  '\t' << dipole::cells_field_array_y[cell] << '\t' << dipole::cells_field_array_z[cell] << std::endl;

   				cell++;
   			}
   		}
   	}
   	Hx.clear();
   	Hy.clear();
   	Hz.clear();
   	Mx.clear();
   	My.clear();
   	Mz.clear();
   //std::cout << cells::x_field_array[cell] << '\t' << cells::y_field_array[cell] << '\t' << cells::z_field_array[cell] <<std::endl;

   }
}
