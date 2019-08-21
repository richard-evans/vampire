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
#include <cstdlib>
#include <iostream>

// C library headers
#include <fenv.h>
#include <signal.h>

// Vampire headers
#include "cells.hpp" // needed for cells::cell_id_array but to be removed
#include "dipole.hpp"
#include "vio.hpp"
#include "vutil.hpp"
#include "sim.hpp"
#include "vutil.hpp"
// dipole module headers
#include "internal.hpp"
#include "../dipole/internal.hpp"


namespace hierarchical{

   namespace internal{



     void update(){
   //     std::cout << "A" << std::endl;

       const double imuB = 1.0/9.27400915e-24;
       hierarchical_mag();

       // instantiate timer
       vutil::vtimer_t timer;

       //   start timer
      //    timer.start();


    	for(int lc=0;lc<dipole::internal::cells_num_local_cells;lc++){
         int cell_i = cells::cell_id_array[lc];
      //   std::cout <<"C\t" << cell_i << '\t' << cells::mag_array_x[cell_i]*1.0/9.27400915e-24 <<'\t' << cells::mag_array_y[cell_i]*1.0/9.27400915e-24 <<'\t' << cells::mag_array_z[cell_i]*1.0/9.27400915e-24 <<std::endl;


         const int start = interaction_list_start_index[lc];
         const int end = interaction_list_end_index[lc];

         // Self demagnetisation factor multiplying m(i)
         const double V = dipole::internal::cells_volume_array[cell_i];
         const double eightPI_three_cell_volume = 8.0*M_PI/(3.0*V);
         double self_demag = eightPI_three_cell_volume;
         //std::cout << self_demag << std::endl;
         // Normalise cell magnetisation by the Bohr magneton
         const double mx_i = cells::mag_array_x[cell_i]*imuB;
         const double my_i = cells::mag_array_y[cell_i]*imuB;
         const double mz_i = cells::mag_array_z[cell_i]*imuB;
   //      std::cout << cell_i << '\t' << mx_i << '\t' << my_i << '\t' << mz_i << std::endl;
         // Add self-demagnetisation as mu_0/4_PI * 8PI*m_cell/3V
         dipole::cells_field_array_x[cell_i] = 0.0;//self_demag * mx_i*0.0; //*0.0
         dipole::cells_field_array_y[cell_i] = 0.0;//self_demag * my_i*0.0; //*0.0
         dipole::cells_field_array_z[cell_i] = 0.0;//self_demag * mz_i*0.0; //*0.0
         // Add self demag to Hdemag --> To get only dipole-dipole contribution comment this and initialise to zero
         dipole::cells_mu0Hd_field_array_x[cell_i] = -0.5*self_demag * mx_i;
         dipole::cells_mu0Hd_field_array_y[cell_i] = -0.5*self_demag * my_i;
         dipole::cells_mu0Hd_field_array_z[cell_i] = -0.5*self_demag * mz_i;

         for(int j = start; j<end;j++){
           int cell_j = interaction_list[j];

           const double mx = mag_array_x[cell_j]*imuB;
           const double my = mag_array_y[cell_j]*imuB;
           const double mz = mag_array_z[cell_j]*imuB;
         //if (cell_i == 0)std::cout<< cell_i << '\t' << mx_i << '\t' << my_i << '\t' << mz_i << "\t" <<  cell_j << '\t' << mx << '\t' << my << '\t' << mz <<std::endl;
           dipole::cells_field_array_x[cell_i]      +=(mx*rij_tensor_xx[j] + my*rij_tensor_xy[j] + mz*rij_tensor_xz[j]);
           dipole::cells_field_array_y[cell_i]      +=(mx*rij_tensor_xy[j] + my*rij_tensor_yy[j] + mz*rij_tensor_yz[j]);
           dipole::cells_field_array_z[cell_i]      +=(mx*rij_tensor_xz[j] + my*rij_tensor_yz[j] + mz*rij_tensor_zz[j]);
           // Demag field
          dipole::cells_mu0Hd_field_array_x[cell_i] +=(mx*rij_tensor_xx[j] + my*rij_tensor_xy[j] + mz*rij_tensor_xz[j]);
          dipole::cells_mu0Hd_field_array_y[cell_i] +=(mx*rij_tensor_xy[j] + my*rij_tensor_yy[j] + mz*rij_tensor_yz[j]);
          dipole::cells_mu0Hd_field_array_z[cell_i] +=(mx*rij_tensor_xz[j] + my*rij_tensor_yz[j] + mz*rij_tensor_zz[j]);
         // std::cout << rij_tensor_xx[j] << '\t' << rij_tensor_xy[j] << '\t' << rij_tensor_xz[j] << '\t' << rij_tensor_yy[j] << '\t' << rij_tensor_yz[j] << '\t' << rij_tensor_zz[j] << '\t' <<std::endl;
         }
      //   std::cout <<"D\t" << cell_i << '\t' << dipole::cells_field_array_x[cell_i] <<'\t' << dipole::cells_field_array_y[cell_i] <<'\t' << dipole::cells_field_array_z[cell_i] <<std::endl;

         dipole::cells_field_array_x[cell_i]       = dipole::cells_field_array_x[cell_i] * 9.27400915e-01;
         dipole::cells_field_array_y[cell_i]       = dipole::cells_field_array_y[cell_i] * 9.27400915e-01;
         dipole::cells_field_array_z[cell_i]       = dipole::cells_field_array_z[cell_i] * 9.27400915e-01;
         dp_fields << sim::time << '\t' << cell_i << '\t' << dipole::cells_field_array_x[cell_i] << '\t' << dipole::cells_field_array_y[cell_i] << '\t' << dipole::cells_field_array_z[cell_i] << '\t' << std::endl;

         // Multiply Hdemg by mu_0/4pi * 1e30 * mu_B to account for normalisation
         // of magnetisation and volume in angstrom
         dipole::cells_mu0Hd_field_array_x[cell_i] = dipole::cells_mu0Hd_field_array_x[cell_i] * 9.27400915e-01;
         dipole::cells_mu0Hd_field_array_y[cell_i] = dipole::cells_mu0Hd_field_array_y[cell_i] * 9.27400915e-01;
         dipole::cells_mu0Hd_field_array_z[cell_i] = dipole::cells_mu0Hd_field_array_z[cell_i] * 9.27400915e-01;
       //if (cell_i == 0) std::cout << sim::time << '\t' << cell_i << '\t' <<  dipole::cells_field_array_x[cell_i] << '\t' << dipole::cells_field_array_y[cell_i] << '\t' << dipole::cells_field_array_z[cell_i] << '\t' << std::endl;

      }
   //   timer.stop();
   //  std::cout << "\tdone! [ " << timer.elapsed_time() << " s ]" << std::endl;
    // zlog << zTs() <<  "\tDIPOLE UPDATE. Time taken: " << timer.elapsed_time() << " s"<< std::endl;







     }
   }
}
