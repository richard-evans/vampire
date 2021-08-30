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

// alias internal hierarchical namespace for brevity
namespace ha = hierarchical::internal;

namespace hierarchical{

//------------------------------------------------------------------------------
// Function to update hierarchical cell magnetization and compute dipole field
//------------------------------------------------------------------------------
void update(std::vector <double>& x_spin_array, // atomic spin directions
            std::vector <double>& y_spin_array,
            std::vector <double>& z_spin_array,
            std::vector <double>& m_spin_array, // atomic spin moment
            std::vector < bool >& magnetic){ // is magnetic

   // inverse Bohr magneton
   const double imuB = 1.0/9.27400915e-24;

   // update hierarchical magnetization in cells
   hierarchical::internal::calculate_hierarchical_magnetisation(x_spin_array, y_spin_array, z_spin_array, m_spin_array, magnetic);

   // instantiate timer
   vutil::vtimer_t timer;

   // start timer
   timer.start();

   // initialise cell fields to a large negative number for MPI reduction
   for (int i = 0 ; i < cells::num_cells; i ++){
      dipole::cells_field_array_x[i] = -100000.0;
      dipole::cells_field_array_y[i] = -100000.0;
      dipole::cells_field_array_z[i] = -100000.0;
   }

   // Compute dipole fields for all cells with atoms (local cells)
	for(int lc = 0; lc < dipole::internal::cells_num_local_cells; lc++){

      // get global cell ID from local cell list
      int cell_i = cells::cell_id_array[lc];

      // Store range of hierarchical cells contributing to dipole field for lc cell
      const int start = ha::interaction_list_start_index[lc];
      const int end = ha::interaction_list_end_index[lc];

      // Self demagnetisation factor multiplying m(i)
      const double V = dipole::internal::cells_volume_array[cell_i];
      const double eightPI_three_cell_volume = 8.0*M_PI/(3.0*V);
      const double self_demag = eightPI_three_cell_volume;
      //std::cout << self_demag << std::endl;

      // Normalise cell magnetisation by the Bohr magneton
      const double mx_i = cells::mag_array_x[cell_i]*imuB;
      const double my_i = cells::mag_array_y[cell_i]*imuB;
      const double mz_i = cells::mag_array_z[cell_i]*imuB;
      // std::cout << cell_i << '\t' << mx_i << '\t' << my_i << '\t' << mz_i << std::endl;

      // Add self-demagnetisation as mu_0/4_PI * 8PI*m_cell/3V
      dipole::cells_field_array_x[cell_i] = self_demag * mx_i; //*0.0
      dipole::cells_field_array_y[cell_i] = self_demag * my_i; //*0.0
      dipole::cells_field_array_z[cell_i] = self_demag * mz_i; //*0.0

      // Add self demag to Hdemag --> To get only dipole-dipole contribution comment this and initialise to zero
      dipole::cells_mu0Hd_field_array_x[cell_i] = -0.5*self_demag * mx_i;
      dipole::cells_mu0Hd_field_array_y[cell_i] = -0.5*self_demag * my_i;
      dipole::cells_mu0Hd_field_array_z[cell_i] = -0.5*self_demag * mz_i;

      // Loop over all cells
      for(int j = start; j<end;j++){

         // get cell ID of neighbouring cell
         int cell_j = ha::interaction_list[j];

         const double mx = ha::mag_array_x[cell_j];//*imuB;
         const double my = ha::mag_array_y[cell_j];//*imuB;
         const double mz = ha::mag_array_z[cell_j];//*imuB;

         //    if (cell_i == 389)  std::cout << cell_j << '\t' << mx << '\t' << my << '\t' << mz<< "\t" << cell_i << '\t' << std::endl;
         //if (cell_i == 0)std::cout<< cell_i << '\t' << mx_i << '\t' << my_i << '\t' << mz_i << "\t" <<  cell_j << '\t' << mx << '\t' << my << '\t' << mz <<std::endl;

         // Compute dipole field contribution using dipole tensor
         dipole::cells_field_array_x[cell_i]      +=(mx*ha::rij_tensor_xx[j] + my*ha::rij_tensor_xy[j] + mz*ha::rij_tensor_xz[j]);
         dipole::cells_field_array_y[cell_i]      +=(mx*ha::rij_tensor_xy[j] + my*ha::rij_tensor_yy[j] + mz*ha::rij_tensor_yz[j]);
         dipole::cells_field_array_z[cell_i]      +=(mx*ha::rij_tensor_xz[j] + my*ha::rij_tensor_yz[j] + mz*ha::rij_tensor_zz[j]);

         // Demag field
         dipole::cells_mu0Hd_field_array_x[cell_i] +=(mx*ha::rij_tensor_xx[j] + my*ha::rij_tensor_xy[j] + mz*ha::rij_tensor_xz[j]);
         dipole::cells_mu0Hd_field_array_y[cell_i] +=(mx*ha::rij_tensor_xy[j] + my*ha::rij_tensor_yy[j] + mz*ha::rij_tensor_yz[j]);
         dipole::cells_mu0Hd_field_array_z[cell_i] +=(mx*ha::rij_tensor_xz[j] + my*ha::rij_tensor_yz[j] + mz*ha::rij_tensor_zz[j]);
         // std::cout << rij_tensor_xx[j] << '\t' << rij_tensor_xy[j] << '\t' << rij_tensor_xz[j] << '\t' << rij_tensor_yy[j] << '\t' << rij_tensor_yz[j] << '\t' << rij_tensor_zz[j] << '\t' <<std::endl;

      }
      //   std::cout <<"D\t" << cell_i << '\t' << dipole::cells_field_array_x[cell_i] <<'\t' << dipole::cells_field_array_y[cell_i] <<'\t' << dipole::cells_field_array_z[cell_i] <<std::endl;

      // Multiply Hdemg by mu_0/4pi * 1e30 * mu_B to account for normalisation of magnetisation and volume in angstrom
      dipole::cells_field_array_x[cell_i] = dipole::cells_field_array_x[cell_i] * 9.27400915e-01;
      dipole::cells_field_array_y[cell_i] = dipole::cells_field_array_y[cell_i] * 9.27400915e-01;
      dipole::cells_field_array_z[cell_i] = dipole::cells_field_array_z[cell_i] * 9.27400915e-01;
      //std::cout << sim::time << '\t' << cell_i << '\t' << dipole::cells_field_array_x[cell_i] << '\t' << dipole::cells_field_array_y[cell_i] << '\t' << dipole::cells_field_array_z[cell_i] << '\t' << std::endl;
      //std::cout << "SP" << dipole::cells_field_array_x[cell_i] << "\t" << dipole::cells_field_array_y[cell_i] << '\t' << dipole::cells_field_array_z[cell_i]<< std::endl;

      // Multiply Hdemg by mu_0/4pi * 1e30 * mu_B to account for normalisation of magnetisation and volume in angstrom
      dipole::cells_mu0Hd_field_array_x[cell_i] = dipole::cells_mu0Hd_field_array_x[cell_i] * 9.27400915e-01;
      dipole::cells_mu0Hd_field_array_y[cell_i] = dipole::cells_mu0Hd_field_array_y[cell_i] * 9.27400915e-01;
      dipole::cells_mu0Hd_field_array_z[cell_i] = dipole::cells_mu0Hd_field_array_z[cell_i] * 9.27400915e-01;
      //if (cell_i == 0) std::cout << sim::time << '\t' << cell_i << '\t' <<  dipole::cells_field_array_x[cell_i] << '\t' << dipole::cells_field_array_y[cell_i] << '\t' << dipole::cells_field_array_z[cell_i] << '\t' << std::endl;

   }

   #ifdef MPICF
      // Reduce fields on all processors so all have correct field values
      MPI_Allreduce(MPI_IN_PLACE, &dipole::cells_field_array_x[0], dipole::internal::cells_num_cells, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &dipole::cells_field_array_y[0], dipole::internal::cells_num_cells, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &dipole::cells_field_array_z[0], dipole::internal::cells_num_cells, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif

   // Check for cells with unrealistic fields from initialisation and zero
   for (int i = 0 ; i < dipole::internal::cells_num_cells; i ++){
      if (dipole::cells_field_array_x[i] < -1000) dipole::cells_field_array_x[i] = 0.0;
      if (dipole::cells_field_array_y[i] < -1000) dipole::cells_field_array_y[i] = 0.0;
      if (dipole::cells_field_array_z[i] < -1000) dipole::cells_field_array_z[i] = 0.0;
      //  std::cout << i << '\t' << dipole::cells_field_array_x[i] << '\t' << dipole::cells_field_array_y[i] << '\t' << dipole::cells_field_array_z[i] <<std::endl;
   }

   timer.stop();

   // Optionally output dipole update time to log file (commented due to call frequency of this function)
   //std::cout << "\tdone! [ " << timer.elapsed_time() << " s ]" << std::endl;
   //zlog << zTs() <<  "\tDIPOLE UPDATE. Time taken: " << timer.elapsed_time() << " s"<< std::endl;
   //std::cout << "dipole update time " << timer.elapsed_time() << " s" << std::endl;

   return;

} // end of update function

} // end of hierarchical namespace
