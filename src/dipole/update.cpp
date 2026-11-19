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
#include "constants.hpp"
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>


// Vampire headers
#include "dipole.hpp"
#include "micromagnetic.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "cells.hpp"
#include "sim.hpp"
#include "errors.hpp"

// dipole module headers
#include "internal.hpp"

namespace dipole{

   //-----------------------------------------------------------------------------
   // Function for updating dipolar / demag fields
   //-----------------------------------------------------------------------------


   void dipole::internal::update_field(){

      if(!dipole::activated) return;

      if(err::check==true){
         terminaltextcolor(RED);
         std::cerr << "dipole::update has been called " << vmpi::my_rank << std::endl;
         terminaltextcolor(WHITE);
      }

      // Define constant imuB = 1/muB to normalise to unitarian values the cell magnetisation
      const double imuB = constants::i_muB;

      // loop over local cells
      for(int lc=0;lc<dipole::internal::cells_num_local_cells;lc++){

         int i = cells::cell_id_array[lc];
         // cells without magnetic atoms (cells:probe-non-magnetic-cells) receive the field but have no volume or self term
         if(dipole::internal::cells_num_atoms_in_cell[i]>0 || cells::probe_non_magnetic_cells){

            const bool magnetic = dipole::internal::cells_num_atoms_in_cell[i]>0;

            // Self demagnetisation factor multiplying m(i)
            const double eightPI_three_cell_volume = magnetic ? 8.0*M_PI/(3.0*dipole::internal::cells_volume_array[i]) : 0.0;
            double self_demag = eightPI_three_cell_volume;

            // Normalise cell magnetisation by the Bohr magneton
            const double mx_i = magnetic ? cells::mag_array_x[i]*imuB : 0.0;
            const double my_i = magnetic ? cells::mag_array_y[i]*imuB : 0.0;
            const double mz_i = magnetic ? cells::mag_array_z[i]*imuB : 0.0;

            // Add self-demagnetisation as mu_0/4_PI * 8PI*m_cell/3V
            dipole::cells_field_array_x[i] = self_demag * mx_i;//*0.0; //*0.0
            dipole::cells_field_array_y[i] = self_demag * my_i;//*0.0; //*0.0
            dipole::cells_field_array_z[i] = self_demag * mz_i;//*0.0; //*0.0
            // Add self demag to Hdemag --> To get only dipole-dipole contribution comment this and initialise to zero
            dipole::cells_mu0Hd_field_array_x[i] = -0.5*self_demag * mx_i;
            dipole::cells_mu0Hd_field_array_y[i] = -0.5*self_demag * my_i;
            dipole::cells_mu0Hd_field_array_z[i] = -0.5*self_demag * mz_i;

            // tensor row of this cell: stored, or computed now with dipole:on-the-fly-tensors
            int row = lc;
            if(internal::on_the_fly_tensors){
               row = 0;
               internal::compute_tensor_row(internal::cells_local_cell_array[lc], row, internal::otf_cutoff, internal::otf_global_atoms_in_cell_count,
                                            internal::otf_cells_pos_and_mom_array, internal::otf_list_of_cells_with_atoms, internal::otf_atoms_in_cells_array);
            }

            // Loop over all other cells to calculate contribution to local cell
            for(int j=0;j<dipole::internal::cells_num_cells;j++){
               if(dipole::internal::cells_num_atoms_in_cell[j]>0){

                  // Normalise the cell magnetisation by the Bohr magneton
                  const double mx = cells::mag_array_x[j]*imuB;
                  const double my = cells::mag_array_y[j]*imuB;
                  const double mz = cells::mag_array_z[j]*imuB;
               // if (i == 0)std::cout<< i << '\t' << mx_i << '\t' << my_i << '\t' << mz_i << "\t" <<  j << '\t' << mx << '\t' << my << '\t' << mz <<std::endl;

                  dipole::cells_field_array_x[i]+=(mx*internal::rij_tensor_xx[row][j] + my*internal::rij_tensor_xy[row][j] + mz*internal::rij_tensor_xz[row][j]);
                  dipole::cells_field_array_y[i]+=(mx*internal::rij_tensor_xy[row][j] + my*internal::rij_tensor_yy[row][j] + mz*internal::rij_tensor_yz[row][j]);
                  dipole::cells_field_array_z[i]+=(mx*internal::rij_tensor_xz[row][j] + my*internal::rij_tensor_yz[row][j] + mz*internal::rij_tensor_zz[row][j]);
                  // Demag field
               // std::cout << "DEMAG\t" << cells::mag_array_x[j] << '\t' << imuB << '\t' <<  mx << '\t' << my << '\t' << mz << '\t' << internal::rij_tensor_xx[lc][j] << '\t' << dipole::cells_field_array_x[i] << std::endl;
                  dipole::cells_mu0Hd_field_array_x[i] +=(mx*internal::rij_tensor_xx[row][j] + my*internal::rij_tensor_xy[row][j] + mz*internal::rij_tensor_xz[row][j]);
                  dipole::cells_mu0Hd_field_array_y[i] +=(mx*internal::rij_tensor_xy[row][j] + my*internal::rij_tensor_yy[row][j] + mz*internal::rij_tensor_yz[row][j]);
                  dipole::cells_mu0Hd_field_array_z[i] +=(mx*internal::rij_tensor_xz[row][j] + my*internal::rij_tensor_yz[row][j] + mz*internal::rij_tensor_zz[row][j]);
               }
            }
            // Multiply the cells B-field by mu_B * mu_0/(4*pi) /1e-30  <-- (constants::muB * 1e-7 / 1e30)
            // where the last term accounts for the fact that the volume was calculated in Angstrom
            dipole::cells_field_array_x[i] = dipole::cells_field_array_x[i] * constants::dipole_prefactor;
            dipole::cells_field_array_y[i] = dipole::cells_field_array_y[i] * constants::dipole_prefactor;
            dipole::cells_field_array_z[i] = dipole::cells_field_array_z[i] * constants::dipole_prefactor;
            //dp_fields << sim::time << '\t' << i << '\t' << dipole::cells_field_array_x[i] << '\t' << dipole::cells_field_array_y[i] << '\t' << dipole::cells_field_array_z[i] << '\t' << std::endl;
            //dp_fields <<sim::time << "\t" <<  dipole::cells_field_array_x[i] << '\t' << dipole::cells_field_array_y[i] << '\t' << dipole::cells_field_array_z[i]  << "\t" << cells::pos_and_mom_array[4*i+0] << '\t' << cells_pos_and_mom_array[4*i+1] << '\t' << cells_pos_and_mom_array[4*i+2] << std::endl;
            // Multiply Hdemg by mu_0/4pi * 1e30 * mu_B to account for normalisation
            // of magnetisation and volume in angstrom
            dipole::cells_mu0Hd_field_array_x[i] = dipole::cells_mu0Hd_field_array_x[i] * constants::dipole_prefactor;
            dipole::cells_mu0Hd_field_array_y[i] = dipole::cells_mu0Hd_field_array_y[i] * constants::dipole_prefactor;
            dipole::cells_mu0Hd_field_array_z[i] = dipole::cells_mu0Hd_field_array_z[i] * constants::dipole_prefactor;
         }
      }

      // Each process has the fields of the cells holding its atoms; micromagnetic cells need all cells
      if(micromagnetic::discretisation_type > 0) dipole::gather_cell_fields();

      dipole::output_dipole_fields();

   } // end of dipole::internal::update_field() function
} // end of dipole namespace
