//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains demag namespace and asssociated functions
///
/// @details Demag between cells uses the following format:
///          For each cell m = SUM(S.mu_s)
///
///         mu_o     (  3.(m.r_hat)r_hat - m   )    [    mu_o m   ])                4*pi e-7
/// H = ---------- . ( ------------------------) -  [ - --------  ]),   prefactor = ---------- = e+23
///         4*pi     (         |r|^3           )    [      3V     ])                4*pi e-30
///
///	An optional performaance optimisation can be made with the following matrix multiplication.
///	This is enabled by setting the flag demag::fast=true
///
///	H = prem * rij_matrix . m
///
///   rij_matrix = [ (3rxrx-1)/rij3 -1/3    3rxry                3rxrz                ]
///                [ 3rxry-1                (3ryry-1)/rij3 -1/3  3ryrz                ]
///                [ 3rxrz-1                3ryrz                (3rzrz-1)/rij3 -1/3) ]
///
///   The matrix is symmetric, and so only 6 numbers are needed:
///
///       rij_matrix[0] = xx
///       rij_matrix[1] = xy = yx
///       rij_matrix[2] = xz = zx
///
///       rij_matrix[3] = yy
///       rij_matrix[4] = yz = zy
///       rij_matrix[5] = zz
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    28/03/2011
/// @internal
///	Created:		18/09/2009
///	Revision:	  ---
///=====================================================================================
///
#include "atoms.hpp"
#include "cells.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "demag.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"


#include <cmath>
#include <iostream>
#include <time.h>

namespace demag{

	bool fast=false;

	int update_rate=100; /// timesteps between updates
	int update_time=-1; /// last update time

	const double prefactor=1.0e+23; // 1e-7/1e30

	std::vector <std::vector < double > > rij_xx;
	std::vector <std::vector < double > > rij_xy;
	std::vector <std::vector < double > > rij_xz;

	std::vector <std::vector < double > > rij_yy;
	std::vector <std::vector < double > > rij_yz;
	std::vector <std::vector < double > > rij_zz;

/// @brief Function to set r_ij matrix values
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    28/03/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		28/03/2011
///	Revision:	  ---
///=====================================================================================
///
void init(){

	std::cout << "Initialising demagnetisation field calculation" << std::endl;
	// check for calling of routine
	if(err::check==true){
		terminaltextcolor(RED);
		std::cerr << "demag::set_rij_matrix has been called " << vmpi::my_rank << std::endl;
		terminaltextcolor(WHITE);
	}
	if(demag::fast==true) {

      // timing function
      #ifdef MPICF
         double t1 = MPI_Wtime();
      #else
         time_t t1;
         t1 = time (NULL);
      #endif
		// Check memory requirements and print to screen
		zlog << zTs() << "Fast demagnetisation field calculation has been enabled and requires " << double(cells::num_cells)*double(cells::num_local_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;
		std::cout << "Fast demagnetisation field calculation has been enabled and requires " << double(cells::num_cells)*double(cells::num_local_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;

		// allocate arrays to store data [nloccell x ncells]
		for(int lc=0;lc<cells::num_local_cells; lc++){

			demag::rij_xx.push_back(std::vector<double>());
			demag::rij_xx[lc].resize(cells::num_cells,0.0);

			demag::rij_xy.push_back(std::vector<double>());
			demag::rij_xy[lc].resize(cells::num_cells,0.0);

			demag::rij_xz.push_back(std::vector<double>());
			demag::rij_xz[lc].resize(cells::num_cells,0.0);

			demag::rij_yy.push_back(std::vector<double>());
			demag::rij_yy[lc].resize(cells::num_cells,0.0);

			demag::rij_yz.push_back(std::vector<double>());
			demag::rij_yz[lc].resize(cells::num_cells,0.0);

			demag::rij_zz.push_back(std::vector<double>());
			demag::rij_zz[lc].resize(cells::num_cells,0.0);

		}

		// calculate matrix prefactors
		zlog << zTs() << "Precalculating rij matrix for demag calculation... " << std::endl;

		// loop over local cells
		for(int lc=0;lc<cells::num_local_cells;lc++){

			// reference global cell ID
			int i = cells::local_cell_array[lc];

			// Loop over all other cells to calculate contribution to local cell
			for(int j=0;j<cells::num_cells;j++){
				if(i!=j){

					const double rx = cells::x_coord_array[j]-cells::x_coord_array[i]; // Angstroms
					const double ry = cells::y_coord_array[j]-cells::y_coord_array[i];
					const double rz = cells::z_coord_array[j]-cells::z_coord_array[i];

					const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz);

					const double ex = rx*rij;
					const double ey = ry*rij;
					const double ez = rz*rij;

					const double rij3 = rij*rij*rij; // Angstroms

					rij_xx[lc][j] = demag::prefactor*((3.0*ex*ex - 1.0)*rij3);
					rij_xy[lc][j] = demag::prefactor*(3.0*ex*ey)*rij3;
					rij_xz[lc][j] = demag::prefactor*(3.0*ex*ez)*rij3;

					rij_yy[lc][j] = demag::prefactor*((3.0*ey*ey - 1.0)*rij3);
					rij_yz[lc][j] = demag::prefactor*(3.0*ey*ez)*rij3;
					rij_zz[lc][j] = demag::prefactor*((3.0*ez*ez - 1.0)*rij3);
		//			std::cout << lc << '\t' << rij_xx[lc][j] << '\t' << rij_xy[lc][j] << '\t' << rij_xz[lc][j] << '\t' << rij_yy[lc][j] << '\t' <<rij_yz[lc][j] << "\t" << rij_zz[lc][j] <<std::endl;
 				}
			}
		}

      #ifdef MPICF
         double t2 = MPI_Wtime();
      #else
         time_t t2;
         t2 = time (NULL);
      #endif
		zlog << zTs() << "Precalculation of rij matrix for demag calculation complete. Time taken: " << t2-t1 << "s."<< std::endl;

	}

	// timing function
   #ifdef MPICF
      double t1 = MPI_Wtime();
   #else
      time_t t1;
      t1 = time (NULL);
   #endif

	// now calculate fields
	demag::update();

	// timing function
   #ifdef MPICF
      double t2 = MPI_Wtime();
   #else
      time_t t2;
      t2 = time (NULL);
   #endif

   zlog << zTs() << "Time required for demag update: " << t2-t1 << "s." << std::endl;

}

/// @brief Function to recalculate demag fields using fast update method
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    28/03/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		28/03/2011
///	Revision:	  ---
///=====================================================================================
///
inline void fast_update(){

	// check for callin of routine
	if(err::check==true) {
		terminaltextcolor(RED);
		std::cerr << "demag::fast_update has been called " << vmpi::my_rank << std::endl;
		terminaltextcolor(WHITE);
	}
	// loop over local cells
	for(int lc=0;lc<cells::num_local_cells;lc++){

		int i = cells::local_cell_array[lc];

      // Calculate inverse volume from number of atoms in macrocell
      // V in A^3 == 1e-30 m3, mu_0 = 4pie-7 -> prefactor = pi*4e23/3V
      const double mu0_three_cell_volume = 8.0e23*M_PI/(3.0*cells::volume_array[i]);

      // Add self-demagnetisation
		cells::x_field_array[i]=mu0_three_cell_volume*cells::x_mag_array[i];
		cells::y_field_array[i]=mu0_three_cell_volume*cells::y_mag_array[i];
		cells::z_field_array[i]=mu0_three_cell_volume*cells::z_mag_array[i];

		// Loop over all other cells to calculate contribution to local cell
		for(int j=0;j<cells::num_cells;j++){

			const double mx = cells::x_mag_array[j];
			const double my = cells::y_mag_array[j];
			const double mz = cells::z_mag_array[j];

			cells::x_field_array[i]+=(mx*rij_xx[lc][j] + my*rij_xy[lc][j] + mz*rij_xz[lc][j]);
			cells::y_field_array[i]+=(mx*rij_xy[lc][j] + my*rij_yy[lc][j] + mz*rij_yz[lc][j]);
			cells::z_field_array[i]+=(mx*rij_xz[lc][j] + my*rij_yz[lc][j] + mz*rij_zz[lc][j]);

		}

		// Output data to vdp file
		//vdp << i << "\t" << cells::num_atoms_in_cell[i] << "\t";
		//vdp << cells::x_coord_array[i] << "\t" << cells::y_coord_array[i] << "\t" << cells::z_coord_array[i] << "\t";
		//vdp << cells::x_field_array[i] << "\t" << cells::y_field_array[i] << "\t" << cells::z_field_array[i] << "\t";
		//vdp << cells::x_mag_array[i] << "\t" << cells::y_mag_array[i] << "\t"<< cells::z_mag_array[i] << "\t" << inv_three_cell_volume*cells::z_mag_array[i] << std::endl;

	}
	//err::vexit();
}

/// @brief Function to recalculate demag fields using standard update method
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    28/03/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		28/03/2011
///	Revision:	  ---
///=====================================================================================
///
inline void std_update(){

	// check for callin of routine
	if(err::check==true){
		terminaltextcolor(RED);
		std::cerr << "demag::std_update has been called " << vmpi::my_rank << std::endl;
		terminaltextcolor(WHITE);
	}
	// loop over local cells
	for(int lc=0;lc<cells::num_local_cells;lc++){

		// get global cell ID
		int i = cells::local_cell_array[lc];

      // Calculate inverse volume from number of atoms in macrocell
      // V in A^3 == 1e-30 m3, mu_0 = 4pie-7 -> prefactor = pi*4e23/3V
      // But multiplied out at the end -> prefactor = pi*4/3V
      const double mu0_three_cell_volume = -4.0*M_PI/(3.0*cells::volume_array[i]);

      // Add self-demagnetisation
		cells::x_field_array[i]=mu0_three_cell_volume*cells::x_mag_array[i];
		cells::y_field_array[i]=mu0_three_cell_volume*cells::y_mag_array[i];
		cells::z_field_array[i]=mu0_three_cell_volume*cells::z_mag_array[i];

		// Loop over all other cells to calculate contribution to local cell
		for(int j=0;j<i;j++){

			const double mx = cells::x_mag_array[j];
			const double my = cells::y_mag_array[j];
			const double mz = cells::z_mag_array[j];

			const double dx = cells::x_coord_array[j]-cells::x_coord_array[i];
			const double dy = cells::y_coord_array[j]-cells::y_coord_array[i];
			const double dz = cells::z_coord_array[j]-cells::z_coord_array[i];

			const double drij = 1.0/sqrt(dx*dx+dy*dy+dz*dz);
			const double drij3 = drij*drij*drij; // Angstroms

			const double ex = dx*drij;
			const double ey = dy*drij;
			const double ez = dz*drij;

			const double s_dot_e = (mx * ex + my * ey + mz * ez);

			cells::x_field_array[i]+=(3.0 * s_dot_e * ex - mx)*drij3;
			cells::y_field_array[i]+=(3.0 * s_dot_e * ey - my)*drij3;
			cells::z_field_array[i]+=(3.0 * s_dot_e * ez - mz)*drij3;

		}

		for(int j=i+1;j<cells::num_cells;j++){

			const double mx = cells::x_mag_array[j];
			const double my = cells::y_mag_array[j];
			const double mz = cells::z_mag_array[j];

			const double dx = cells::x_coord_array[j]-cells::x_coord_array[i];
			const double dy = cells::y_coord_array[j]-cells::y_coord_array[i];
			const double dz = cells::z_coord_array[j]-cells::z_coord_array[i];

			const double drij = 1.0/sqrt(dx*dx+dy*dy+dz*dz);
			const double drij3 = drij*drij*drij;

			const double ex = dx*drij;
			const double ey = dy*drij;
			const double ez = dz*drij;

			const double s_dot_e = (mx * ex + my * ey + mz * ez);

			cells::x_field_array[i]+=(3.0 * s_dot_e * ex - mx)*drij3;
			cells::y_field_array[i]+=(3.0 * s_dot_e * ey - my)*drij3;
			cells::z_field_array[i]+=(3.0 * s_dot_e * ez - mz)*drij3;

		}

		cells::x_field_array[i]*=demag::prefactor;
		cells::y_field_array[i]*=demag::prefactor;
		cells::z_field_array[i]*=demag::prefactor;

		// Output data to vdp file
		//vdp << i << "\t" << cells::num_atoms_in_cell[i] << "\t";
		//vdp << cells::x_coord_array[i] << "\t" << cells::y_coord_array[i] << "\t" << cells::z_coord_array[i] << "\t";
		//vdp << cells::x_field_array[i] << "\t" << cells::y_field_array[i] << "\t" << cells::z_field_array[i] << "\t";
		//vdp << cells::x_mag_array[i] << "\t" << cells::y_mag_array[i] << "\t"<< cells::z_mag_array[i] << "\t" << inv_three_cell_volume*demag::prefactor*cells::z_mag_array[i] << std::endl;

	}
	//err::vexit();
}

/// @brief Wrapper Function to update demag fields
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    28/03/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		28/03/2011
///	Revision:	  ---
///=====================================================================================
///
void update(){

	if(err::check==true){
		terminaltextcolor(RED);
		std::cerr << "demag::update has been called " << vmpi::my_rank << std::endl;
		terminaltextcolor(WHITE);
	}
	// prevent double calculation for split integration (MPI)
	if(demag::update_time!=sim::time){

		// Check if update required
	  if(sim::time%demag::update_rate==0){

		//if updated record last time at update
		demag::update_time=sim::time;

		// update cell magnetisations
		cells::mag();

		// recalculate demag fields
		if(demag::fast==true) fast_update();
		else std_update();

		// For MPI version, only add local atoms
		#ifdef MPICF
			const int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
		#else
			const int num_local_atoms = atoms::num_atoms;
		#endif

		// Update Atomistic Dipolar Field Array
		for(int atom=0;atom<num_local_atoms;atom++){
			const int cell = atoms::cell_array[atom];

			// Copy field from macrocell to atomistic spin
			atoms::x_dipolar_field_array[atom]=cells::x_field_array[cell];
			atoms::y_dipolar_field_array[atom]=cells::y_field_array[cell];
			atoms::z_dipolar_field_array[atom]=cells::z_field_array[cell];
		}

		} // End of check for update rate
	} // end of check for update time

}

} // end of namespace demag
