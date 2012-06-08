///
/// @file
/// @brief Contains demag namespace and asssociated functions
///
/// @details Demag between cells uses the following format:
///          For each cell m = SUM(S.mu_s)
///
///         mu_o     (  3.(m.r_hat)r_hat - m    [    m  ])                4*pi e-7
/// H = ---------- . ( ------------------------ [ - --- ]),   prefactor = ---------- = e+23
///       a^3        (         |r|^3            [    3  ])                4*pi e-30
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

	const double prefactor=1.0e+23;

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

	// check for calling of routine
	if(err::check==true) std::cerr << "demag::set_rij_matrix has been called " << vmpi::my_rank << std::endl;
	
	if(demag::fast==true) {
		
		// timing function
		time_t t1;
		t1 = time (NULL);
	
		// Check memory requirements and print to screen
		zlog << zTs() << "Fast demagnetisation field calculation has been enabled and requires " << double(cells::num_cells*cells::num_local_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;
		
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
				
					const double rx = cells::x_coord_array[j]-cells::x_coord_array[i];
					const double ry = cells::y_coord_array[j]-cells::y_coord_array[i];
					const double rz = cells::z_coord_array[j]-cells::z_coord_array[i];

					const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz);
					
					const double ex = rx*rij;
					const double ey = ry*rij;
					const double ez = rz*rij;

					const double rij3 = rij*rij*rij;

					rij_xx[lc][j] = demag::prefactor*((3.0*ex*ex - 1.0)*rij3);
					rij_xy[lc][j] = demag::prefactor*(3.0*ex*ey)*rij3;
					rij_xz[lc][j] = demag::prefactor*(3.0*ex*ez)*rij3;

					rij_yy[lc][j] = demag::prefactor*((3.0*ey*ey - 1.0)*rij3);
					rij_yz[lc][j] = demag::prefactor*(3.0*ey*ez)*rij3;
					rij_zz[lc][j] = demag::prefactor*((3.0*ez*ez - 1.0)*rij3);

				}
			}
		}
		
		time_t t2;
		t2 = time (NULL);
		zlog << zTs() << "Precalculation of rij matrix for demag calculation complete. Time taken: " << t2-t1 << "s."<< std::endl;
		
	}
	
	// timing function
	time_t t1;
	t1 = time (NULL);
		
	// now calculate fields
	demag::update();

	// timing function
	time_t t2;
	t2 = time (NULL);
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
	if(err::check==true) std::cerr << "demag::fast_update has been called " << vmpi::my_rank << std::endl;

	const double inv_three_cell_volume = -demag::prefactor*4.0*M_PI/(3.0*cells::size*cells::size*cells::size); // 1.0/(3*a*a*a)

	// loop over local cells
	for(int lc=0;lc<cells::num_local_cells;lc++){
		
		int i = cells::local_cell_array[lc];
		
		cells::x_field_array[i]=inv_three_cell_volume*cells::x_mag_array[i];
		cells::y_field_array[i]=inv_three_cell_volume*cells::y_mag_array[i];
		cells::z_field_array[i]=inv_three_cell_volume*cells::z_mag_array[i];

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
	if(err::check==true) std::cerr << "demag::std_update has been called " << vmpi::my_rank << std::endl;
	
	const double inv_three_cell_volume = -4.0*M_PI/(3.0*cells::size*cells::size*cells::size); // 1.0/(3*a*a*a)

	//std::cout << cells::num_local_cells << std::endl;

	// loop over local cells
	for(int lc=0;lc<cells::num_local_cells;lc++){
		
		// get global cell ID
		int i = cells::local_cell_array[lc];

		//std::cout << i << std::endl;
		// zero field arrays
		cells::x_field_array[i]=inv_three_cell_volume*cells::x_mag_array[i];
		cells::y_field_array[i]=inv_three_cell_volume*cells::y_mag_array[i];
		cells::z_field_array[i]=inv_three_cell_volume*cells::z_mag_array[i];
		
		// Loop over all other cells to calculate contribution to local cell 
		for(int j=0;j<i;j++){
				
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

	if(err::check==true) std::cerr << "demag::update has been called " << vmpi::my_rank << std::endl;

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

