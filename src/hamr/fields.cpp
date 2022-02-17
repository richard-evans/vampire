//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo 2022. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>
#include <algorithm>

// Vampire headers
#include "atoms.hpp"
#include "errors.hpp"
#include "hamr.hpp"
#include "material.hpp"
#include "random.hpp"
#include "vio.hpp"

// hamr headers
#include "internal.hpp"

namespace hamr{

	void fields(const int start_index,
					const int end_index,
					double H_applied,
					const double temperature,
					const double Tmin,
					const double Tmax,
					const double Hvecx,
					const double Hvecy,
					const double Hvecz,
					std::vector<double>& x_total_external_field_array,
					std::vector<double>& y_total_external_field_array,
					std::vector<double>& z_total_external_field_array
					){
		///======================================================
		/// 		Function to calculate HAMR fields
		///		Richard Evans 2011
		///		Revision: Andrea Meo 2022
		///======================================================

		if(err::check==true){std::cout << "calculate_hamr_fields has been called" << std::endl;}

		// Define useful variables
		const double DeltaT = Tmax - Tmin;

		// declare head-field variables hamr::internal::H_bounds_min[0]
		const double H_osc_amplit=hamr::internal::H_osc_amplit; 
		const double Hloc_parity_field=H_applied*(-1.0)*double(2*(int(hamr::internal::head_position[0]/H_osc_amplit)%2)-1);
		//const double Hloc_parity_field=H_applied;

		// Add localised thermal field
		generate (x_total_external_field_array.begin()+start_index,x_total_external_field_array.begin()+end_index, mtrandom::gaussian);
		generate (y_total_external_field_array.begin()+start_index,y_total_external_field_array.begin()+end_index, mtrandom::gaussian);
		generate (z_total_external_field_array.begin()+start_index,z_total_external_field_array.begin()+end_index, mtrandom::gaussian);

		if(hamr::head_laser_on){
			for(int atom=start_index;atom<end_index;atom++){
				const int imaterial=hamr::internal::atom_type_array[atom];
				const double sqrt_T = hamr::internal::calculate_gaussian_profile(atom, Tmin, DeltaT);
				const double H_th_sigma = sqrt_T*mp::material[imaterial].H_th_sigma;
				x_total_external_field_array[atom] *= H_th_sigma; 
				y_total_external_field_array[atom] *= H_th_sigma; 
				z_total_external_field_array[atom] *= H_th_sigma; 
			}

			// Add localised applied field
			hamr::internal::apply_field_spatial_box(start_index, end_index, 
																Hloc_parity_field,
																Hvecx, Hvecy, Hvecz, 
																x_total_external_field_array, 
																y_total_external_field_array, 
																z_total_external_field_array);
		}
		else{
			// Otherwise just use global temperature
			double sqrt_T=sqrt(temperature);
			for(int atom=start_index;atom<end_index;atom++){
				const int imaterial=hamr::internal::atom_type_array[atom];
				const double H_th_sigma = sqrt_T*mp::material[imaterial].H_th_sigma;
				x_total_external_field_array[atom] *= H_th_sigma; 
				y_total_external_field_array[atom] *= H_th_sigma; 
				z_total_external_field_array[atom] *= H_th_sigma; 
			}
		} // end of global temperature and field

		return;
	} // end of fields() functions

}