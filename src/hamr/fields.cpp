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
					const double Hvecx,
					const double Hvecy,
					const double Hvecz,
					std::vector<double>& x_total_external_field_array,
					std::vector<double>& y_total_external_field_array,
					std::vector<double>& z_total_external_field_array
					){

		if(err::check==true){std::cout << "calculate_hamr_fields has been called" << std::endl;}

		// Define useful variables
		const double DeltaT = hamr::internal::Tmax - hamr::internal::Tmin;
		const double Hloc_parity_field=H_applied;

		// Add localised thermal field
		generate (hamr::internal::x_field_array.begin()+start_index,hamr::internal::x_field_array.begin()+end_index, mtrandom::gaussian);
		generate (hamr::internal::y_field_array.begin()+start_index,hamr::internal::y_field_array.begin()+end_index, mtrandom::gaussian);
		generate (hamr::internal::z_field_array.begin()+start_index,hamr::internal::z_field_array.begin()+end_index, mtrandom::gaussian);

		if(hamr::head_laser_on){

			// Apply local temperature field
			hamr::internal::apply_temperature_profile(start_index, end_index, hamr::internal::Tmin, DeltaT);

			// Add localised applied field
			hamr::internal::apply_field_spatial_box(start_index, end_index,
																Hloc_parity_field,
																Hvecx, Hvecy, Hvecz,
																x_total_external_field_array,
																y_total_external_field_array,
																z_total_external_field_array);
		}
		// Otherwise just use global temperature
		else{
			// unroll sigma for speed
			std::vector<double> sigma_prefactor(0);
			sigma_prefactor.reserve(mp::material.size());

			// Calculate material temperature (with optional rescaling)
			for(unsigned int mat=0;mat<mp::material.size();mat++){

			   // Calculate temperature rescaling
			   double alpha = mp::material[mat].temperature_rescaling_alpha;
			   double Tc = mp::material[mat].temperature_rescaling_Tc;
			   // if T<Tc T/Tc = (T/Tc)^alpha else T = T
			   double rescaled_temperature = temperature < Tc ? Tc*pow(temperature/Tc,alpha) : temperature;
			   double sqrt_T=sqrt(rescaled_temperature);
			   sigma_prefactor.push_back(sqrt_T*mp::material[mat].H_th_sigma);
			}

			for(int atom=start_index;atom<end_index;atom++){

				const int imaterial=hamr::internal::atom_type_array[atom];
				const double H_th_sigma = sigma_prefactor[imaterial];

				x_total_external_field_array[atom] = hamr::internal::x_field_array[atom] * H_th_sigma;
				y_total_external_field_array[atom] = hamr::internal::y_field_array[atom] * H_th_sigma;
				z_total_external_field_array[atom] = hamr::internal::z_field_array[atom] * H_th_sigma;
			}
		} // end of global temperature and field

		return;
	} // end of fields() functions

}
