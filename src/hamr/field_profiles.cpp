//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo 2022. 
//
// All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>
// #include "math.h"

// Vampire headers
#include "hamr.hpp"
#include "material.hpp"
#include "vmpi.hpp"

// hamr headers
#include "internal.hpp"

namespace hamr{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Function to calculate the external field with trapezoidal temporal profile
      //-----------------------------------------------------------------------------
		double update_field_time_trapz_profile(const uint64_t current_time,  
                                          const uint64_t ramp_time,
                                          const uint64_t bit_time
                                          ){

		   const double Hmax = hamr::internal::Hmax; // max field
		   const double Hmin = hamr::internal::Hmin;
			double H_inc = 0.0;

			// Determine max magnitude of external field during initial ramp
			if(current_time <= ramp_time){
				H_inc = (Hmax-Hmin) / hamr::internal::H_ramp_time * mp::dt_SI;
			}
			// In the central region the max magnitude of the field is Hmax 
			// else if(current_time > ramp_time && current_time < bit_time-ramp_time){
			// 	H_inc = 0.0;
			// }
			// Decrease field in the final region
			else if(current_time >= bit_time-ramp_time && current_time <= bit_time){
				H_inc = -1.0*(Hmax-Hmin) / hamr::internal::H_ramp_time * mp::dt_SI;
			}

         return H_inc; 
      } // end of trapezoidal profile


      //-----------------------------------------------------------------------------
      // Function to apply field only to atoms within box around centre of writer coil
      //-----------------------------------------------------------------------------
      void apply_field_spatial_box(const int start_index,
					                  const int end_index,
											const double Hloc_parity_field,
											const double Hvecx,
											const double Hvecy,
											const double Hvecz,
					                  std::vector<double>& x_total_external_field_array,
					                  std::vector<double>& y_total_external_field_array,
					                  std::vector<double>& z_total_external_field_array
      									){
			// Add localised applied field
			for(int atom=start_index;atom<end_index;atom++){
				const double cx = hamr::internal::atom_coords_x[atom];
				const double cy = hamr::internal::atom_coords_y[atom];
				const double Hloc_min_x = hamr::internal::head_position_x - hamr::internal::H_bounds_x - hamr::internal::NPS;  // Shift field box in downtrack of NPS
				const double Hloc_max_x = hamr::internal::head_position_x + hamr::internal::H_bounds_x - hamr::internal::NPS;  // Shift field box in downtrack of NPS
				const double Hloc_min_y = hamr::internal::head_position_y - hamr::internal::H_bounds_y;
				const double Hloc_max_y = hamr::internal::head_position_y + hamr::internal::H_bounds_y;

				// If atoms within field box, add contribution from external field
				if((cx >= Hloc_min_x) && (cx <= Hloc_max_x) && (cy >= Hloc_min_y) && (cy <= Hloc_max_y)){
					hamr::internal::x_field_array[atom] += Hvecx*Hloc_parity_field;
					hamr::internal::y_field_array[atom] += Hvecy*Hloc_parity_field;
					hamr::internal::z_field_array[atom] += Hvecz*Hloc_parity_field;
				}
				x_total_external_field_array[atom] += hamr::internal::x_field_array[atom];
				y_total_external_field_array[atom] += hamr::internal::y_field_array[atom];
				z_total_external_field_array[atom] += hamr::internal::z_field_array[atom];
			}

			return;
      } // end of spatial box profile func



   } // end of namespace internal
} // end of namespace hamr
