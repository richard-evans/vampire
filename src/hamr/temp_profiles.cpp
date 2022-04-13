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
      // Function to calculate the atomic temperature due to a laser Gaussian profile
      //-----------------------------------------------------------------------------
      double calculate_gaussian_profile(const int atom,
                                       const double Tmin,
                                       const double DeltaT 
                                       ){

			const double laser_sigma_x2 = hamr::internal::laser_sigma_x * hamr::internal::laser_sigma_x;
			const double laser_sigma_y2 = hamr::internal::laser_sigma_y * hamr::internal::laser_sigma_y;
			const double px = hamr::internal::head_position_x;
			const double py = hamr::internal::head_position_y;
			const double cx = hamr::internal::atom_coords_x[atom];
			const double cy = hamr::internal::atom_coords_y[atom];
			const double cx2 = (cx-px)*(cx-px);
			const double cy2 = (cy-py)*(cy-py);

			const double denx = 2.0 * laser_sigma_x2;
			const double one_over_denx = 1.0/denx;
			const double deny = 2.0 * laser_sigma_y2;
			const double one_over_deny = 1.0/deny;
			double exp_x =  exp(-cx2 * one_over_denx);
			double exp_y =  exp(-cy2 * one_over_deny);

			double temp = Tmin + DeltaT * exp_x * exp_y;

			return temp;

      }


		/* ------------------------------------------------------------------  /
		/  Continuous HAMR process                                             /
		/  T(x,y,t) = (Tmax-Tmin)* exp( -(x-v*t)*(x-v*t)/(2*sigmax*sigmax) )   /
		/                        * exp( -(y-Yc)*(y-Yc)  /(2*sigmay*sigmay) )   /
		/                                                                      /
		/  v=Head speed; Yc=Head y-coordinate                                  /
		/  sigmax,sigmay = FWHM in x,y-directions                              /
		/ ------------------------------------------------------------------- */
      void apply_temperature_profile(const int start_index, const int end_index, const double Tmin, const double DeltaT)
		{
			// Define constants
			const double laser_sigma_x2 = hamr::internal::laser_sigma_x * hamr::internal::laser_sigma_x;
			const double laser_sigma_y2 = hamr::internal::laser_sigma_y * hamr::internal::laser_sigma_y;
			const double px = hamr::internal::head_position_x;
			const double py = hamr::internal::head_position_y;
			const double denx = 2.0 * laser_sigma_x2;
			const double one_over_denx = 1.0/denx;
			const double deny = 2.0 * laser_sigma_y2;
			const double one_over_deny = 1.0/deny;

			for(int atom=start_index;atom<end_index;atom++){

				const int imaterial=hamr::internal::atom_type_array[atom];
				double alpha = mp::material[imaterial].temperature_rescaling_alpha;
				double Tc = mp::material[imaterial].temperature_rescaling_Tc;

				const double cx = hamr::internal::atom_coords_x[atom];
				const double cy = hamr::internal::atom_coords_y[atom];
				const double cx2 = (cx-px)*(cx-px);
				const double cy2 = (cy-py)*(cy-py);

				// Get local temperature filed from application of heat profile
				const double exp_x =  exp(-cx2 * one_over_denx);
				const double exp_y =  exp(-cy2 * one_over_deny);
				const double temp = Tmin + DeltaT * exp_x * exp_y;

				// if T<Tc T/Tc = (T/Tc)^alpha else T = T
				double rescaled_temperature = temp < Tc ? Tc*pow(temp/Tc,alpha) : temp;
				const double sqrt_T=sqrt(rescaled_temperature);

				const double H_th_sigma = sqrt_T*mp::material[imaterial].H_th_sigma;

				hamr::internal::x_field_array[atom] *= H_th_sigma;
				hamr::internal::y_field_array[atom] *= H_th_sigma;
				hamr::internal::z_field_array[atom] *= H_th_sigma;
			}
			return;
      } // end of apply_temperature_profile

   } // end of namespace internal
} // end of namespace hamr
