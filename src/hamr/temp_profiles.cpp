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

			// const double denominator = 1.0/ sqrt(8.0*log(2.0));
			// const double laser_sigma_x = hamr::internal::fwhm_x * denominator;
			// const double laser_sigma_y = hamr::internal::fwhm_y * denominator;
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
			// double sqrt_T = sqrt(temp);
			// // const double sqrt_T = sqrt( Tmin + DeltaT * exp(-cx2/(2.0 * laser_sigma_x2)) * exp(-cy2/(2.0 * laser_sigma_y2)) );

         // return sqrt_T;
			return temp;

      }

   } // end of namespace internal
} // end of namespace hamr
