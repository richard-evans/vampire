//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "atoms.hpp" // for exchange list type defs
#include "exchange.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

namespace internal{

   //-----------------------------------------------------------------------------------------
   // Function to calculate biquadratic exchange fields for spins between start and end index
   //-----------------------------------------------------------------------------------------
   // Biquadratic exchange given by
   //
   //    E_bq = -J_bq (Si . Sj)^2 = -J_bq Sxi * Sxj ( Sxi * Sxj + Syi * Syj + Szi * Szj) ...
   //
   // Field given by:
   //
   //    H_bq^x = -dE/dSix = -2 Jbq Sxj ( Sxi * Sxj + Syi * Syj + Szi * Szj)
   //    H_bq^y = -dE/dSiy = -2 Jbq Syj ( Sxi * Sxj + Syi * Syj + Szi * Szj)
   //    H_bq^z = -dE/dSiz = -2 Jbq Szj ( Sxi * Sxj + Syi * Syj + Szi * Szj)
   //
   //-----------------------------------------------------------------------------------------
   void biquadratic_exchange_fields(const int start_index, // first atom for exchange interactions to be calculated
                                    const int end_index, // last +1 atom to be calculated
                                    const std::vector<int>& neighbour_list_start_index,
                                    const std::vector<int>& neighbour_list_end_index,
                                    const std::vector<int>& type_array, // type for atom
                                    const std::vector<int>& neighbour_list_array, // list of interactions between atoms
                                    const std::vector<int>& neighbour_interaction_type_array, // list of interaction type for each pair of atoms with value given in exchange list
                                    const std::vector<value_t>&  bq_i_exchange_list, // list of isotropic biquadratic exchange constants
                                    const std::vector<vector_t>& bq_v_exchange_list, // list of vectorial biquadratic exchange constants
                                    const std::vector<tensor_t>& bq_t_exchange_list, // list of tensorial biquadratic exchange constants
                                    const std::vector<double>& spin_array_x, // spin vectors for atoms
                                    const std::vector<double>& spin_array_y,
                                    const std::vector<double>& spin_array_z,
                                    std::vector<double>& field_array_x, // field vectors for atoms
                                    std::vector<double>& field_array_y,
                                    std::vector<double>& field_array_z){

      // loop over all atoms
		for(int atom = start_index; atom < end_index; ++atom){

         // temporary variables (registers) to calculate intermediate sum
			double hx = 0.0;
			double hy = 0.0;
			double hz = 0.0;

         // temporary constants for loop start and end indices
			const int start = neighbour_list_start_index[atom];
			const int end   = neighbour_list_end_index[atom]+1;

         // load spin Si components into temporary contants
         const double six = spin_array_x[atom];
         const double siy = spin_array_y[atom];
         const double siz = spin_array_z[atom];

         // loop over all neighbours
			for(int nn = start; nn < end; ++nn){

            // get neighbouring atom number
				const int natom = neighbour_list_array[nn];

            // get exchange constant between atoms
				const double twoJbq = 2.0*bq_i_exchange_list[ neighbour_interaction_type_array[nn] ].Jij;

            // load spin Sj components
            const double sjx = spin_array_x[natom];
            const double sjy = spin_array_y[natom];
            const double sjz = spin_array_z[natom];

            const double si_dot_sj = six*sjx + siy*sjy + siz*sjz;

				hx += twoJbq * sjx*si_dot_sj; // add exchange fields
				hy += twoJbq * sjy*si_dot_sj;
				hz += twoJbq * sjz*si_dot_sj;

			}

			field_array_x[atom] += hx; // save total field to field array
			field_array_y[atom] += hy;
			field_array_z[atom] += hz;

		}

		return;

	}

} // end of internal namespace

} // end of exchange namespace
