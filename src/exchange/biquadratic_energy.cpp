//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2018. All rights reserved.
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

   //---------------------------------------------------------------------------
   // Calculate isotropic biquadratic exchange energy for a single spin
   //
   //    E_bq = -J_bq (Si . Sj)^2 = -J_bq Sxi * Sxj ( Sxi * Sxj + Syi * Syj + Szi * Szj) ...
   //
   //---------------------------------------------------------------------------
   double spin_biquadratic_exchange_energy_isotropic(const int atom, const double sx, const double sy, const double sz){

   	// energy
   	double energy=0.0;

   	// Loop over neighbouring spins to calculate exchange
   	for(int nn = internal::biquadratic_neighbour_list_start_index[atom]; nn <= internal::biquadratic_neighbour_list_end_index[atom]; ++nn){

   		const int natom = internal::biquadratic_neighbour_list_array[nn];

         // get exchange constant between atoms
         const double Jbq = internal::bq_i_exchange_list[ internal::biquadratic_neighbour_interaction_type_array[nn] ].Jij;

         // load spin Sj components
         const double sjx = atoms::x_spin_array[natom];
         const double sjy = atoms::y_spin_array[natom];
         const double sjz = atoms::z_spin_array[natom];

         const double si_dot_sj = sx*sjx + sy*sjy + sz*sjz;

         // note: sum over j only (not sum over i for j) leads to a silent factor 1/2 in exchange energy value
         //       - must be normalised in statistics to account for double sum
   		energy -= Jbq * si_dot_sj * si_dot_sj;

   	}

   	return energy;

   }

   //-----------------------------------------------------------------------------------------
   // Function to calculate biquadratic exchange energy for spin atom
   //-----------------------------------------------------------------------------------------
   double single_spin_biquadratic_energy(const int atom, const double sx, const double sy, const double sz){

      // if biquadratic exchange not enabled do nothing
      if(!exchange::biquadratic) return 0.0;

      // select calculation based on exchange type
      switch(internal::biquadratic_exchange_type){

   		case exchange::isotropic:
            return spin_biquadratic_exchange_energy_isotropic(atom, sx, sy, sz);
            break;


         //case internal::vectorial:
            //return spin_biquadratic_exchange_energy_vectorial(atom, sx, sy, sz);
            //break;


         //case internal::tensorial:
            //return spin_biquadratic_exchange_energy_tensorial(atom, sx, sy, sz);
            //break;

         default:
            return 0.0;

   	}

      return 0.0;

   }

} // end of exchange namespace
