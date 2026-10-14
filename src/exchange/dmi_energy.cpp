//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
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
   // Calculate Dzyaloshinskii-Moriya energy for a single spin
   //
   // The DMI is stored as the antisymmetric part of the tensorial exchange,
   // J^A = (J - J^T)/2, so that
   //
   //       E = - Si . J^A . Sj = Dij . (Si x Sj)
   //
   // This includes both the built-in DMI and any antisymmetric components of
   // user-defined tensorial interactions, but excludes symmetric anisotropic
   // terms such as Kitaev exchange. The DMI energy is a component of the
   // exchange energy and is not added separately to the total energy.
   //---------------------------------------------------------------------------
   double single_spin_dmi_energy(const int atom, const double sx, const double sy, const double sz){

      // DMI is only possible with tensorial exchange
      if(internal::exchange_type != exchange::tensorial) return 0.0;

      // energy
      double energy=0.0;

      // Loop over neighbouring spins to calculate dmi
      for(int nn = atoms::neighbour_list_start_index[atom]; nn <= atoms::neighbour_list_end_index[atom]; ++nn){

         const int natom = atoms::neighbour_list_array[nn];
         const zten_t& J = atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]];

         // antisymmetric components of exchange tensor
         const double Axy = 0.5 * (J.Jij[0][1] - J.Jij[1][0]);
         const double Axz = 0.5 * (J.Jij[0][2] - J.Jij[2][0]);
         const double Ayz = 0.5 * (J.Jij[1][2] - J.Jij[2][1]);

         const double S[3]={atoms::x_spin_array[natom],atoms::y_spin_array[natom],atoms::z_spin_array[natom]};

         // note: sum over j only (not sum over i for j) leads to a silent factor 1/2 in dmi energy value
         //       - must be normalised in statistics to account for double sum
         energy -= ( Axy * (sx * S[1] - sy * S[0]) +
                     Axz * (sx * S[2] - sz * S[0]) +
                     Ayz * (sy * S[2] - sz * S[1]) );

      }

      return energy;

   }

} // end of exchange namespace
