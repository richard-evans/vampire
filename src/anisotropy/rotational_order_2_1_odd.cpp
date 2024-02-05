//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Jack Collings 2024. All rights reserved.
//
//   Email: jbc525@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ Standard Library Headers

// Vampire Headers
#include "anisotropy.hpp"

// Anisotropy Module Headers
#include "internal.hpp"

namespace anisotropy{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{
      //---------------------------------------------------------------------------------
      // Function to add second order magnetocrystalline second order rotational
      // anisotropy based on vector e for the easy/hard/z axis and another vector for the
      // x axis.
      //
      // Higher order anisotropies generally need to be described using orthogonal
      // functions. The usual form, a series in S, leads to cross pollution of terms,
      // giving strange temperature dependencies.
      //
      // The anisotropies are described with a minimal orthogonal set expansion,
      // preserving the orthogonality of different orders while being simple to
      // implement and understand. Explicity the energies are described by normalising
      // the inner summation of the 2,4,6 order spherical harmonics to the prefactor
      // of the highest order term with an abritrary shift so that E(0) = 0.
      //
      // The rotational term here is given by
      // E_{2r-1} = -k_{2r-1}cos{theta}sin{theta}sin{phi}
      //
      // The field is found by taking the negative gradient w.r.t. the magnetic moment
      // basis.
      //
      //--------------------------------------------------------------------------------------------------------------

      void second_order_theta_first_order_phi_odd_fields( std::vector< double >& spin_array_x,
                                       std::vector< double >& spin_array_y,
                                       std::vector< double >& spin_array_z,
                                       std::vector< int >&    atom_material_array,
                                       std::vector< double >& field_array_x,
                                       std::vector< double >& field_array_y,
                                       std::vector< double >& field_array_z,
                                       const int start_index,
                                       const int end_index )
      {

         // if not enabled then do nothing
         if( !internal::enable_rotational_2_1_order_odd ) return;

         // Loop over all atoms between start and end index
         for( int atom = start_index; atom < end_index; ++atom )
         {

            // get atom material
            const int mat = atom_material_array[ atom ];

            // store spin direction in temporary variables
            const double sx = spin_array_x[ atom ];
            const double sy = spin_array_y[ atom ];
            const double sz = spin_array_z[ atom ];

            const double ex = internal::ku_vector[ mat ].x;
            const double ey = internal::ku_vector[ mat ].y;
            const double ez = internal::ku_vector[ mat ].z;

            const double gx = internal::kl_vector[ mat ].x;
            const double gy = internal::kl_vector[ mat ].y;
            const double gz = internal::kl_vector[ mat ].z;

            // calculate S_y
            const double Sy = sx * gx + sy * gy + sz * gz;

            // calculate S_z
            const double Sz = sx * ex + sy * ey + sz * ez;

            // get reduced anisotropy constant ku/mu_s
            const double k2r1_odd = internal::k2r1_odd[ mat ];

            // calculate field terms
            const double y_component = k2r1_odd * Sz;
            const double z_component = k2r1_odd * Sy;

            // sum in components of field
            field_array_x[atom] += y_component * gx + z_component * ex;
            field_array_y[atom] += y_component * gy + z_component * ey;
            field_array_z[atom] += y_component * gz + z_component * ez;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add 2-theta-1-phi odd anisotropy
      //---------------------------------------------------------------------------------

      double second_order_theta_first_order_phi_odd_energy( const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz )
      {

         const double ex = internal::ku_vector[ mat ].x;
         const double ey = internal::ku_vector[ mat ].y;
         const double ez = internal::ku_vector[ mat ].z;

         const double gx = internal::kl_vector[ mat ].x;
         const double gy = internal::kl_vector[ mat ].y;
         const double gz = internal::kl_vector[ mat ].z;

         // calculate cos{theta}sin{theta}sin{phi}
         const double Sz = sx * ex + sy * ey + sz * ez;
         const double Sy = sx * gx + sy * gy + sz * gz;

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double k2r1_odd = internal::k2r1_odd[ mat ];

         return - k2r1_odd * Sy * Sz;

      }

   }

}
