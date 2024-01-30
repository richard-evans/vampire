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
      // E_{2r-2} = -k_{2r-2}sin^2{theta}sin{2phi}
      //
      // The field is found by taking the negative gradient w.r.t. the magnetic moment
      // basis.
      //
      //--------------------------------------------------------------------------------------------------------------

      // Define useful constants
      const double two = 2.0;

      void second_order_theta_second_order_phi_odd_fields( std::vector< double >& spin_array_x,
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
         if( !internal::enable_rotational_2_2_order_odd ) return;

         // Loop over all atoms between start and end index
         for( int atom = start_index; atom < end_index; ++atom )
         {

            // get atom material
            const int mat = atom_material_array[ atom ];

            // store spin direction in temporary variables
            const double sx = spin_array_x[ atom ];
            const double sy = spin_array_y[ atom ];
            const double sz = spin_array_z[ atom ];

            const double fx = internal::kr_vector[ mat ].x;
            const double fy = internal::kr_vector[ mat ].y;
            const double fz = internal::kr_vector[ mat ].z;

            const double gx = internal::kl_vector[ mat ].x;
            const double gy = internal::kl_vector[ mat ].y;
            const double gz = internal::kl_vector[ mat ].z;

            // calculate S_x
            const double Sx = sx * fx + sy * fy + sz * fz;

            // calculate S_y
            const double Sy = sx * gx + sy * gy + sz * gz;

            // get reduced anisotropy constant ku/mu_s
            const double two_k2r2_odd = two * internal::k2r2_odd[mat];

            // calculate field terms
            const double x_component = two_k2r2_odd * Sy;
            const double y_component = two_k2r2_odd * Sx;

            // sum in components of field
            field_array_x[atom] += x_component * fx + y_component * gx;
            field_array_y[atom] += x_component * fy + y_component * gy;
            field_array_z[atom] += x_component * fz + y_component * gz;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add 2-theta-2-phi odd anisotropy
      //---------------------------------------------------------------------------------

      double second_order_theta_second_order_phi_odd_energy( const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz )
      {

         const double fx = internal::kr_vector[ mat ].x;
         const double fy = internal::kr_vector[ mat ].y;
         const double fz = internal::kr_vector[ mat ].z;

         const double gx = internal::kl_vector[ mat ].x;
         const double gy = internal::kl_vector[ mat ].y;
         const double gz = internal::kl_vector[ mat ].z;

         // calculate sin^2{theta}sin{2phi} = 2sin^2{theta}sin{phi}cos{phi} = 2 * S_x * S_y
         const double Sx = sx * fx + sy * fy + sz * fz;
         const double Sy = sx * gx + sy * gy + sz * gz;

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double two_k2r2_odd = two * internal::k2r2_odd[ mat ];

         return - two_k2r2_odd * Sx * Sy;

      }

   }

}
