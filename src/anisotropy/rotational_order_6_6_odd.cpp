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

namespace anisotropy
{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal
   {
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
      // E_{6r-6} = -k_{6r-6}sin^6{theta}sin{6phi}
      //
      // The field is found by taking the negative gradient w.r.t. the magnetic moment
      // basis.
      //
      //--------------------------------------------------------------------------------------------------------------

      // Define useful constants
      const double ten = 10.0;
      const double two = 2.0;
      const double three = 3.0;
      const double six = 6.0;
      const double five = 5.0;

      void sixth_order_theta_sixth_order_phi_odd_fields( std::vector< double >& spin_array_x,
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
         if( !internal::enable_rotational_6_6_order_odd ) return;

         // Loop over all atoms between start and end index
         for( int atom = start_index; atom < end_index; ++atom )
         {

            // get atom material
            const int mat = atom_material_array[ atom ];

            const double sx = spin_array_x[ atom ]; // store spin direction in temporary variables
            const double sy = spin_array_y[ atom ];
            const double sz = spin_array_z[ atom ];

            const double fx = internal::kr_vector[ mat ].x;
            const double fy = internal::kr_vector[ mat ].y;
            const double fz = internal::kr_vector[ mat ].z;

            const double gx = internal::kl_vector[ mat ].x;
            const double gy = internal::kl_vector[ mat ].y;
            const double gz = internal::kl_vector[ mat ].z;

            // calculate S_x and S_x^3 parts
            const double Sx = sx * fx + sy * fy + sz * fz;
            const double Sx2 = Sx * Sx;
            const double Sx4 = Sx2 * Sx2;

            // calculate S_y and S_y^3 parts
            const double Sy = sx * gx + sy * gy + sz * gz;
            const double Sy2 = Sy * Sy;
            const double Sy4 = Sy2 * Sy2;

            const double Sx2Sy2 = Sx2 * Sy2;

            // get reduced anisotropy constant ku/mu_s
            const double six_k6r6_odd = six * internal::k6r6_odd[ mat ];

            // calculate full form to add to field
            const double x_component = six_k6r6_odd * Sy * ( five * Sx4 - ten * Sx2Sy2 + Sy4 );
            const double y_component = six_k6r6_odd * Sx * ( Sx4 - ten * Sx2Sy2 + five * Sy4 );

            field_array_x[atom] += x_component * fx + y_component * gx;
            field_array_y[atom] += x_component * fy + y_component * gy;
            field_array_z[atom] += x_component * fz + y_component * gz;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add 6-theta-6-phi odd anisotropy
      //---------------------------------------------------------------------------------

      double sixth_order_theta_sixth_order_phi_odd_energy( const int atom,
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

         // Calculate sin^6{theta}sin{6phi}
         //          = sin^6{theta}( 6cos^5{phi}sin{phi} - 20cos^3{phi}sin^3{phi} + 6cos{phi}sin^5{phi} )
         //          = 6 * Sx^5 * Sy - 20 * Sx^3 * Sy^3 + 6 * Sx * Sy^5

         const double Sx = sx * fx + sy * fy + sz * fz;
         const double Sx2 = Sx * Sx;
         const double Sx4 = Sx2 * Sx2;

         const double Sy = sx * gx + sy * gy + sz * gz;
         const double Sy2 = Sy * Sy;
         const double Sy4 = Sy2 * Sy2;

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double two_k6r6_odd = two * internal::k6r6_odd[ mat ];

         return - two_k6r6_odd * Sx * Sy * ( three * Sx4 - ten * Sx2 * Sy2 + three * Sy4 );

      }

   }

}
