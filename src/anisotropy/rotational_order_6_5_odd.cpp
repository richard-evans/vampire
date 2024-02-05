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
      // E_{6r-5} = -k_{6r-5}sin^5{theta}cos(theta)sin{5phi}
      // 
      // The field is found by taking the negative gradient w.r.t. the magnetic moment
      // basis.
      //
      //--------------------------------------------------------------------------------------------------------------
      
      // Define useful constants
		const double two = 2.0;
		const double three = 3.0;
      const double five = 5.0;
		const double six = 6.0;
		const double ten = 10.0;
		const double twelve = 12.0;
		const double sixteen = 16.0;
		const double twenty = 20.0;
		const double twentyfive = 25.0;
		const double sixty = 60.0;

      void sixth_order_theta_fifth_order_phi_odd_fields(	std::vector< double >& spin_array_x,
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
         if( !internal::enable_rotational_6_5_order_odd ) return;

         // Loop over all atoms between start and end index
         for( int atom = start_index; atom < end_index; ++atom )
         {

            // get atom material
            const int mat = atom_material_array[ atom ];

            const double sx = spin_array_x[ atom ]; // store spin direction in temporary variables
            const double sy = spin_array_y[ atom ];
            const double sz = spin_array_z[ atom ];

            const double ex = internal::ku_vector[ mat ].x;
            const double ey = internal::ku_vector[ mat ].y;
            const double ez = internal::ku_vector[ mat ].z;

            const double gx = internal::kl_vector[ mat ].x;
            const double gy = internal::kl_vector[ mat ].y;
            const double gz = internal::kl_vector[ mat ].z;

            // calculate S_z and S_z^3 parts
            const double Sz = sx * ex + sy * ey + sz * ez;
            const double Sz2 = Sz * Sz;
            const double Sz4 = Sz2 * Sz2;

            // calculate S_y and S_y^3 parts
            const double Sy = sx * gx + sy * gy + sz * gz;
            const double Sy2 = Sy * Sy;
            const double Sy4 = Sy2 * Sy2;

				const double Sy2Sz2 = Sy2 * Sz2;

            // get reduced anisotropy constant ku/mu_s
            const double k6r5_odd = internal::k6r5_odd[ mat ];

            // calculate full form to add to field
            const double y_component = k6r5_odd * Sz * five * ( sixteen * Sy4 + twelve * Sy2Sz2 + Sz4 - two * ( six * Sy2 + Sz2 ) + five );
            const double z_component = k6r5_odd * Sy * ( sixteen * Sy4 + sixty * Sy2Sz2 + twentyfive * Sz4 - ten * ( two * Sy2 + three * Sz2 ) + five );
            
            field_array_x[ atom ] += z_component * ex + y_component * gx;
            field_array_y[ atom ] += z_component * ey + y_component * gy;
            field_array_z[ atom ] += z_component * ez + y_component * gz;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add 6-theta-5-phi-odd anisotropy
      //---------------------------------------------------------------------------------

      double sixth_order_theta_fifth_order_phi_odd_energy(	const int atom,
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

         // calculate - k_{6r-5} (  )
         const double Sz = sx * ex + sy * ey + sz * ez;
         const double Sz2 = Sz * Sz;

         const double Sy = sx * gx + sy * gy + sz * gz;
         const double Sy2 = Sy * Sy;

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double k6r5_odd = internal::k6r5_odd[ mat ];

         return - k6r5_odd * Sy * Sz * ( sixteen * Sy2 * Sy2 + twenty * Sy2 * Sz2 + five * Sz2 * Sz2 - ten * ( two * Sy2 + Sz2 ) + five );

      }

   }

}
