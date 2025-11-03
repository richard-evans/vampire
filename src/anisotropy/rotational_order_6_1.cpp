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
      // E_{6r1} = -k_{6r1}sin{theta} ( cos^5{theta} - ( 10 / 11 ) cos^3(theta) + ( 5 / 33 ) cos(theta) ) cos{phi}
      // 
      // The field is found by taking the negative gradient w.r.t. the magnetic moment
      // basis.
      //
      //--------------------------------------------------------------------------------------------------------------
      
      // Define useful constants
		const double five = 5.0;
		const double oneothirtythree = 1.0 / 33.0;
		const double fiveothirtythree = 5.0 / 33.0;
		const double sixoeleven = 6.0 / 11.0;
		const double tenoeleven = 10.0 / 11.0;

      void sixth_order_theta_first_order_phi_fields(	std::vector< double >& spin_array_x,
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
         if( !internal::enable_rotational_6_1_order ) return;

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

            const double fx = internal::kr_vector[ mat ].x;
            const double fy = internal::kr_vector[ mat ].y;
            const double fz = internal::kr_vector[ mat ].z;

            // calculate S_z and S_z^3 parts
            const double Sz = sx * ex + sy * ey + sz * ez;
            const double Sz2 = Sz * Sz;
            const double Sz4 = Sz2 * Sz2;

            // calculate S_y and S_y^3 parts
            const double Sx = sx * fx + sy * fy + sz * fz;

            // get reduced anisotropy constant ku/mu_s
            const double k6r1 = internal::k6r1[ mat ];

            // calculate full form to add to field
            const double x_component = k6r1 * Sz * ( Sz4 - tenoeleven * Sz2 + fiveothirtythree );
            const double z_component = k6r1 * Sx * five * ( Sz4 - sixoeleven * Sz2 + oneothirtythree );
            
            field_array_x[ atom ] += z_component * ex + x_component * fx;
            field_array_y[ atom ] += z_component * ey + x_component * fy;
            field_array_z[ atom ] += z_component * ez + x_component * fz;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add 6-theta-1-phi anisotropy
      //---------------------------------------------------------------------------------

      double sixth_order_theta_first_order_phi_energy(	const int atom,
                                          					const int mat,
                                          					const double sx,
                                          					const double sy,
                                          					const double sz )
      {

         const double ex = internal::ku_vector[ mat ].x;
         const double ey = internal::ku_vector[ mat ].y;
         const double ez = internal::ku_vector[ mat ].z;

         const double fx = internal::kr_vector[ mat ].x;
         const double fy = internal::kr_vector[ mat ].y;
         const double fz = internal::kr_vector[ mat ].z;

         // calculate - k_{6r1} Sx Sz ( Sz^4 - ( 10 / 11 ) Sz^2 + 5 / 33 )
         const double Sz = sx * ex + sy * ey + sz * ez;
         const double Sz2 = Sz * Sz;

         const double Sx = sx * fx + sy * fy + sz * fz;

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double k6r1 = internal::k6r1[ mat ];

         return - k6r1 * Sx * Sz * ( Sz2 * Sz2 - tenoeleven * Sz2 + fiveothirtythree );

      }

   }

}
