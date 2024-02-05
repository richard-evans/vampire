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
      // E_{4r-1} = -k_{4r-1}sin{theta}sin{phi}(cos^3{theta} - (3/7)cos{theta})
      //
      // The field is found by taking the negative gradient w.r.t. the magnetic moment
      // basis.
      //
      //--------------------------------------------------------------------------------------------------------------
      
      //Define useful constants
      const double three = 3.0;
		const double threeoseven = 3.0 / 7.0;

      void fourth_order_theta_first_order_phi_odd_fields(	std::vector< double >& spin_array_x,
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
         if( !internal::enable_rotational_4_1_order_odd ) return;

         // Loop over all atoms between start and end index
         for( int atom = start_index; atom < end_index; ++atom )
         {

            // store spin direction in temporary variables
            const double sx = spin_array_x[ atom ];
            const double sy = spin_array_y[ atom ];
            const double sz = spin_array_z[ atom ];

            // get atom material
            const int mat = atom_material_array[ atom ];

				const double ex = internal::ku_vector[ mat ].x;
				const double ey = internal::ku_vector[ mat ].y;
				const double ez = internal::ku_vector[ mat ].z;

         	const double gx = internal::kl_vector[ mat ].x;
         	const double gy = internal::kl_vector[ mat ].y;
         	const double gz = internal::kl_vector[ mat ].z;

				// calculate S_z and S_z^2
				const double Sz = sx * ex + sy * ey + sz * ez;
         	const double Sz2 = Sz * Sz;

            // calculate S_y and S_y^3 parts
            const double Sy = sx * gx + sy * gy + sz * gz;
            const double Sy2 = Sy * Sy;

            // get reduced anisotropy constant ku/mu_s
            const double k4r1_odd = internal::k4r1_odd[ mat ];

            // calculate full form to add to field
            const double y_component = k4r1_odd * Sz * ( Sz2 - threeoseven );
				const double z_component = k4r1_odd *  Sy * ( three * Sz2 - threeoseven );

            // Sum components into field
            field_array_x[ atom ] += y_component * gx + z_component * ex;
            field_array_y[ atom ] += y_component * gy + z_component * ey;
            field_array_z[ atom ] += y_component * gz + z_component * ez;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add 4-theta-3-phi-odd anisotropy
      //---------------------------------------------------------------------------------
      double fourth_order_theta_first_order_phi_odd_energy(	const int atom,
                                          						const int mat,
                                          						const double sx,
                                          						const double sy,
                                          						const double sz )
      {

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double k4r1_odd = internal::k4r1_odd[ mat ];

			const double ex = internal::ku_vector[ mat ].x;
			const double ey = internal::ku_vector[ mat ].y;
			const double ez = internal::ku_vector[ mat ].z;

         const double gx = internal::kl_vector[ mat ].x;
         const double gy = internal::kl_vector[ mat ].y;
         const double gz = internal::kl_vector[ mat ].z;

         // Calculate - k{4r-1} sin{theta}sin{phi}(cos^3{theta} - (3/7)cos{theta})	= - k{4r-1} * Sy * ( Sz^3 - ( 3 / 7 ) * Sz )
         const double Sy = sx * gx + sy * gy + sz * gz;

			const double Sz = sx * ex + sy * ey + sz * ez;

         return - k4r1_odd * Sy * Sz * ( Sz * Sz - threeoseven );

      }

   }

}
