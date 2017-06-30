//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: sw766@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <fstream>
#include <cmath>

// Vampire headers
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //---------------------------------------------------------------------------------
      // Function to add fourth order uniaxial anisotropy along vector e
      //
      //  Higher order anisotropies generally need to be described using spherical harmonics. The usual form (a
      //  series in S leads to cross pollution of terms, giving strange temperature dependencies.
      //
      //  The harmonics are described with Legendre polynomials with even order, which for 2nd, 4th and 6th are:
      //  ( http://en.wikipedia.org/wiki/Legendre_polynomials )
      //
      //  k2(sz) = - (1/2)  * (3sz^2 - 1)
      //  k4(sz) = - (1/8)  * (35sz^4 - 30sz^2 + 3)
      //  k6(sz) = - (1/16) * (231sz^6 - 315*sz^4 + 105sz^2 - 5)
      //
      //  The harmonics feature an arbritrary 2/3 factor compared with the usual form, and so in VAMPIRE these are
      //  renormalised to maintain consistency for the 2nd order terms.
      //
      //  The field induced by the harmonics is given by the first derivative w.r.t. sz. This can be projected onto
      //  any arbritrary direction ex,ey,ez allowing higher order anisotropy terms along any direction. This
      //  direction is shared with the other uniaxial anisotropy coefficients since they should not be used
      //  simultaneously.
      //
      //--------------------------------------------------------------------------------------------------------------
      //
      // Example 1: e == [0,0,1] || z
      //
      //       Energy (scalar) = -ku4 (1/8) * (35 Sz^4 - 30 Sz^2 + 3)
      //
      //       Tensor T[ij] = e[i] e[j] = [  0   0   0  ]
      //                                  [  0   0   0  ]
      //                                  [  0   0   1  ]
      //
      //       Energy (tensor) = S^2 . T4 . S^2 + S . T2 . S
      //
      //                       = [ Sx^2 Sy^2 Sz^2 ] [  0   0          0      ] [ Sx^2 ] + ...
      //                                            [  0   0          0      ] [ Sy^2 ]
      //                                            [  0   0   - 35 ku4 / 8  ] [ Sz^2 ]
      //
      //                       = Sx Txx Sx + Sx Txy Sy + Sx Txz Sz +
      //                         Sy Tyx Sx + Sy Tyy Sy + Sy Tyz Sz +
      //                         Sz Tzx Sx + Sz Tzy Sy + Sz Tzz Sz
      //
      //                       = Sz^2 Tzz Sz^2
      //
      // Need another example of off-diagonal anisotropy ... TBD
      //
      //---------------------------------------------------------------------------------
      void uniaxial_fourth_order(const unsigned int num_atoms, std::vector<int>& atom_material_array, std::vector<double>& inverse_mu_s){

         // constant factors
         const double fourth_order_prefactor = (  35.0 / 8.0) * (2.0 / 3.0 );
         const double second_order_prefactor = ( -30.0 / 8.0) * (2.0 / 3.0 );

         //----------------------------------------------------------------------------------
         // Loop over all atoms and calculate second order tensor components
         //----------------------------------------------------------------------------------
         for (int atom=0; atom < num_atoms; ++atom){

            // get atom material
            const unsigned int mat = atom_material_array[atom];

            const double i_mu_s = inverse_mu_s[mat];
            const double ku4 = internal::mp[mat].ku4;

            const double e[3] = { internal::mp[mat].ku_vector[0],
                                  internal::mp[mat].ku_vector[1],
                                  internal::mp[mat].ku_vector[2] };

            // Loop over tensor components and store anisotropy values in Tesla
            for (int i = 0; i < 3; ++i){
               for (int j = 0; j < 3; ++j){
                  internal::second_order_tensor[ index(atom, i, j) ] += second_order_prefactor * ku4 * e[i] * e[j] * i_mu_s;
                  internal::fourth_order_tensor[ index(atom, i, j) ] += fourth_order_prefactor * ku4 * e[i] * e[j] * i_mu_s;
               }
            }

         }

         return;

      }

   } // end of internal namespace

} // end of anisotropy namespace
