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
      // Function to add second order uniaxial anisotropy along vector e
      //
      // Example 1: e == [0,0,1] || z
      //
      //       Energy (scalar) = - ku Sz * Sz
      //
      //       Tensor T[ij] = e[i] e[j] = [  0   0   0  ]
      //                                  [  0   0   0  ]
      //                                  [  0   0   1  ]
      //
      //       Energy (tensor) = S . T . S
      //
      //                       = [ Sx Sy Sz ] [  0   0   0  ] [ Sx ]
      //                                      [  0   0   0  ] [ Sy ]
      //                                      [  0   0   1  ] [ Sz ]
      //
      //                       = Sx Txx Sx + Sx Txy Sy + Sx Txz Sz +
      //                         Sy Tyx Sx + Sy Tyy Sy + Sy Tyz Sz +
      //                         Sz Tzx Sx + Sz Tzy Sy + Sz Tzz Sz
      //
      //                       = Sz Tzz Sz
      //
      // Example 2: e == [0.707,0.707,0] || xy
      //
      //       Energy (scalar) = - ku (Sx ex + Sy ey + Sz ez)**2
      //                       = Sx ex Sx ex + 2 Sx Sy ex ey + Sy ey Sy ey
      //                       = 0.5 Sx Sx + Sx Sy + 0.5 Sy Sy
      //
      //       Tensor T[ij] = e[i] e[j] = [  0.5   0.5   0  ]
      //                                  [  0.5   0.5   0  ]
      //                                  [   0     0    0  ]
      //
      //       Energy (tensor) = S . T . S
      //
      //                       = [ Sx Sy Sz ] [  0.5   0.5   0  ] [ Sx ]
      //                                      [  0.5   0.5   0  ] [ Sy ]
      //                                      [   0     0    0  ] [ Sz ]
      //
      //                       = Sx Txx Sx + Sx Txy Sy + Sx Txz Sz +
      //                         Sy Tyx Sx + Sy Tyy Sy + Sy Tyz Sz +
      //                         Sz Tzx Sx + Sz Tzy Sy + Sz Tzz Sz
      //
      //                       = Sx Txx Sx + Sx Txy Sy +
      //                         Sy Tyx Sx + Sy Tyy Sy
      //
      //                       = 0.5 Sx Sx + 0.5 Sx Sy + 0.5 Sy Sx + 0.5 Sy Sy
      //
      //                       Q.E.D
      //
      //---------------------------------------------------------------------------------
      void uniaxial_second_order(const unsigned int num_atoms, std::vector<int>& atom_material_array, std::vector<double>& inverse_mu_s){

         //----------------------------------------------------------------------------------
         // Loop over all atoms and calculate second order tensor components
         //----------------------------------------------------------------------------------
         for (int atom=0; atom < num_atoms; ++atom){

            // get atom material
            const unsigned int mat = atom_material_array[atom];

            const double i_mu_s = inverse_mu_s[mat];
            const double ku2 = internal::mp[mat].ku2;

            const double e[3] = { internal::mp[mat].ku_vector[0],
                                  internal::mp[mat].ku_vector[1],
                                  internal::mp[mat].ku_vector[2] };

            // Loop over tensor components and store anisotropy values in Tesla (-dE/dS)
            for (int i = 0; i < 3; ++i){
               for (int j = 0; j < 3; ++j){
                  internal::second_order_tensor[ index(atom, i, j) ] += ku2 * e[i] * e[j] * i_mu_s;
               }
            }

         }

         return;

      }

   } // end of internal namespace

} // end of anisotropy namespace
