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
      // Function to add Neel pair anisotropy
      //
      // Example 1:
      //                               o -- ø -- o
      //                           y        |
      //                           ^        o
      //                           |__ > x
      //
      //
      //       Energy (scalar) = -ks/2 sum_j ( S_i . e_ij )^2            (1)
      //                       = -ks/2 ( S_i sum_j e_ij )^2
      //                       = -ks/2 ( S_i S_i \sum_j e_ij e_ij )      (2)
      //
      //       Vectors e_ij = [-1,0,0],[1,0,0],[0,-1,0]
      //
      //       Energy from Eq. 1 = -ks/2 [ ( -Sx . - Sx )^2 + ( Sx . Sx )^2 + ( -Sy . - Sy )^2 ]
      //                         = -ks/2 [ 2 Sx^2 + Sy^2 ]
      //
      //       Energy from Eq. 2 = -ks/2 [ Sx^2 ( -1.-1 + 1.1 ) + Sy^2 (-1.-1) ]
      //                         = -ks/2 [ 2 Sx^2 + Sy^2 ]
      //
      //       Tensor T[ij] = sum_j e_ij[i] e_ij[j] = [  2   0   0  ]
      //                                              [  0   1   0  ]
      //                                              [  0   0   0  ]
      //
      //       Energy (tensor) = S . T . S
      //
      //                       = [ Sx Sy Sz ] [  2   0   0  ] [ Sx ]
      //                                      [  0   1   0  ] [ Sy ]
      //                                      [  0   0   0  ] [ Sz ]
      //
      //                       = Sx Txx Sx + Sx Txy Sy + Sx Txz Sz +
      //                         Sy Tyx Sx + Sy Tyy Sy + Sy Tyz Sz +
      //                         Sz Tzx Sx + Sz Tzy Sy + Sz Tzz Sz
      //
      //                       = Sx 2 Sx + Sy Sy
      //
      //---------------------------------------------------------------------------------
      void neel_fields(std::vector<double>& spin_array_x,
                       std::vector<double>& spin_array_y,
                       std::vector<double>& spin_array_z,
                       std::vector<int>&    atom_material_array,
                       std::vector<double>& field_array_x,
                       std::vector<double>& field_array_y,
                       std::vector<double>& field_array_z,
                       const int start_index,
                       const int end_index){

         // if surface anisotropy is not used, then do nothing
         if(!internal::enable_neel_anisotropy) return;

         // loop over all atoms
         for(int atom = start_index; atom<end_index; atom++){

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];

            const int index = 9*atom; // get atom index in tensor array

            // Second order
            double hx = 2.0 * ( internal::neel_tensor[index + 0] * sx +
                                internal::neel_tensor[index + 1] * sy +
                                internal::neel_tensor[index + 2] * sz );

            double hy = 2.0 * ( internal::neel_tensor[index + 3] * sx +
                                internal::neel_tensor[index + 4] * sy +
                                internal::neel_tensor[index + 5] * sz );

            double hz = 2.0 * ( internal::neel_tensor[index + 6] * sx +
                                internal::neel_tensor[index + 7] * sy +
                                internal::neel_tensor[index + 8] * sz );

            //store net field
            field_array_x[atom] += hx;
            field_array_y[atom] += hy;
            field_array_z[atom] += hz;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add neel anisotropy energy
      //---------------------------------------------------------------------------------
      double neel_energy(const int atom,
                         const int mat,
                         const double sx,
                         const double sy,
                         const double sz){

         // get index for tensor
         const unsigned int index = 9*atom;

         double energy = 0.0;

         energy =  (sx * internal::neel_tensor[index + 0] * sx
                  + sx * internal::neel_tensor[index + 1] * sy
                  + sx * internal::neel_tensor[index + 2] * sz
                  + sy * internal::neel_tensor[index + 3] * sx
                  + sy * internal::neel_tensor[index + 4] * sy
                  + sy * internal::neel_tensor[index + 5] * sz
                  + sz * internal::neel_tensor[index + 6] * sx
                  + sz * internal::neel_tensor[index + 7] * sy
                  + sz * internal::neel_tensor[index + 8] * sz);

         // return energy after multiplying by -1
         return -1.0*energy;

      }

   } // end of internal namespace

} // end of anisotropy namespace
