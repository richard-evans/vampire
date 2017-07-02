//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
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
      // Function to add fourth order cubic anisotropy
      //---------------------------------------------------------------------------------
      void cubic_fourth_order(const unsigned int num_atoms, std::vector<int>& atom_material_array, std::vector<double>& inverse_mu_s){

         //----------------------------------------------------------------------------------
         // Loop over all atoms and calculate second order tensor components
         //
         // conventional cubic anisotropy given by
         // E = + (mx^2 my^2 + my^2 mz^2 + mx_2 mz^2)
         //   = -1/2 ( 1 + mx^4 + my^4 + mz^4 )
         //
         // factor -1/2 and irrelevant constant into calculation
         //
         //----------------------------------------------------------------------------------
         for (int atom=0; atom < num_atoms; ++atom){

            // get atom material
            const unsigned int mat = atom_material_array[atom];

            const double i_mu_s = inverse_mu_s[mat];

            // Strore constant including prefactor and conversion to Tesla (-dE/dS)
            const double kc4 = 0.5 * internal::mp[mat].kc4;

            // Loop over tensor components and store anisotropy values in Tesla
            for (int i = 0; i < 3; ++i){
               for (int j = 0; j < 3; ++j){
                  if( i == j ){
                     internal::fourth_order_tensor[ index(atom, i, j) ] += kc4 * i_mu_s;
                  }
               }
            }

         }

         std::ofstream ofile("energy.txt");

         for(int i = 0; i <= 180; i+=2){
            for(int j = 0; j <= 360; j+=2){
               double sx = sin(i*3.14159/180.0)*cos(j*3.14159/180.0);
               double sy = sin(i*3.14159/180.0)*sin(j*3.14159/180.0);
               double sz = cos(i*3.14159/180.0);

               double energy = /*- sx * internal::second_order_tensor[0] * sx
                               - sx * internal::second_order_tensor[1] * sy
                               - sx * internal::second_order_tensor[2] * sz
                               - sy * internal::second_order_tensor[3] * sx
                               - sy * internal::second_order_tensor[4] * sy
                               - sy * internal::second_order_tensor[5] * sz
                               - sz * internal::second_order_tensor[6] * sx
                               - sz * internal::second_order_tensor[7] * sy
                               - sz * internal::second_order_tensor[8] * sz*/

                               - sx * sx * internal::fourth_order_tensor[0] * sx * sx
                               - sx * sx * internal::fourth_order_tensor[1] * sy * sy
                               - sx * sx * internal::fourth_order_tensor[2] * sz * sz
                               - sy * sy * internal::fourth_order_tensor[3] * sx * sx
                               - sy * sy * internal::fourth_order_tensor[4] * sy * sy
                               - sy * sy * internal::fourth_order_tensor[5] * sz * sz
                               - sz * sz * internal::fourth_order_tensor[6] * sx * sx
                               - sz * sz * internal::fourth_order_tensor[7] * sy * sy
                               - sz * sz * internal::fourth_order_tensor[8] * sz * sz;

               ofile << i << "\t" << j << "\t" << energy << std::endl;
            }
            ofile << std::endl;
         }

         return;

      }

   } // end of internal namespace

} // end of anisotropy namespace
