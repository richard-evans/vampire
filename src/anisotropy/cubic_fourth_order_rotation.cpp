//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Roberto Moreno Ortega and Richard Evans 2018. All rights reserved.
//
//   Email: roberto.moreno@york.ac.uk
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
      void cubic_fourth_order_rotation_fields(std::vector<double>& spin_array_x,
                                              std::vector<double>& spin_array_y,
                                              std::vector<double>& spin_array_z,
                                              std::vector<int>&    atom_material_array,
                                              std::vector<double>& field_array_x,
                                              std::vector<double>& field_array_y,
                                              std::vector<double>& field_array_z,
                                              const int start_index,
                                              const int end_index){

         // if not enabled then do nothing
         if(!internal::enable_cubic_fourth_order_rotation) return;

         // scale factor from derivative of E = -1/2 (sx^4 + sy^4 + sz^4)
         const double scale = 0.5*4.0;

         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; atom++){

            // get atom material
            const int mat = atom_material_array[atom];

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];

            // Unit vectors defining the easy axes for the cubic anisotropy in the new basis
            double e1[3] = { internal::mp[mat].kc_vector1[0],
                             internal::mp[mat].kc_vector1[1],
                             internal::mp[mat].kc_vector1[2] };

            // e2 *must* be orthogonal to e1
            double e2[3] = { internal::mp[mat].kc_vector2[0],
                             internal::mp[mat].kc_vector2[1],
                             internal::mp[mat].kc_vector2[2] };

            // es3 is orthogonal to e1 and e2
            double e3[3] = { internal::mp[mat].kc_vector3[0],
                             internal::mp[mat].kc_vector3[1],
                             internal::mp[mat].kc_vector3[2] };

            // get reduced anisotropy constant ku/mu_s
            const double kc4 = internal::kc4[mat];

            // calculate field (double negative from scale factor and negative derivative)
            const double k4 = scale*kc4;

            // calculate s.e for each vector
            const double sdote1 = sx * e1[0] + sy * e1[1] + sz * e1[2];
            const double sdote2 = sx * e2[0] + sy * e2[1] + sz * e2[2];
            const double sdote3 = sx * e3[0] + sy * e3[1] + sz * e3[2];

            // calculate (s.e)^3
            const double sdote1_3 = sdote1 * sdote1 * sdote1;
            const double sdote2_3 = sdote2 * sdote2 * sdote2;
            const double sdote3_3 = sdote3 * sdote3 * sdote3;

            // add field for cubic anisotropy in new basis
            field_array_x[atom] += k4*(sdote1_3*e1[0] + sdote2_3*e2[0] + sdote3_3*e3[0]);
            field_array_y[atom] += k4*(sdote1_3*e1[1] + sdote2_3*e2[1] + sdote3_3*e3[1]);
            field_array_z[atom] += k4*(sdote1_3*e1[2] + sdote2_3*e2[2] + sdote3_3*e3[2]);

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add fourth order cubic anisotropy in rotated basis (see manual)
      //---------------------------------------------------------------------------------
      double cubic_fourth_order_rotation_energy(const int atom,
                                                const int mat,
                                                const double sx,
                                                const double sy,
                                                const double sz){

         // Unit vectors defining the easy axes for the cubic anisotropy in the new basis
         double e1[3] = { internal::mp[mat].kc_vector1[0],
                          internal::mp[mat].kc_vector1[1],
                          internal::mp[mat].kc_vector1[2] };

         // e2 *must* be orthogonal to e1
         double e2[3] = { internal::mp[mat].kc_vector2[0],
                          internal::mp[mat].kc_vector2[1],
                          internal::mp[mat].kc_vector2[2] };

         // es3 is orthogonal to e1 and e2
         double e3[3] = { internal::mp[mat].kc_vector3[0],
                          internal::mp[mat].kc_vector3[1],
                          internal::mp[mat].kc_vector3[2] };

         // get reduced anisotropy constant ku/mu_s
         const double kc4 = internal::kc4[mat];

         // calculate s.e for each vector
         const double sdote1 = sx * e1[0] + sy * e1[1] + sz * e1[2];
         const double sdote2 = sx * e2[0] + sy * e2[1] + sz * e2[2];
         const double sdote3 = sx * e3[0] + sy * e3[1] + sz * e3[2];

         // calculate (s.e)^4
         const double sx4 = sdote1 * sdote1 * sdote1 * sdote1;
         const double sy4 = sdote2 * sdote2 * sdote2 * sdote2;
         const double sz4 = sdote3 * sdote3 * sdote3 * sdote3;

         return -0.5*kc4*(sx4 + sy4 + sz4);

      }

   } // end of internal namespace

} // end of anisotropy namespace
