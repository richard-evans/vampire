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
      void cubic_fourth_order_fields(std::vector<double>& spin_array_x,
                                     std::vector<double>& spin_array_y,
                                     std::vector<double>& spin_array_z,
                                     std::vector<int>&    atom_material_array,
                                     std::vector<double>& field_array_x,
                                     std::vector<double>& field_array_y,
                                     std::vector<double>& field_array_z,
                                     const int start_index,
                                     const int end_index){

         // if not enabled then do nothing
         if(!internal::enable_cubic_fourth_order) return;

         // scale factor from derivative of E = -1/2 (sx^4 + sy^4 + sz^4)
         const double scale = 0.5*4.0;

         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; atom++){

            // get atom material
            const int mat = atom_material_array[atom];

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];

            // get reduced anisotropy constant ku/mu_s
            const double kc4 = internal::kc4[mat];

            // calculate field (double negative from scale factor and negative derivative)
            const double k4 = scale*kc4;

            field_array_x[atom] += sx*sx*sx*k4;
            field_array_y[atom] += sy*sy*sy*k4;
            field_array_z[atom] += sz*sz*sz*k4;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add fourth order cubic anisotropy
      // E = -1/2 k4 (sx^4 + sy^4 + sz^4)
      //---------------------------------------------------------------------------------
      double cubic_fourth_order_energy(const int atom,
                                       const int mat,
                                       const double sx,
                                       const double sy,
                                       const double sz){

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double kc4 = internal::kc4[mat];

         const double sx4 = sx*sx*sx*sx;
         const double sy4 = sy*sy*sy*sy;
         const double sz4 = sz*sz*sz*sz;

         return -0.5*kc4*(sx4 + sy4 + sz4);

      }

   } // end of internal namespace

} // end of anisotropy namespace
