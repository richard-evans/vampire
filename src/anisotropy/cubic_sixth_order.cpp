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
      // Function to add sixth order cubic anisotropy
      //---------------------------------------------------------------------------------
      void cubic_sixth_order_fields( std::vector<double>& spin_array_x,
                                     std::vector<double>& spin_array_y,
                                     std::vector<double>& spin_array_z,
                                     std::vector<int>&    atom_material_array,
                                     std::vector<double>& field_array_x,
                                     std::vector<double>& field_array_y,
                                     std::vector<double>& field_array_z,
                                     const int start_index,
                                     const int end_index){

         // if not enabled then do nothing
         if(!internal::enable_cubic_sixth_order) return;

         // scale factor from derivative of E = (sx^2 sy^2 sz^2)
         const double scale = -2.0;

         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; atom++){

            // get atom material
            const int mat = atom_material_array[atom];

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];

            // get reduced anisotropy constant ku/mu_s
            const double kc6 = internal::kc6[mat];

            // calculate intermediate expansions
            const double sx2 = sx*sx;
            const double sy2 = sy*sy;
            const double sz2 = sz*sz;

            // calculate field (double negative from scale factor and negative derivative)
            const double k6 = scale*kc6;

            field_array_x[atom] += sx*sy2*sz2*k6;
            field_array_y[atom] += sy*sz2*sx2*k6;
            field_array_z[atom] += sz*sx2*sy2*k6;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add fourth order cubic anisotropy
      // E = + k6 (sx^2 sy^2 sz^2)
      //---------------------------------------------------------------------------------
      double cubic_sixth_order_energy(const int atom,
                                      const int mat,
                                      const double sx,
                                      const double sy,
                                      const double sz){

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double kc6 = internal::kc6[mat];

         return kc6*sx*sx*sy*sy*sz*sz;

      }

   } // end of internal namespace

} // end of anisotropy namespace
