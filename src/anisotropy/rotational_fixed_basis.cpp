//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard Evans 2020. All rights reserved.
//
//   Email: sarah.jenkins@york.ac.uk richard.evans@york.ac.uk
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
      // Function to add fourth order rotational anisotropy of the form
      //
      // E_4r = sin^3 theta cos (4 phi)
      //
      // In cartesian coordinates this expands to
      //
      // E_4r = 1 + Sz^4 - 8Sx^2 + 8Sx^2Sz^2 + 8Sx^4 - 2Sz^2
      //      = 1 + Sz^4 - 8Sy^2 + 8Sy^2Sz^2 + 8Sy^4 - 2Sz^2
      // E_4r = 1 - 8*Sx^2  + 8*Sx^4

      // The associated internal field (-dE/dS) is then
      //
      // Hx = 16 Sx (1 - Sz^2 - 2Sx^2)
      // Hy = 16 Sy (1 - Sz^2 - 2Sy^2)
      // Hz =  4 Sz (1 - Sz^2 - 4Sx^2)
      //--------------------------------------------------------------------------------------------------------------
      void rotational_fourth_order_fields_fixed_basis(std::vector<double>& spin_array_x,
                                                      std::vector<double>& spin_array_y,
                                                      std::vector<double>& spin_array_z,
                                                      std::vector<int>&    atom_material_array,
                                                      std::vector<double>& field_array_x,
                                                      std::vector<double>& field_array_y,
                                                      std::vector<double>& field_array_z,
                                                      const int start_index,
                                                      const int end_index){

         //if not enabled then do nothing
         if(!internal::enable_fourth_order_rotational) return;

         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; atom++){

            // get atom material
            const int mat = atom_material_array[atom];

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];

            // get reduced anisotropy constant ku/mu_s
            const double k4r = internal::k4r[mat];

            const double sx2 = sx*sx;
            const double sy2 = sy*sy;
            const double sz2 = sz*sz;

            field_array_x[atom] += k4r * 8.0 * sx * (1.0 - sz2 - 2.0 * sx2);
            field_array_y[atom] += k4r * 8.0 * sy * (1.0 - sz2 - 2.0 * sy2);
            field_array_z[atom] += k4r * 2.0 * sz * (1.0 - 2.0 * sz2 - 4.0 * sx2 - 4.0 * sy2);

        }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add second order uniaxial anisotropy in x,y and z
      // E = 2/3 * - ku2 (1/2)  * (3sz^2 - 1) == -ku2 sz^2 + const
      //---------------------------------------------------------------------------------
      double rotational_fourth_order_energy_fixed_basis(
         const int atom,
         const int mat,
         const double sx,
         const double sy,
         const double sz){

         // get reduced anisotropy constant ku/mu_s
         const double k4r = internal::k4r[mat];

         const double sx2 = sx*sx;
         const double sz2 = sz*sz;
         const double sx4 = sx2*sx2;
         const double sz4 = sz2*sz2;

         const double energy = k4r*(1.0 + sz4 - 8.0*sx2 + 8.0*sx2*sz2 + 8.0*sx4-2.0*sz2);

         return energy;

      }

   } // end of internal namespace

} // end of anisotropy namespace
