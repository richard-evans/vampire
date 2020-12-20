//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2020. All rights reserved.
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
      // Function to add fourth order uniaxial anisotropy along vector e
      //
      //  Higher order anisotropies generally need to be described using orthogonal
      //  functions. The usual form (a series in S leads to cross pollution of terms,
      //  giving strange temperature dependencies.
      //
      //  The anisotropies are described with a minimal orthogonal set expansion,
      //  preserving the orthogonality of different orders while being simple to
      //  implement and understand. Explicity the energies are described by normalising
      //  the inner summation of the 2,4,6 order spherical harmonics to the prefactor
      //  of the highest order term with an abritrary shift so that E(0) = 0.
      //
      //  k2(sz) = k2 (1 - sz^2)
      //  k4(sz) = k4 (sz^4 - 30sz^2 / 35 - 5/35)
      //  k6(sz) = k6 (sz^6 - 315*sz^4/231 + 105sz^2/231 - 5/231)
      //
      //  The field induced by the harmonics is given by the first derivative w.r.t. sz.
      //  This can be projected onto any arbritrary direction ex,ey,ez allowing higher
      //  order anisotropy terms along any direction. This direction is shared with the
      //  other uniaxial anisotropy coefficients since they should not be used
      //  simultaneously.
      //
      //--------------------------------------------------------------------------------------------------------------
      void uniaxial_fourth_order_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index){

         // if not enabled then do nothing
         if(!internal::enable_uniaxial_fourth_order) return;

         // constant factors
         const double sixtyothirtyfive = 60.0/35.0;

         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; atom++){

            // get atom material
            const int mat = atom_material_array[atom];

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];

            const double ex = internal::ku_vector[mat].x;
            const double ey = internal::ku_vector[mat].y;
            const double ez = internal::ku_vector[mat].z;

            // get reduced anisotropy constant ku/mu_s
            const double ku4 = internal::ku4[mat];

            const double sdote  = (sx*ex + sy*ey + sz*ez);
            const double sdote3 = sdote*sdote*sdote;

            // calculate field (double negative from scale factor and negative derivative)
            const double k4 = -ku4*(4.0*sdote3 - sixtyothirtyfive*sdote);

            field_array_x[atom] += ex*k4;
            field_array_y[atom] += ey*k4;
            field_array_z[atom] += ez*k4;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add fourth order uniaxial anisotropy
      //---------------------------------------------------------------------------------
      // file scope constant
      const double fiveothirtyfive  = 5.0  / 35.0;
      const double thirtyothirtyfive = 30.0 / 35.0;

      double uniaxial_fourth_order_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz){

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double ku4 = internal::ku4[mat];


         const double ex = internal::ku_vector[mat].x;
         const double ey = internal::ku_vector[mat].y;
         const double ez = internal::ku_vector[mat].z;

         const double sdote  = (sx*ex + sy*ey + sz*ez);
         const double sdote2 = sdote*sdote;

         return ku4*(sdote2*sdote2 - thirtyothirtyfive*sdote2 - fiveothirtyfive);

      }

   } // end of internal namespace

} // end of anisotropy namespace
