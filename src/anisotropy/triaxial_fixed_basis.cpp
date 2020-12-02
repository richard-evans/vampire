//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins, Andrew Naden and Richard Evans 2020. All rights reserved.
//
//   Email: sarah.jenkins@york.ac.uk ajn521@york.ac.uk richard.evans@york.ac.uk
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
      // Function to add second order uniaxial anisotropy along vectors x,y and z
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
      //  The harmonics feature an arbritrary +2/3 factor compared with the usual form, and so in VAMPIRE these are
      //  renormalised to maintain consistency for the 2nd order terms.
      //
      //  The field induced by the harmonics is given by the first derivative w.r.t. sz. This can be projected onto
      //  any arbritrary direction ex,ey,ez allowing higher order anisotropy terms along any direction. This
      //  direction is shared with the other uniaxial anisotropy coefficients since they should not be used
      //  simultaneously.
      //
      //--------------------------------------------------------------------------------------------------------------
      void triaxial_second_order_fields_fixed_basis(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index){

         //if not enabled then do nothing
        if(!internal::enable_triaxial_anisotropy) return;

         // rescaling prefactor
         // E = 2/3 * - ku2 (1/2)  * (3sz^2 - 1) == -ku2 sz^2 + const
         // H = -dE/dS = +2ku2 sz
          const double scale = 2.0; // 2*2/3 = 2 Factor to rescale anisotropies to usual scale
         //
         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; atom++){

            // get atom material
            const int mat = atom_material_array[atom];

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];

            // get reduced anisotropy constant ku/mu_s
            const double kx = internal::ku_triaxial_vector_x[mat];
            const double ky = internal::ku_triaxial_vector_y[mat];
            const double kz = internal::ku_triaxial_vector_z[mat];

            field_array_x[atom] += scale*kx*sx;
            field_array_y[atom] += scale*ky*sy;
            field_array_z[atom] += scale*kz*sz;
        }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add second order uniaxial anisotropy in x,y and z
      // E = 2/3 * - ku2 (1/2)  * (3sz^2 - 1) == -ku2 sz^2 + const
      //---------------------------------------------------------------------------------
       double triaxial_second_order_energy_fixed_basis(const int atom,

                                           const int mat,
                                           const double sx,
                                           const double sy,
                                           const double sz){
      
          // get reduced anisotropy constant ku/mu_s (Tesla)
      
          const double kx = internal::ku_triaxial_vector_x[mat];
          const double ky = internal::ku_triaxial_vector_y[mat];
          const double kz = internal::ku_triaxial_vector_z[mat];
          
	      const double energy = (sx*sx*kx + sy*sy*ky + sz*sz*kz);
      
          return -(energy);
      
       }

      void triaxial_fourth_order_fields_fixed_basis(std::vector<double>& spin_array_x,
                                                    std::vector<double>& spin_array_y,
                                                    std::vector<double>& spin_array_z,
                                                    std::vector<int>&    atom_material_array,
                                                    std::vector<double>& field_array_x,
                                                    std::vector<double>& field_array_y,
                                                    std::vector<double>& field_array_z,
                                                    const int start_index,
                                                    const int end_index){

         // if not enabled then do nothing
         if(!internal::enable_triaxial_fourth_order) return;

         // constant factors
         const double oneo8 = 1.0/8.0;

         // rescaling prefactor
         const double scale = oneo8*2.0/3.0; // Factor to rescale anisotropies to usual scale

         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; atom++){

            // get atom material
            const int mat = atom_material_array[atom];

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];


            // get reduced anisotropy constant ku/mu_s
            const double kx = internal::ku4_triaxial_vector_x[mat];
            const double ky = internal::ku4_triaxial_vector_y[mat];
            const double kz = internal::ku4_triaxial_vector_z[mat];

            const double sdotk  = (sx*kx + sy*ky + sz*kz);
            const double sdotk3 = sdotk*sdotk*sdotk;

            // calculate field (double negative from scale factor and negative derivative)
            //const double k4 = scale*(140.0*sdotk3 - 60.0*sdotk);
            const double k4 = scale*(sdotk3 - (60.0/35.0)*sdotk);

            field_array_x[atom] += kx*k4;
            field_array_y[atom] += ky*k4;
            field_array_z[atom] += kz*k4;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add fourth order uniaxial anisotropy
      // E = 2/3 * - (1/8)  * (35sz^4 - 30sz^2 + 3)
      //---------------------------------------------------------------------------------
      double triaxial_fourth_order_energy_fixed_basis(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz){

         // get reduced anisotropy constant ku/mu_s (Tesla)

         const double kx = internal::ku4_triaxial_vector_x[mat];
         const double ky = internal::ku4_triaxial_vector_y[mat];
         const double kz = internal::ku4_triaxial_vector_z[mat];

         const double sdotk  = (sx*kx + sy*ky + sz*kz);
         const double sdotk2 = sdotk*sdotk;

         // factor = 2/3 * -1/8 = -1/12 = -0.08333333333
         //return -0.08333333333*(35.0*sdotk2*sdotk2 - 30.0*sdotk2);
	     const double thirty_over_thirtyfive = 30.0/35.0;
         return -0.08333333333*(sdotk2*sdotk2 -  thirty_over_thirtyfive*sdotk2);

      }

   } // end of internal namespace

} // end of anisotropy namespace
