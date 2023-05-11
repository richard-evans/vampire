//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Jack Collings and Richard Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ Standard Library Headers

// Vampire Headers
#include "anisotropy.hpp"

// Anisotropy Module Headers
#include "internal.hpp"

namespace anisotropy{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{
      //---------------------------------------------------------------------------------
      // Function to add second order magnetocrystalline second order rotational
      // anisotropy based on vector e for the easy/hard/z axis and another vector for the
      // x axis.
      //
      // Higher order anisotropies generally need to be described using orthogonal
      // functions. The usual form, a series in S, leads to cross pollution of terms,
      // giving strange temperature dependencies.
      //
      // The anisotropies are described with a minimal orthogonal set expansion,
      // preserving the orthogonality of different orders while being simple to
      // implement and understand. Explicity the energies are described by normalising
      // the inner summation of the 2,4,6 order spherical harmonics to the prefactor
      // of the highest order term with an abritrary shift so that E(0) = 0.
      //
      // The rotational term here is given by
      // E_{62} = -k_{6r2}sin^2{theta}(cos^4{theta} - (9/11)cos^2{theta} + (1/33))cos{2phi}
      //
      // The field is found by taking the negative gradient w.r.t. the magnetic moment 
      // basis and is detailed in an as yet unpublished paper.
      //
      //--------------------------------------------------------------------------------------------------------------

      // define useful consts
      const double thirtytwooeleven = 32.0 / 11.0;
      const double sixteenothirtythree = 16.0 / 33.0;
      const double two = 2.0;
      const double three = 3.0;

      void sixth_order_theta_second_order_phi_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index){

         // if not enabled then do nothing
         if(!internal::enable_rotational_6_2_order) return;

         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; ++atom){

            // get atom material
            const int mat = atom_material_array[atom];

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];

            const double fx = internal::kr_vector[mat].x;
            const double fy = internal::kr_vector[mat].y;
            const double fz = internal::kr_vector[mat].z;

            const double gx = internal::kl_vector[mat].x;
            const double gy = internal::kl_vector[mat].y;
            const double gz = internal::kl_vector[mat].z;

            // calculate S_x and S_x^3 parts
            const double Sx = sx * fx + sy * fy + sz * fz;
            const double Sx2 = Sx * Sx;

            // calculate S_y and S_y^3 parts
            const double Sy = sx * gx + sy * gy + sz * gz;
            const double Sy2 = Sy * Sy;

            const double Sx2pSy2 = Sx2 + Sy2;

            // get reduced anisotropy constant ku/mu_s
            const double two_k6r2 = two * internal::k6r2[mat];

            // calculate full form to add to field
            const double fullx = two_k6r2 * ( Sx * Sx2pSy2 * ( three * Sx2 - Sy2 ) - Sx * ( thirtytwooeleven * Sx2 - sixteenothirtythree ) );
            const double fully = two_k6r2 * ( Sy * Sx2pSy2 * ( Sx2 - three * Sy2 ) + Sy * ( thirtytwooeleven * Sy2 - sixteenothirtythree ) );

            field_array_x[atom] += fullx * fx + fully * gx;
            field_array_y[atom] += fullx * fy + fully * gy;
            field_array_z[atom] += fullx * fz + fully * gz;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add 6-theta-2-phi anisotropy
      //---------------------------------------------------------------------------------

      // Define useful constants
      const double sixteenoeleven = 16.0 / 11.0;
      const double one = 1.0;

      double sixth_order_theta_second_order_phi_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz){

         const double fx = internal::kr_vector[mat].x;
         const double fy = internal::kr_vector[mat].y;
         const double fz = internal::kr_vector[mat].z;

         const double gx = internal::kl_vector[mat].x;
         const double gy = internal::kl_vector[mat].y;
         const double gz = internal::kl_vector[mat].z;

         // calculate sin^2{theta} * ( cos^4{theta} - ( 6 / 11 ) * cos^2{theta} + ( 1 / 33 ) ) * cos{2phi}
         //          = ( sin^4{theta} - ( 16 / 11 ) * sin^2{theta} + 16 / 33 ) * ( 2 * sin^2{theta} cos^2{phi} - sin^2{theta} )
         const double Sx = sx * fx + sy * fy + sz * fz;
         const double Sx2 = Sx * Sx;

         const double Sy = sx * gx + sy * gy + sz * gz;
         const double Sy2 = Sy * Sy;

         const double sintheta2 = Sx2 + Sy2;

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double k6r2 = internal::k6r2[mat];

         return - k6r2 * (sintheta2 * sintheta2 - sixteenoeleven * sintheta2 + sixteenothirtythree ) * (Sx2 - Sy2);

      }
   }
}
