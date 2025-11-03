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
      // E_{44} = -k_{4r4}sin^4{theta}cos{4phi}
      //
      // The field is found by taking the negative gradient w.r.t. the magnetic moment
      // basis and is detailed in an as yet unpublished paper.
      //
      //--------------------------------------------------------------------------------------------------------------
      //Define useful constants
      const double four = 4.0;
      void fourth_order_theta_fourth_order_phi_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index){

         // if not enabled then do nothing
         if(!internal::enable_rotational_4_4_order) return;

         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; ++atom){

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];

            // get atom material
            const int mat = atom_material_array[atom];

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

            // get reduced anisotropy constant ku/mu_s
            const double four_k4r4 = four * internal::k4r4[mat];

            // calculate full form to add to field
            const double fullx = four_k4r4 * Sx * (Sx2 - 3 * Sy2);
            const double fully = four_k4r4 * Sy * (Sy2 - 3 * Sx2);

            field_array_x[atom] += fullx * fx;
            field_array_y[atom] += fullx * fy;
            field_array_z[atom] += fullx * fz;

            // sum y-component of field, where y-direction is represented by gx, gy, gz
            field_array_x[atom] += fully * gx;
            field_array_y[atom] += fully * gy;
            field_array_z[atom] += fully * gz;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add 6-theta-2-phi anisotropy
      //---------------------------------------------------------------------------------
      // Define useful constants
      const double six = 6.0;

      double fourth_order_theta_fourth_order_phi_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz){

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double k4r4 = internal::k4r4[mat];

         const double fx = internal::kr_vector[mat].x;
         const double fy = internal::kr_vector[mat].y;
         const double fz = internal::kr_vector[mat].z;

         const double gx = internal::kl_vector[mat].x;
         const double gy = internal::kl_vector[mat].y;
         const double gz = internal::kl_vector[mat].z;

         // calculate sin^4{theta}cos{4phi} = sin^4{theta} * ( 8 * cos^4{phi} - 8 * cos^2{phi} + 1 )
         //                                 = 8 * Sx^4 - 8 * sin^2{theta} * Sx^2 + sin^4{theta}
         //                                 = Sx^4 - 6 Sx^2 * Sy^2 + Sy^4
         const double Sx = sx * fx + sy * fy + sz * fz;
         const double Sx2 = Sx * Sx;

         const double Sy = sx * gx + sy * gy + sz * gz;
         const double Sy2 = Sy * Sy;

         return - k4r4 * (Sx2 * Sx2 - six * Sx2 * Sy2 + Sy2 * Sy2 );

      }
   }
}
