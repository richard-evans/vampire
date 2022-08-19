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
      // E_{42} = -k_{42}sin^2{theta}(cos^2{theta} - 1/7)cos{2phi}
      // 
      // The field is found by taking the negative gradient w.r.t. the magnetic moment 
      // basis and is detailed in an as yet unpublished paper.
      //
      //--------------------------------------------------------------------------------------------------------------
      
      // Define useful constants
      const double twelve_o_seven = 12.0/7.0;
      const double four = 4.0;
      
      void fourth_order_theta_second_order_phi_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index){

         // if not enabled then do nothing
         if(!internal::enable_rotational_4_2_order) return;

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
            const double Sxpart = twelve_o_seven * Sx;
            
            const double Sx3part = four * Sx * Sx * Sx;

            // calculate S_y and S_y^3 parts
            const double Sy = sx * gx + sy * gy + sz * gz;
            const double Sypart = twelve_o_seven * Sy;
            
            const double Sy3part = four * Sy * Sy * Sy;

            // get reduced anisotropy constant ku/mu_s
            const double k4r2 = internal::k4r2[mat];

            // calculate full form to add to field
            const double fullSx = k4r2 * (Sx3part - Sxpart);
            const double fullSy = k4r2 * (Sy3part - Sypart);
            
            field_array_x[atom] += - fullSx * fx + fullSy * gx;
            field_array_y[atom] += - fullSx * fy + fullSy * gy;
            field_array_z[atom] += - fullSx * fz + fullSy * gz;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add 4-theta-2-phi anisotropy
      //---------------------------------------------------------------------------------

      const double one_o_seven = 1.0 / 7.0;

      double fourth_order_theta_second_order_phi_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz){

         const double ex = internal::ku_vector[mat].x;
         const double ey = internal::ku_vector[mat].y;
         const double ez = internal::ku_vector[mat].z;

         const double fx = internal::kr_vector[mat].x;
         const double fy = internal::kr_vector[mat].y;
         const double fz = internal::kr_vector[mat].z;

         const double gx = internal::kl_vector[mat].x;
         const double gy = internal::kl_vector[mat].y;
         const double gz = internal::kl_vector[mat].z;

         // calculate cos^2{theta} = S_z
         const double costheta = sx * ex + sy * ey + sz * ez;
         const double costheta2 = costheta * costheta;

         // calculate sin^2{theta}cos{2phi} = sin^2{theta}cos^2{phi} - sin^2{theta}sin^2{phi}
         //          = S_x^2 - S_y^2
         const double Sx = sx * fx + sy * fy + sz * fz;
         const double Sy = sx * gx + sy * gy + sz * gz;
         const double sintheta2cos2phi = Sx * Sx - Sy * Sy;

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double k4r2 = internal::k4r2[mat];

         return - k4r2 * (costheta2 - one_o_seven) * sintheta2cos2phi;

      }
   }
}

