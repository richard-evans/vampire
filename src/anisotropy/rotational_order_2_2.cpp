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
      // E_{22} = -k_{22}sin^2{theta}cos{2phi}
      // 
      // The field is found by taking the negative gradient w.r.t. the magnetic moment 
      // basis and is detailed in an as yet unpublished paper.
      //
      //--------------------------------------------------------------------------------------------------------------
      void second_order_theta_second_order_phi_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index){

         // if not enabled then do nothing
         if(!internal::enable_rotational_2_2_order) return;

         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; atom++){

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

            // calculate S_x
            const double Sx = sx*fx + sy*fy + sz*fz;
            
            // calculate S_y
            const double Sy = sx*gx + sy*gy + sz*gz;

            // get reduced anisotropy constant ku/mu_s
            const double k2r2 = internal::k2r2[mat];

            // calculate field terms
            const double full_Sx = 2.0*k2r2*Sx;
            const double full_Sy = 2.0*k2r2*Sy;

            // sum x-component of field, where x-direction is represented by fx, fy, fz
            field_array_x[atom] += full_Sx*fx;
            field_array_y[atom] += full_Sx*fy;
            field_array_z[atom] += full_Sx*fz;

            // sum y-component of field, where y-direction is represented by gx, gy, gz
            field_array_x[atom] -= full_Sy*gx;
            field_array_y[atom] -= full_Sy*gy;
            field_array_z[atom] -= full_Sy*gz;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add 2-theta-2-phi anisotropy
      //---------------------------------------------------------------------------------

      double second_order_theta_second_order_phi_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz){

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double k2r2 = internal::k2r2[mat];

         const double fx = internal::kr_vector[mat].x;
         const double fy = internal::kr_vector[mat].y;
         const double fz = internal::kr_vector[mat].z;

         const double gx = internal::kl_vector[mat].x;
         const double gy = internal::kl_vector[mat].y;
         const double gz = internal::kl_vector[mat].z;

         // calculate sin^2{theta}cos{2phi} = sin^2{theta}cos^2{phi} - sin^2{theta}sin^2{phi} = S_x^2 - S_y^2
         const double Sx = sx*fx + sy*fy + sz*fz;
         const double Sy = sx*gx + sy*gy + sz*gz;
         const double sintheta2cos2phi = Sx*Sx - Sy*Sy;

         return -k2r2*sintheta2cos2phi;

      }
   }
}

