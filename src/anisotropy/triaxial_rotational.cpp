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
      void triaxial_second_order_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index){

         //if not enabled then do nothing
        if(!internal::enable_triaxial_anisotropy_rotated) return;

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

            double eA[3] = {internal::ku_triaxial_basis1x[mat],internal::ku_triaxial_basis1y[mat],internal::ku_triaxial_basis1z[mat]};
            double eB[3] = {internal::ku_triaxial_basis2x[mat],internal::ku_triaxial_basis2y[mat],internal::ku_triaxial_basis2z[mat]};
            double eC[3] = {internal::ku_triaxial_basis3x[mat],internal::ku_triaxial_basis3y[mat],internal::ku_triaxial_basis3z[mat]};

            // get reduced anisotropy constant ku/mu_s
            const double kA = internal::ku_triaxial_vector_x[mat];
            const double kB = internal::ku_triaxial_vector_y[mat];
            const double kC = internal::ku_triaxial_vector_z[mat];

            field_array_x[atom] += scale*(kA*eA[0]*sx + kB*eB[0]*sx + kC*eC[0]*sx);
            field_array_y[atom] += scale*(kA*eA[1]*sy + kB*eB[1]*sy + kC*eC[1]*sy);
            field_array_z[atom] += scale*(kA*eA[2]*sz + kB*eB[2]*sz + kC*eC[2]*sz);

      //      std::cout << "second" << scale*(kA*eA[0]*sx + kB*eB[0]*sx + kC*eC[0]*sx) << '\t' << scale*(kA*eA[1]*sy + kB*eB[1]*sy + kC*eC[1]*sy) << '\t' <<  scale*(kA*eA[2]*sz + kB*eB[2]*sz + kC*eC[2]*sz) <<std::endl;
        }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add second order uniaxial anisotropy in x,y and z
      // E = 2/3 * - ku2 (1/2)  * (3sz^2 - 1) == -ku2 sz^2 + const
      //---------------------------------------------------------------------------------
      double triaxial_second_order_energy(const int atom,
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

      void triaxial_fourth_order_fields(std::vector<double>& spin_array_x,
                                                    std::vector<double>& spin_array_y,
                                                    std::vector<double>& spin_array_z,
                                                    std::vector<int>&    atom_material_array,
                                                    std::vector<double>& field_array_x,
                                                    std::vector<double>& field_array_y,
                                                    std::vector<double>& field_array_z,
                                                    const int start_index,
                                                    const int end_index){

         // if not enabled then do nothing
         if(!internal::enable_triaxial_fourth_order_rotated) return;

         // constant factors
         const double oneo8 = 1.0/8.0;

         // rescaling prefactor
         const double scale = oneo8*2.0/3.0; // Factor to rescale anisotropies to usual scale


         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; atom++){

            // get atom material
            const int mat = atom_material_array[atom];

            double eA[3] = {internal::ku4_triaxial_basis1x[mat],internal::ku4_triaxial_basis1y[mat],internal::ku4_triaxial_basis1z[mat]};
            double eB[3] = {internal::ku4_triaxial_basis2x[mat],internal::ku4_triaxial_basis2y[mat],internal::ku4_triaxial_basis2z[mat]};
            double eC[3] = {internal::ku4_triaxial_basis3x[mat],internal::ku4_triaxial_basis3y[mat],internal::ku4_triaxial_basis3z[mat]};

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];


            // get reduced anisotropy constant ku/mu_s
            const double kA = internal::ku4_triaxial_vector_x[mat];
            const double kB = internal::ku4_triaxial_vector_y[mat];
            const double kC = internal::ku4_triaxial_vector_z[mat];

            const double sdoteA = eA[0]*sx + eA[1]*sy + eA[2]*sz;
            const double sdoteB = eB[0]*sx + eB[1]*sy + eB[2]*sz;
            const double sdoteC = eC[0]*sx + eC[1]*sy + eC[2]*sz;

            const double sdoteA3 = sdoteA*sdoteA*sdoteA;
            const double sdoteB3 = sdoteB*sdoteB*sdoteB;
            const double sdoteC3 = sdoteC*sdoteC*sdoteC;

            // calculate field (double negative from scale factor and negative derivative)
            //const double k4 = scale*(140.0*sdotk3 - 60.0*sdotk);
            const double k4A = scale*(sdoteA3 - (60.0/35.0)*sdoteA);
            const double k4B = scale*(sdoteB3 - (60.0/35.0)*sdoteB);
            const double k4C = scale*(sdoteC3 - (60.0/35.0)*sdoteC);

            field_array_x[atom] += kA*eA[0]*k4A +kB*eB[0]*k4B +kC*eC[0]*k4C;
            field_array_y[atom] += kA*eA[1]*k4A +kB*eB[1]*k4B +kC*eC[1]*k4C;
            field_array_z[atom] += kA*eA[2]*k4A +kB*eB[2]*k4B +kC*eC[2]*k4C;

          //  std::cout <<"fourth" <<  atom << '\t' << kA*eA[0]*k4A +kB*eB[0]*k4B +kC*eC[0]*k4C << '\t' <<  kA*eA[1]*k4A +kB*eB[1]*k4B +kC*eC[1]*k4C << '\t' << kA*eA[2]*k4A +kB*eB[2]*k4B +kC*eC[2]*k4C << '\t' << kA << '\t' << kB << '\t' << kC << "\t" << k4A << '\t' << k4B << '\t' << k4C << "\t" << std::endl;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add fourth order uniaxial anisotropy
      // E = 2/3 * - (1/8)  * (35sz^4 - 30sz^2 + 3)
      //---------------------------------------------------------------------------------
      double triaxial_fourth_order_energy(const int atom,
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
