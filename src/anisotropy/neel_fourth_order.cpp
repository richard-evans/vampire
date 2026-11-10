//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Lennon Sivaruban 2026. All rights reserved.
//
//   Email: lennon.sivaruban@manchester.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

//---------------------------------------------------------------------------------
// Fourth order Neel pair anisotropy for collinear spins, from Eq. 6 of
// P Nieves, J Tranchida, S Arapan and D Legut, Phys. Rev. B 103, 094437 (2021)
//
//    E_i = -1/2 sum_j { l [ (s.e)^2 - 1/3 ] + q [ (s.e)^4 - 6/7 (s.e)^2 + 3/35 ] }
//
// where e is the unit vector from atom i to neighbour j and l, q are the pair
// constants kij and qij. Grouping powers of (s.e) gives
//
//    E_i = - sum_j (l - 6/7 q)/2 (s.e)^2  -  sum_j q/2 (s.e)^4  +  sum_j ( l/6 - 3q/70 )
//
// The quadratic term is the standard Neel tensor with pair constant (l - 6/7 q)/2,
// calculated in initialize_neel.cpp and applied by neel_fields and neel_energy.
// This file adds the quartic and constant terms. The quartic term is a fully
// symmetric rank 4 tensor with 15 independent components, stored as the
// coefficients C_k of the monomials sx^a sy^b sz^c (a + b + c = 4):
//
//    - sum_j q/2 (s.e)^4 = - sum_k C_k sx^a sy^b sz^c
//
//    C_k = sum_j q/2 4!/(a! b! c!) ex^a ey^b ez^c
//
// so that the field (Tesla) H = -dE/ds = sum_k C_k d(sx^a sy^b sz^c)/ds.
//---------------------------------------------------------------------------------

namespace anisotropy{

   namespace internal{

      namespace{

         // exponents (a, b, c) of the monomials sx^a sy^b sz^c
         const int exponents[15][3] = { {4,0,0}, {0,4,0}, {0,0,4},
                                        {3,1,0}, {3,0,1}, {1,3,0}, {0,3,1}, {1,0,3}, {0,1,3},
                                        {2,2,0}, {2,0,2}, {0,2,2},
                                        {2,1,1}, {1,2,1}, {1,1,2} };

         // multinomial coefficients 4!/(a! b! c!) for each monomial
         const double multinomial[15] = { 1.0, 1.0, 1.0,
                                          4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
                                          6.0, 6.0, 6.0,
                                          12.0, 12.0, 12.0 };

      }

      //---------------------------------------------------------------------------------
      // Function to add the quartic and constant terms for one neighbour of atom, with
      // pair constants l and q (Tesla) and unit vector eij from atom to the neighbour
      //---------------------------------------------------------------------------------
      void add_neel_fourth_order_pair(const int atom, const double l, const double q, const double eij[3]){

         // powers e^n of the unit vector components (n = 0-4)
         const double ex[5] = {1.0, eij[0], eij[0]*eij[0], eij[0]*eij[0]*eij[0], eij[0]*eij[0]*eij[0]*eij[0]};
         const double ey[5] = {1.0, eij[1], eij[1]*eij[1], eij[1]*eij[1]*eij[1], eij[1]*eij[1]*eij[1]*eij[1]};
         const double ez[5] = {1.0, eij[2], eij[2]*eij[2], eij[2]*eij[2]*eij[2], eij[2]*eij[2]*eij[2]*eij[2]};

         for(int k = 0; k < 15; k++){
            internal::neel_fourth_order_tensor[15*atom + k] += 0.5 * q * multinomial[k] * ex[exponents[k][0]] * ey[exponents[k][1]] * ez[exponents[k][2]];
         }

         // constant terms -1/2 ( -l/3 + 3q/35 )
         internal::neel_fourth_order_energy_constant[atom] += l/6.0 - 3.0*q/70.0;

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add the quartic fourth order Neel anisotropy field
      //---------------------------------------------------------------------------------
      void neel_fourth_order_fields(std::vector<double>& spin_array_x,
                                    std::vector<double>& spin_array_y,
                                    std::vector<double>& spin_array_z,
                                    std::vector<int>&    atom_material_array,
                                    std::vector<double>& field_array_x,
                                    std::vector<double>& field_array_y,
                                    std::vector<double>& field_array_z,
                                    const int start_index,
                                    const int end_index){

         // if fourth order Neel anisotropy is not used, then do nothing
         if(!internal::enable_neel_fourth_order_anisotropy) return;

         // loop over all atoms
         for(int atom = start_index; atom < end_index; atom++){

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];

            // powers s^(n-1) of the spin components in element n, with element 0 = 0 so that
            // the derivative a s^(a-1) = a * sp[a] is zero for a = 0 and s^a = sp[a+1]
            const double sxp[6] = {0.0, 1.0, sx, sx*sx, sx*sx*sx, sx*sx*sx*sx};
            const double syp[6] = {0.0, 1.0, sy, sy*sy, sy*sy*sy, sy*sy*sy*sy};
            const double szp[6] = {0.0, 1.0, sz, sz*sz, sz*sz*sz, sz*sz*sz*sz};

            const int index = 15*atom; // get atom index in tensor array

            double hx = 0.0;
            double hy = 0.0;
            double hz = 0.0;

            // H = sum_k C_k d(sx^a sy^b sz^c)/ds
            for(int k = 0; k < 15; k++){
               const int a = exponents[k][0];
               const int b = exponents[k][1];
               const int c = exponents[k][2];
               const double ck = internal::neel_fourth_order_tensor[index + k];
               hx += ck * a * sxp[a] * syp[b+1] * szp[c+1];
               hy += ck * b * sxp[a+1] * syp[b] * szp[c+1];
               hz += ck * c * sxp[a+1] * syp[b+1] * szp[c];
            }

            // store net field
            field_array_x[atom] += hx;
            field_array_y[atom] += hy;
            field_array_z[atom] += hz;

         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add the quartic and constant fourth order Neel anisotropy energy
      //---------------------------------------------------------------------------------
      double neel_fourth_order_energy(const int atom,
                                      const int mat,
                                      const double sx,
                                      const double sy,
                                      const double sz){

         // powers s^n of the spin components (n = 0-4)
         const double sxp[5] = {1.0, sx, sx*sx, sx*sx*sx, sx*sx*sx*sx};
         const double syp[5] = {1.0, sy, sy*sy, sy*sy*sy, sy*sy*sy*sy};
         const double szp[5] = {1.0, sz, sz*sz, sz*sz*sz, sz*sz*sz*sz};

         const int index = 15*atom; // get atom index in tensor array

         // E = - sum_k C_k sx^a sy^b sz^c + constant
         double energy = 0.0;
         for(int k = 0; k < 15; k++){
            energy += internal::neel_fourth_order_tensor[index + k] * sxp[exponents[k][0]] * syp[exponents[k][1]] * szp[exponents[k][2]];
         }

         return -energy + internal::neel_fourth_order_energy_constant[atom];

      }

   } // end of internal namespace

} // end of anisotropy namespace
