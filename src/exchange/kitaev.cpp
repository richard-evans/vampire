//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and S Jenkins 2021. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "exchange.hpp"
#include "material.hpp"
#include "vio.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

namespace internal{

   //------------------------------------------------------------------------------
   // Function to calculate Kitaev exchange interactions. The Kitaev interactions
   // are directional exchange (2-ion) anisotropy along the bond direction e_ij.
   //
   //                                O --------- O
   //                                i    e_ij   j
   //
   // For neighbours i and j, the Kitaev interaction is reflected as an
   // anisotropy with an easy axis along e_ij. This is done by computing the
   // direction vector between neighbour atoms and then adding an appropriate
   // contribution to the exchange tensor. A standard cutoff range is applied to
   // limit the interaction range to nearest neighbours.
   //
   //------------------------------------------------------------------------------
   void calculate_kitaev(std::vector<std::vector <neighbours::neighbour_t> >& cneighbourlist){

      // if kitaev is not needed then do nothing
      if(!internal::enable_kitaev) return;

      // Print informative message to log file
      zlog << zTs() << "Calculating Kitaev interactions" << std::endl;

      // temporary tensor for calculating sum
      std::vector<double> tmp_tensor(9,0.0);

      // get cutoff range
      const double kcor = internal::kitaev_cutoff_range;
      const double cutoff_sq = kcor * kcor;

      unsigned int counter = 0;

      //std::ofstream ofile("kitaev.txt");

      //	Loop over all atoms i
      for(uint64_t i=0; i < static_cast<uint64_t>(atoms::num_atoms); i++){

         // get material id for atom
         const unsigned int imat = atoms::type_array[i];

         // get inverse moment (for normalisation)
         const double i_mu_s = 1.0/mp::material[imat].mu_s_SI;

         // loop over all neighbours j
         for(unsigned int j = 0; j < cneighbourlist[i].size(); j++){

            // get atom number for neighbour i
            const unsigned int nj = cneighbourlist[i][j].nn;

            // get material id for j atom
            const unsigned int jmat = atoms::type_array[nj];

            // ignore self interaction
            if(i != nj){

               // get atomic position vector i->j
               double eij[3]={cneighbourlist[i][j].vx, cneighbourlist[i][j].vy, cneighbourlist[i][j].vz};

               // get length of vector
               const double mod_eij_sq = eij[0]*eij[0] + eij[1]*eij[1] + eij[2]*eij[2];

               //ofile << i << "\t" << nj << std::endl <<
               //         "\t\teij: " << eij[0]/sqrt(mod_eij_sq) << "\t" << eij[1]/sqrt(mod_eij_sq) << "\t" << eij[2]/sqrt(mod_eij_sq) << std::endl;

               // check if that the bond length is less than cutoff range
               if( mod_eij_sq <= cutoff_sq){

                  // normalise to unit vector
                  const double inv_rij=1.0/sqrt(mod_eij_sq);

                  // reinitialise exchange tensor for each i-j interaction
                  for(int idx = 0; idx < 9; idx++) tmp_tensor[idx] = 0.0;

                  // normalise components to unit vector
                  for(int idx = 0; idx < 3; idx++){
                     eij[idx] = eij[idx] * inv_rij;
                  }

                  // Set pair Kitaev constant normalised to mu_s_i
                  const double Kij = exchange::internal::mp[imat].kitaev[jmat] * i_mu_s;

                  //---------------------------------------------------------------
                  //        Expression of Kitaev within the exchange tensor
                  //---------------------------------------------------------------
                  //
                  //       E = Kij (Si . e_ij)(Sj . e_ij)
                  //
                  //       Si . e_ij = Si[0]*eij[0] + Si[1]*eij[1] + Si[2]*eij[2]
                  //       Sj . e_ij = Sj[0]*eij[0] + Sj[1]*eij[1] + Sj[2]*eij[2]
                  //
                  // Exchange tensor format
                  //
                  //     E = Si eijT Sj = eijx*eijx + ...
                  //
                  //---------------------------------------------------------------
                  //
                  tmp_tensor[0] -= Kij * eij[0]*eij[0]; // Kxx
                  tmp_tensor[1] -= Kij * eij[0]*eij[1]; // Kxy
                  tmp_tensor[2] -= Kij * eij[0]*eij[2]; // Kxz
                  tmp_tensor[3] -= Kij * eij[1]*eij[0]; // Kyx
                  tmp_tensor[4] -= Kij * eij[1]*eij[1]; // Kyx
                  tmp_tensor[5] -= Kij * eij[1]*eij[2]; // Kyz
                  tmp_tensor[6] -= Kij * eij[2]*eij[0]; // Kzx
                  tmp_tensor[7] -= Kij * eij[2]*eij[1]; // Kzy
                  tmp_tensor[8] -= Kij * eij[2]*eij[2]; // Kzy

               } // end of cutoff range if

            } // end of i!= j if

            //ofile << "-------------------" << std::endl;
            //std::cout << i << "\t" << nj << "\t" << imat << '\t' << jmat << '\t'<<
            //         tmp_tensor [0] << "\t" << tmp_tensor [1] << "\t" << tmp_tensor [2] << "\t" <<
            //         tmp_tensor [3] << "\t" << tmp_tensor [4] << "\t" << tmp_tensor [5] << "\t" <<
            //         tmp_tensor [6] << "\t" << tmp_tensor [7] << "\t" << tmp_tensor [8] << "\n" << std::endl;

            // save tensor for interaction i-j
            atoms::t_exchange_list[counter].Jij[0][0] += tmp_tensor [0];
            atoms::t_exchange_list[counter].Jij[0][1] += tmp_tensor [1];
            atoms::t_exchange_list[counter].Jij[0][2] += tmp_tensor [2];

            atoms::t_exchange_list[counter].Jij[1][0] += tmp_tensor [3];
            atoms::t_exchange_list[counter].Jij[1][1] += tmp_tensor [4];
            atoms::t_exchange_list[counter].Jij[1][2] += tmp_tensor [5];

            atoms::t_exchange_list[counter].Jij[2][0] += tmp_tensor [6];
            atoms::t_exchange_list[counter].Jij[2][1] += tmp_tensor [7];
            atoms::t_exchange_list[counter].Jij[2][2] += tmp_tensor [8];

            counter++; // increment interaction counter

         }

      } // end of atom loop

      //ofile.close();

      return;

   } // end of kitaev initialisation

} // end of internal namespace

} // end of exchange namespace
