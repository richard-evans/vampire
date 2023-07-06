//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
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
   // Function to calculate dmi vectors
   //
   //                                    k
   //                                    o
   //                                  /   \
   //                                 /     \
   //                                O ----- O
   //                                i       j
   //
   // For neighbours i and j, the DMI is a mediated interaction between i and j
   // via k. The strength of the interaction is the sum of Dik and Djk where they
   // may be different if i and j are of different material types.
   //
   // The DMI vector is defined as the cross prouct rik x rjk when both are
   // within their respective cutoff ranges for i-k and j-k interactions.
   //
   //------------------------------------------------------------------------------
   void calculate_dmi(std::vector<std::vector <neighbours::neighbour_t> >& cneighbourlist){

      // if dmi is not needed then do nothing
      if(!internal::enable_dmi) return;

      // Print informative message to log file
      zlog << zTs() << "Calculating Dzyaloshinskii-Moriya interactions" << std::endl;

      // temporary tensor for calculating sum
      std::vector<double> tmp_tensor(9,0.0);

      // get cutoff range
      const double cutoff_sq = internal::dmi_cutoff_range*internal::dmi_cutoff_range;

      // counter for number of interactions i-j
      unsigned int counter = 0;

      // counter for total number of interactions
      unsigned int total_counter = 0;

      //std::ofstream ofile("dmi.txt");

      //	Loop over all atoms i
      for(unsigned int i=0; i < static_cast<unsigned int>(atoms::num_atoms); i++){

         // get material id for atom
         const unsigned int imat = atoms::type_array[i];

         // get inverse moment
         const double i_mu_s = 1.0/mp::material[imat].mu_s_SI;

         // loop over all neighbours j
         for(unsigned int j = 0; j < cneighbourlist[i].size(); j++){

            // get atom number for neighbour i
            const unsigned int nj = cneighbourlist[i][j].nn;

            // get material id for j atom
            const unsigned int jmat = atoms::type_array[nj];

            // reinitialise dmi tensor for each i-j interaction
            for(int idx = 0; idx < 9; idx++) tmp_tensor[idx] = 0.0;

            // ignore self interaction
            if(i != nj){
               // for each interaction j loop over all neighbours k to calculate
               // mediated interactions within cutoff range
               for(unsigned int k = 0; k < cneighbourlist[i].size(); k++){

                  // get atom number for neighbour k
                  const unsigned int nk = cneighbourlist[i][k].nn;

                  // ignore self interaction
                  if(nj != nk){

                     // get material id for k atom
                     const unsigned int kmat = atoms::type_array[nk];

                     // get atomic position vector i->k
                     double eik[3]={cneighbourlist[i][k].vx, cneighbourlist[i][k].vy, cneighbourlist[i][k].vz};
                     const double mod_eik_sq = eik[0]*eik[0] + eik[1]*eik[1] + eik[2]*eik[2];

                     // get atomic position vector i->j
                     double eij[3]={cneighbourlist[i][j].vx, cneighbourlist[i][j].vy, cneighbourlist[i][j].vz};

                     // calculate ejk from vector addition eik - eij
                     double ejk[3]={eik[0] - eij[0], eik[1] - eij[1], eik[2] - eij[2]};
                     const double mod_ejk_sq = ejk[0]*ejk[0] + ejk[1]*ejk[1] + ejk[2]*ejk[2];

                     // check if both ik and jk are closer than cutoff range and i != j
                     // before computing DMI
                     if( mod_eik_sq <= cutoff_sq && mod_ejk_sq <= cutoff_sq ){

                        // normalise to unit vector
                        const double inv_rik=1.0/sqrt(mod_eik_sq);
                        const double inv_rjk=1.0/sqrt(mod_ejk_sq);

                        // normalise components to unit vector
                        for(int idx = 0; idx < 3; idx++){
                           eik[idx] = eik[idx] * inv_rik;
                           ejk[idx] = ejk[idx] * inv_rjk;
                        }

                        // Set pair dmi constant average between i=->k and j->k normalised to mu_s_i
                        const double Dij = 0.5*( exchange::internal::mp[imat].dmi[kmat] +
                                                 exchange::internal::mp[jmat].dmi[kmat]) * i_mu_s;

                        // compute direction of Dij
                        const double Dx = eik[1]*ejk[2] - eik[2]*ejk[1];
                        const double Dy = eik[2]*ejk[0] - eik[0]*ejk[2];
                        const double Dz = eik[0]*ejk[1] - eik[1]*ejk[0];

                        //ofile << i << "\t" << nj << "\t" << nk << "\t" << sqrt(mod_eik_sq) << "\t" << sqrt(mod_ejk_sq) << std::endl <<
                        //           "\t\teik: " << eik[0] << "\t" << eik[1] << "\t" << eik[2] << std::endl <<
                        //           "\t\teij: " << eij[0] << "\t" << eij[1] << "\t" << eij[2] << std::endl <<
                        //           "\t\tejk: " << ejk[0] << "\t" << ejk[1] << "\t" << ejk[2] << std::endl;
                        //ofile << i << "\t" << nj << "\t" << nk << "\t" << sqrt(mod_eik_sq) << "\t" << sqrt(mod_ejk_sq) << "\t" << Dx << "\t" << Dy << "\t" << Dz << std::endl;

                        // increment total interaction counter
                        total_counter++;

                        //---------------------------------------------------------------
                        //        Expression of DMI within the exchange tensor
                        //---------------------------------------------------------------
                        //
                        //       E = Dij . (Si x Sj)
                        //
                        // Si x Sj = Si[1]*Sj[2] - Si[2]*Sj[1] = [ Siy Sjz - Siz Sjy ]
                        //           Si[2]*Sj[0] - Si[0]*Sj[2]   [ Siz Sjx - Six Sjz ]
                        //           Si[0]*Sj[1] - Si[1]*Sj[0]   [ Six Sjy - Siy Sjx ]
                        //
                        //       E = Dx [ Siy Sjz - Siz Sjy ] +
                        //           Dy [ Siz Sjx - Six Sjz ] +
                        //           Dz [ Six Sjy - Siy Sjx ]
                        //
                        // Exchange tensor format
                        //
                        //     E = Si J Sj = Six Jxx Sjx + Six Jxy Sjy + Six Jxz Sjz +
                        //                   Siy Jyx Sjx + Siy Jyy Sjy + Siy Jyz Sjz +
                        //                   Siz Jzx Sjx + Siz Jzy Sjy + Siz Jzz Sjz
                        //
                        //  Adding dmi components to the exchange tensor
                        //
                        //       Jxy = J[0][1] = +Dz
                        //       Jyx = J[1][0] = -Dz      [  0  +Dz  -Dy ]
                        //       Jzx = J[2][0] = +Dy   =  [ -Dz  0   +Dx ]
                        //       Jxz = J[0][2] = -Dy      [ +Dy  -Dx  0  ]
                        //       Jyz = J[1][2] = +Dx
                        //       Jzy = J[2][1] = -Dx
                        //
                        //---------------------------------------------------------------
                        //
                        tmp_tensor[1] -= Dij * +Dz; // Jxy
                        tmp_tensor[2] -= Dij * -Dy; // Jxz
                        tmp_tensor[3] -= Dij * -Dz; // Jyx
                        tmp_tensor[5] -= Dij * +Dx; // Jyz
                        tmp_tensor[6] -= Dij * +Dy; // Jzx
                        tmp_tensor[7] -= Dij * -Dx; // Jzy

                        //std::cout << i << "\t" << nj << "\t" << nk << "\t|\t" << Dx << "\t" << Dy << "\t" << Dz << "\t" << Dij << "\t|\t" <<
                        //   " eik: " << eik[0] << "\t" << eik[1] << "\t" << eik[2] << "\t|\t" <<
                        //   " ejk: " << ejk[0] << "\t" << ejk[1] << "\t" << ejk[2] << "\t|\t" <<
                        //   tmp_tensor [0] << "\t" << Dij * +Dz      << "\t" << Dij * -Dy << "\t" <<
                        //   Dij * -Dz      << "\t" << tmp_tensor [4] << "\t" << Dij * +Dx << "\t" <<
                        //   Dij * +Dy << "\t" << Dij * -Dx << "\t" << tmp_tensor [8] << "\n" << std::flush;

                     } // end of cutoff range if

                  } // end of k != j if

               } // end of k neighbour loop

            } // end of i!= j if

            // zero almost zero components for symmetry
            for(int i=0 ; i<9; i++){
               if(fabs(tmp_tensor[i]) < 1.0e-12) tmp_tensor[i] = 0.0;
            }

            // dmi exchange tensor for each atom pair
            //ofile << i << "\t" << nj << "\t" << imat << '\t' << jmat << '\t'<<
            //         tmp_tensor [0] << "\t" << tmp_tensor [1] << "\t" << tmp_tensor [2] << "\t" <<
            //         tmp_tensor [3] << "\t" << tmp_tensor [4] << "\t" << tmp_tensor [5] << "\t" <<
            //         tmp_tensor [6] << "\t" << tmp_tensor [7] << "\t" << tmp_tensor [8] << "\t" << std::endl;

            // net DMI vectors for each atom pair
            //ofile << i << "\t" << nj << "\t" << imat << '\t' << jmat << '\t'<< tmp_tensor [5] << "\t" << tmp_tensor [6] << "\t" << tmp_tensor [1] << std::endl;

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

      zlog << zTs() << "Generated " << total_counter << " dmi interactions with an average of " << double(total_counter) / double(atoms::num_atoms) << " interactions per atom" << std::endl;

      if(total_counter == 0){
         std::cerr << "Warning! DMI is enabled but not interactions were found. The built-in DMI is therefore not working - try increasing the dmi interaction range." << std::endl;
         zlog << zTs() << "Warning! DMI is enabled but not interactions were found. The built-in DMI is therefore not working - try increasing the dmi interaction range." << std::endl;
      }

      return;

   } // end of dmi initialisation

} // end of internal namespace

} // end of exchange namespace
