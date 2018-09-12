//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>
#include <sstream>

// Vampire headers
#include "atoms.hpp" // to be removed
#include "create.hpp" // to be removed

#include "anisotropy.hpp"
#include "errors.hpp"
#include "units.hpp"
#include "vio.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

namespace internal{

   //---------------------------------------------------------------------------
   // Function to calculate surface anisotropy tensor
   //---------------------------------------------------------------------------
   void initialise_neel_anisotropy_tensor(std::vector <std::vector <bool> >& nearest_neighbour_interactions_list,
                                          std::vector<std::vector <neighbours::neighbour_t> >& cneighbourlist){

      // Print informative message to log file
      zlog << zTs() << "Using Néel pair anisotropy for atoms with < threshold number of neighbours." << std::endl;

      // allocate memory for neel anisotropy tensor
      internal::neel_tensor.resize( 9 * atoms::num_atoms, 0.0 );

      // temporary tensor for calculating sum
      std::vector<double> tmp_tensor(9,0.0);

      //	Populate surface atom and 1D nearest neighbour list and index arrays
      for(int atom=0;atom<atoms::num_atoms;atom++){


         // Only calculate parameters for atoms with less than full nn coordination
         if(atoms::surface_array[atom]){

            // get material id for atom
            const unsigned int imat = atoms::type_array[atom];

            // get inverse moment
            const double i_mu_s = 1.0/mp::material[imat].mu_s_SI;

            // zero tensor
            for(int idx = 0; idx < 9; idx++) tmp_tensor[idx] = 0.0;

            // loop over all neighbours
            for(unsigned int nn = 0; nn < cneighbourlist[atom].size(); nn++){

               // only add nearest neighbours to list
               if(nearest_neighbour_interactions_list[atom][nn]==true){

                  // get atom number for neighbour
                  const unsigned int natom = cneighbourlist[atom][nn].nn;

                  // get material id for j atom
                  const unsigned int jmat = atoms::type_array[natom];

                  // get atomic position vector i->j
                  double eij[3]={cneighbourlist[atom][nn].vx, cneighbourlist[atom][nn].vy, cneighbourlist[atom][nn].vz};

                  // normalise to unit vector
                  const double rij = sqrt(eij[0]*eij[0]+eij[1]*eij[1]+eij[2]*eij[2]);
                  const double invrij = 1.0/rij;

                  // normalise components to unit vector
                  eij[0] = eij[0] * invrij;
                  eij[1] = eij[1] * invrij;
                  eij[2] = eij[2] * invrij;

                  // get pair anisotropy constant between atom i and j in Tesla
                  double kij = anisotropy::internal::mp[imat].kij[jmat]*i_mu_s;

                  // Adjust kij by exponential factor if needed
                  if(internal::neel_range_dependent){
                     // Normalised exponential e(r0) = 1.0
                     kij = kij*exp(-neel_exponential_factor*( rij - neel_exponential_range )/neel_exponential_range);
                  };

                  // loop over tensor components and sum total including local anisotropy constant
                  // note inclusion of factor - 2 . 1/2 = -1 in tensor field compared to energy due to derivative
                  for(int i = 0; i < 3; ++i){
                     for(int j = 0; j < 3; ++j){
                        tmp_tensor[ 3*i + j ] += kij * eij[i] * eij[j];
                     }
                  }

                  //int natom = cneighbourlist[atom][nn].nn;
                  //std::cout << "nn_id: " << nn << " j: " << cneighbourlist[atom][nn].nn << "\trange: " << 1.0/invrij << " ";
                  //std::cout << "eij: " << eij[0] << " " << eij[1] << " " << eij[2] << "\tatomi: ";
                  //std::cout << catom_array[atom].x << " " << catom_array[atom].y << " " << catom_array[atom].z << "\tatomj: ";
                  //std::cout << catom_array[natom].x << " " << catom_array[natom].y << " " << catom_array[natom].z << std::endl;
                  //std::cout << atom << "\t" << natom << "\t" << kij << "\t" << anisotropy::internal::mp[imat].kij[jmat] << "\t" << i_mu_s << std::endl;

               }

            } // end of neighbour loop

            // save tensor for atom
            anisotropy::internal::neel_tensor[ 9 * atom + 0 ] = tmp_tensor [0];
            anisotropy::internal::neel_tensor[ 9 * atom + 1 ] = tmp_tensor [1];
            anisotropy::internal::neel_tensor[ 9 * atom + 2 ] = tmp_tensor [2];

            anisotropy::internal::neel_tensor[ 9 * atom + 3 ] = tmp_tensor [3];
            anisotropy::internal::neel_tensor[ 9 * atom + 4 ] = tmp_tensor [4];
            anisotropy::internal::neel_tensor[ 9 * atom + 5 ] = tmp_tensor [5];

            anisotropy::internal::neel_tensor[ 9 * atom + 6 ] = tmp_tensor [6];
            anisotropy::internal::neel_tensor[ 9 * atom + 7 ] = tmp_tensor [7];
            anisotropy::internal::neel_tensor[ 9 * atom + 8 ] = tmp_tensor [8];

         }
         //std::cout << std::endl;
      } // end of atom loop

      //---------------------------------------------------------------------------------
      // Output Neel tensors to file
      //---------------------------------------------------------------------------------

      /*std::ofstream ofile("Neel_tensor.txt");

      for(int atom=0; atom < atoms::num_atoms; atom++){
         // only output tensors for surface atoms
         if(atoms::surface_array[atom]){
            ofile << atom << "\t" << atoms::x_coord_array[atom] << "\t" << atoms::y_coord_array[atom] << "\t" << atoms::z_coord_array[atom] << "\t";
            for(int i=0; i<9; i++) ofile << anisotropy::internal::neel_tensor[ 9 * atom + i ] << "\t";
            ofile << std::endl;
         }
      }

      ofile.close();*/

   } // end of surface anisotropy initialisation

} // end of internal namespace

} // end of anisotropy namespace
