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
//#include "atoms.hpp" // for exchange list type defs
#include "exchange.hpp"
#include "dipole.hpp"
#include <iostream>
// exchange module headers
#include "internal.hpp"

namespace exchange{

namespace internal{

   //-----------------------------------------------------------------------------
   // Function to calculate exchange fields for spins between start and end index
   //-----------------------------------------------------------------------------
   void exchange_fields(const int start_index, // first atom for exchange interactions to be calculated
                        const int end_index, // last +1 atom to be calculated
                        const std::vector<int>& neighbour_list_start_index,
                        const std::vector<int>& neighbour_list_end_index,
                        const std::vector<int>& type_array, // type for atom
                        const std::vector<int>& neighbour_list_array, // list of interactions between atoms
                        const std::vector<int>& neighbour_interaction_type_array, // list of interaction type for each pair of atoms with value given in exchange list
                        const std::vector <zval_t>& i_exchange_list, // list of isotropic exchange constants
                        const std::vector <zvec_t>& v_exchange_list, // list of vectorial exchange constants
                        const std::vector <zten_t>& t_exchange_list, // list of tensorial exchange constants
                        const std::vector<double>& spin_array_x, // spin vectors for atoms
                        const std::vector<double>& spin_array_y,
                        const std::vector<double>& spin_array_z,
                        std::vector<double>& field_array_x, // field vectors for atoms
                        std::vector<double>& field_array_y,
                        std::vector<double>& field_array_z){

   	// Use appropriate function for exchange calculation
   	switch(internal::exchange_type){

   		case exchange::isotropic:

            // loop over all atoms
   			for(int atom = start_index; atom < end_index; ++atom){

               // temporary variables (registers) to calculate intermediate sum
   				double hx = 0.0;
   				double hy = 0.0;
   				double hz = 0.0;

               // temporray constants for loop start and end indices
   				const int start = neighbour_list_start_index[atom];
   				const int end   = neighbour_list_end_index[atom]+1;

               // loop over all neighbours
   				for(int nn = start; nn < end; ++nn){

   					const int natom = neighbour_list_array[nn]; // get neighbouring atom number
   					const double Jij = i_exchange_list[ neighbour_interaction_type_array[nn] ].Jij; // get exchange constant between atoms

   					hx += Jij * spin_array_x[natom]; // add exchange fields
   					hy += Jij * spin_array_y[natom];
   					hz += Jij * spin_array_z[natom];

   				}

   				field_array_x[atom] += hx; // save total field to field array
   				field_array_y[atom] += hy;
   				field_array_z[atom] += hz;

   			}
   			break;

   		case exchange::vectorial: // vector

            // loop over all atoms
            for(int atom = start_index; atom < end_index; ++atom){

               // temporary variables (registers) to calculate intermediate sum
               double hx = 0.0;
               double hy = 0.0;
               double hz = 0.0;

               // temporray constants for loop start and end indices
               const int start = neighbour_list_start_index[atom];
               const int end   = neighbour_list_end_index[atom]+1;

               // loop over all neighbours
               for(int nn = start; nn < end; ++nn){

                  const int natom = neighbour_list_array[nn]; // get neighbouring atom number
                  const int iid = neighbour_interaction_type_array[nn]; // interaction id

   					const double Jij[3]={v_exchange_list[iid].Jij[0],
   												v_exchange_list[iid].Jij[1],
   												v_exchange_list[iid].Jij[2]};

                  hx += Jij[0] * spin_array_x[natom]; // add exchange fields
   					hy += Jij[1] * spin_array_y[natom];
   					hz += Jij[2] * spin_array_z[natom];

   				}

               field_array_x[atom] += hx; // save total field to field array
   				field_array_y[atom] += hy;
   				field_array_z[atom] += hz;

   			}
   			break;

   		case exchange::tensorial: // tensor

            // loop over all atoms
            for(int atom = start_index; atom < end_index; ++atom){

               // temporary variables (registers) to calculate intermediate sum
               double hx = 0.0;
               double hy = 0.0;
               double hz = 0.0;

               // temporray constants for loop start and end indices
               const int start = neighbour_list_start_index[atom];
               const int end   = neighbour_list_end_index[atom]+1;

               // loop over all neighbours
               for(int nn = start; nn < end; ++nn){

                  const int natom = neighbour_list_array[nn]; // get neighbouring atom number
                  const int iid = neighbour_interaction_type_array[nn]; // interaction id

   					const double Jij[3][3]={ {t_exchange_list[iid].Jij[0][0],
   													  t_exchange_list[iid].Jij[0][1],
   												     t_exchange_list[iid].Jij[0][2]},

   													 {t_exchange_list[iid].Jij[1][0],
   													  t_exchange_list[iid].Jij[1][1],
   													  t_exchange_list[iid].Jij[1][2]},

   													 {t_exchange_list[iid].Jij[2][0],
   													  t_exchange_list[iid].Jij[2][1],
   													  t_exchange_list[iid].Jij[2][2]} };

   					const double S[3]={spin_array_x[natom], spin_array_y[natom], spin_array_z[natom]};

   					hx += ( Jij[0][0] * S[0] + Jij[0][1] * S[1] + Jij[0][2] * S[2]);
   					hy += ( Jij[1][0] * S[0] + Jij[1][1] * S[1] + Jij[1][2] * S[2]);
   					hz += ( Jij[2][0] * S[0] + Jij[2][1] * S[1] + Jij[2][2] * S[2]);
   				}

               field_array_x[atom] += hx; // save total field to field array
   				field_array_y[atom] += hy;
   				field_array_z[atom] += hz;

   			}
   			break;

   		}

         // if (dipole::atomsitic_tensor_enabled){
         //
         //    for(int atom = start_index; atom < end_index; ++atom){
         //       // temporray constants for loop start and end indices
         //       const int start = dipole::atomistic_dd_neighbourlist_start[atom];
         //       const int end   = dipole::atomistic_dd_neighbourlist_end[atom]+1;
         //       std::cout << atom << '\t' << start << '\t' << end << std::endl;
         //       // loop over all neighbours
         //       for(int nn = start; nn < end; ++nn){
         //          int atomj = dipole::atomistic_dd_neighbourlist[nn];
         //       //   std::cout <<"NN:\t" << atom << '\t' << atomj << '\t' <<  std::endl;
         //       }
         //    }
         //
         //
         //    }




   		return;

   	}

} // end of internal namespace

} // end of exchange namespace
