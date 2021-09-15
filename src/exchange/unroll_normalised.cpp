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
#include "errors.hpp"
#include "exchange.hpp"
#include "unitcell.hpp"
#include "vio.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

   namespace internal{

   //----------------------------------------------------------------------------
   // Function to initialize exchange interactions assuming normalised values
   // of the exchange coupling from the unit cell which are then multiplied
   // by the material exchange matrix
   //
   // This requires additional memory since each interaction is potentially
   // unique, requiring that the whole exchange list be unrolled
   //----------------------------------------------------------------------------
   void unroll_normalised_exchange_interactions(std::vector<std::vector <neighbours::neighbour_t> >& bilinear2){

   	// temporary class variables
   	zval_t tmp_zval;
   	zvec_t tmp_zvec;
   	zten_t tmp_zten;

   	switch(exchange::internal::exchange_type){
   		case exchange::isotropic:
   			// unroll material calculations
   			std::cout << "Using generic/normalised form of exchange interaction with " << cs::unit_cell.bilinear.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 1.0*double(atoms::neighbour_list_array.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			atoms::i_exchange_list.reserve(atoms::neighbour_list_array.size());
   			// loop over all interactions
            for(int atom = 0; atom < atoms::num_atoms; atom++){
   				const int imaterial = atoms::type_array[atom];
               const double imus = 1.0 / mp::material[imaterial].mu_s_SI; // get inverse spin moment
   				for(int nn = atoms::neighbour_list_start_index[atom]; nn <= atoms::neighbour_list_end_index[atom]; nn++){
   					const int natom = atoms::neighbour_list_array[nn];
   					const int jmaterial = atoms::type_array[natom];
   					atoms::i_exchange_list.push_back(tmp_zval);
                  // get unit cell interaction id
                  int i = atoms::neighbour_interaction_type_array[nn];
                  // get shell ID for interaction
                  const int shell = cs::unit_cell.bilinear.interaction[i].shell;
                  // get exchange value from 4D exchange matrix
                  std::vector<double> Jij = internal::bilinear_exchange_constants.get_exchange_values(imaterial, jmaterial, shell);
                  // set exchange field, normalising to mu_s^i
                  atoms::i_exchange_list[nn].Jij = cs::unit_cell.bilinear.interaction[i].Jij[0][0] * Jij[0] * imus;
   					// reset interation id to neighbour number - causes segfault if nn out of range
   					atoms::neighbour_interaction_type_array[nn] = nn;
   				}
   			}
   			break;

         case exchange::vectorial: // normalised vectorial exchange
      			// unroll material calculations
      			std::cout << "Using normalised vectorial form of exchange interaction with " << cs::unit_cell.bilinear.interaction.size() << " total interactions." << std::endl;
      			zlog << zTs() << "Unrolled exchange template requires " << 3.0*double(atoms::neighbour_list_array.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
      			atoms::v_exchange_list.reserve(atoms::neighbour_list_array.size());
      			// loop over all interactions
      			for(int atom = 0; atom < atoms::num_atoms; atom++){
      				const int imaterial = atoms::type_array[atom];
                  const double imus = 1.0 / mp::material[imaterial].mu_s_SI; // get inverse spin moment
      		      for(int nn = atoms::neighbour_list_start_index[atom];nn <= atoms::neighbour_list_end_index[atom]; nn++){
      					const int natom = atoms::neighbour_list_array[nn];
      					const int jmaterial = atoms::type_array[natom];
      					atoms::v_exchange_list.push_back(tmp_zvec);
                     // get unit cell interaction id
                     int i = atoms::neighbour_interaction_type_array[nn];
                     // get shell ID for interaction
                     const int shell = cs::unit_cell.bilinear.interaction[i].shell;
                     // get exchange value from 4D exchange matrix
                     std::vector<double> Jij = internal::bilinear_exchange_constants.get_exchange_values(imaterial, jmaterial, shell);
                     // set exchange field, normalising to mu_s^i
                     if( Jij.size() == 3 ){
                        atoms::v_exchange_list[nn].Jij[0] = cs::unit_cell.bilinear.interaction[i].Jij[0][0] * Jij[0] * imus;
                        atoms::v_exchange_list[nn].Jij[1] = cs::unit_cell.bilinear.interaction[i].Jij[1][1] * Jij[1] * imus;
                        atoms::v_exchange_list[nn].Jij[2] = cs::unit_cell.bilinear.interaction[i].Jij[2][2] * Jij[2] * imus;
                     }
                     else if( Jij.size() == 1 ){
                        atoms::v_exchange_list[nn].Jij[0] = cs::unit_cell.bilinear.interaction[i].Jij[0][0] * Jij[0] * imus;
                        atoms::v_exchange_list[nn].Jij[1] = cs::unit_cell.bilinear.interaction[i].Jij[1][1] * Jij[0] * imus;
                        atoms::v_exchange_list[nn].Jij[2] = cs::unit_cell.bilinear.interaction[i].Jij[2][2] * Jij[0] * imus;
                     }
                     else{
                        std::cerr     << "Programmer error! Exchange values size of " << Jij.size() << " must be 1 or 3 values. Exiting" << std::endl;
                        zlog << zTs() << "Programmer error! Exchange values size of " << Jij.size() << " must be 1 or 3 values. Exiting" << std::endl;
                        err::vexit();
                     }
                   //  if (imaterial == 4) std::cout << imaterial << '\t' << jmaterial << "\t" << d << '\t' << J  << "\t" << atoms::v_exchange_list[nn].Jij[0] <<  std::endl;
      					// reset interaction id to neighbour number - causes segfault if nn out of range
      					atoms::neighbour_interaction_type_array[nn] = nn;
      				}
      			}
      			break;

         case exchange::tensorial: // normalised tensorial exchange
         {
   			std::cout << "Using normalised tensorial form of exchange interaction with " << cs::unit_cell.bilinear.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 9.0*double(atoms::neighbour_list_array.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			// unroll isotopic interactions
   			atoms::t_exchange_list.reserve(atoms::neighbour_list_array.size());

            // loop over all interactions
            for(int atom = 0; atom < atoms::num_atoms; atom++){
               const int imaterial = atoms::type_array[atom];
               const double imus = 1.0 / mp::material[imaterial].mu_s_SI; // get inverse spin moment
               for(int nn = atoms::neighbour_list_start_index[atom]; nn <= atoms::neighbour_list_end_index[atom]; nn++){

                  const int natom = atoms::neighbour_list_array[nn]; // atom id of neighbour atom
                  const int jmaterial = atoms::type_array[natom]; // material of neighbour atom

                  atoms::t_exchange_list.push_back(tmp_zten);

                  // get unit cell interaction id
                  int i = atoms::neighbour_interaction_type_array[nn];
                  // get shell ID for interaction
                  const int shell = cs::unit_cell.bilinear.interaction[i].shell;
                  // get exchange value from 4D exchange matrix
                  std::vector<double> Jij = internal::bilinear_exchange_constants.get_exchange_values(imaterial, jmaterial, shell);
                  // set exchange field, normalising to mu_s^i
                  // future development may allow for generic inclusion of DMI parameter from exchange tensor, but not currently enabled
                  if( Jij.size() == 3 ){
                     atoms::t_exchange_list[nn].Jij[0][0] = cs::unit_cell.bilinear.interaction[i].Jij[0][0] * Jij[0] * imus;
                     atoms::t_exchange_list[nn].Jij[0][1] = cs::unit_cell.bilinear.interaction[i].Jij[0][1] * 0.0;
                     atoms::t_exchange_list[nn].Jij[0][2] = cs::unit_cell.bilinear.interaction[i].Jij[0][2] * 0.0;
                     atoms::t_exchange_list[nn].Jij[1][0] = cs::unit_cell.bilinear.interaction[i].Jij[1][0] * 0.0;
                     atoms::t_exchange_list[nn].Jij[1][1] = cs::unit_cell.bilinear.interaction[i].Jij[1][1] * Jij[1] * imus;
                     atoms::t_exchange_list[nn].Jij[1][2] = cs::unit_cell.bilinear.interaction[i].Jij[1][2] * 0.0;
                     atoms::t_exchange_list[nn].Jij[2][0] = cs::unit_cell.bilinear.interaction[i].Jij[2][0] * 0.0;
                     atoms::t_exchange_list[nn].Jij[2][1] = cs::unit_cell.bilinear.interaction[i].Jij[2][1] * 0.0;
                     atoms::t_exchange_list[nn].Jij[2][2] = cs::unit_cell.bilinear.interaction[i].Jij[2][2] * Jij[2] * imus;
                  }
                  else if( Jij.size() == 1 ){
                     atoms::t_exchange_list[nn].Jij[0][0] = cs::unit_cell.bilinear.interaction[i].Jij[0][0] * Jij[0] * imus;
                     atoms::t_exchange_list[nn].Jij[0][1] = cs::unit_cell.bilinear.interaction[i].Jij[0][1] * 0.0;
                     atoms::t_exchange_list[nn].Jij[0][2] = cs::unit_cell.bilinear.interaction[i].Jij[0][2] * 0.0;
                     atoms::t_exchange_list[nn].Jij[1][0] = cs::unit_cell.bilinear.interaction[i].Jij[1][0] * 0.0;
                     atoms::t_exchange_list[nn].Jij[1][1] = cs::unit_cell.bilinear.interaction[i].Jij[1][1] * Jij[0] * imus;
                     atoms::t_exchange_list[nn].Jij[1][2] = cs::unit_cell.bilinear.interaction[i].Jij[1][2] * 0.0;
                     atoms::t_exchange_list[nn].Jij[2][0] = cs::unit_cell.bilinear.interaction[i].Jij[2][0] * 0.0;
                     atoms::t_exchange_list[nn].Jij[2][1] = cs::unit_cell.bilinear.interaction[i].Jij[2][1] * 0.0;
                     atoms::t_exchange_list[nn].Jij[2][2] = cs::unit_cell.bilinear.interaction[i].Jij[2][2] * Jij[0] * imus;
                  }
                  else{
                     std::cerr     << "Programmer error! Exchange values size of " << Jij.size() << " must be 1 or 3 values. Exiting" << std::endl;
                     zlog << zTs() << "Programmer error! Exchange values size of " << Jij.size() << " must be 1 or 3 values. Exiting" << std::endl;
                     err::vexit();
                  }

                  // reset interation id to neighbour number - causes segfault if nn out of range
                  atoms::neighbour_interaction_type_array[nn] = nn;

               } // end of neighbour loop
   			} // end of atom loop
         }
   		break;

   		default:
   			terminaltextcolor(RED);
   			std::cerr << "Error! - Unknown unit cell exchange type " << internal::exchange_type << "; unable to unroll exchenge template. Exiting" << std::endl;
   			terminaltextcolor(WHITE);
   			err::vexit();
   			break;
   	}

      return;

   }

} // end of internal namespace

} // end of exchange namespace
