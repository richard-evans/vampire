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
   void unroll_normalised_biquadratic_exchange_interactions(){

      // check that biquadratic interactions are enabled, if not do nothing
      if(!exchange::biquadratic) return;

   	// temporary class variables
   	value_t tmp_zval;
   	vector_t tmp_zvec;
   	tensor_t tmp_zten;

   	switch(exchange::internal::biquadratic_exchange_type){
   		case exchange::isotropic:
   			// unroll material calculations
   			std::cout << "Using generic/normalised form of exchange interaction with " << cs::unit_cell.biquadratic.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 1.0*double(exchange::internal::biquadratic_neighbour_list_array.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			exchange::internal::bq_i_exchange_list.reserve(exchange::internal::biquadratic_neighbour_list_array.size());
   			// loop over all interactions
   			for(int atom = 0; atom < atoms::num_atoms; atom++){
   				const int imaterial = atoms::type_array[atom];
               const double imus = 1.0 / mp::material[imaterial].mu_s_SI; // get inverse spin moment
   				for(int nn = exchange::internal::biquadratic_neighbour_list_start_index[atom]; nn <= exchange::internal::biquadratic_neighbour_list_end_index[atom]; nn++){
   					const int natom = exchange::internal::biquadratic_neighbour_list_array[nn];
   					const int jmaterial = atoms::type_array[natom];
   					exchange::internal::bq_i_exchange_list.push_back(tmp_zval);
                  // get unit cell interaction id
                  int i = exchange::internal::biquadratic_neighbour_interaction_type_array[nn];
                  // get shell ID for interaction
                  const int shell = cs::unit_cell.bilinear.interaction[i].shell;
                  // get exchange value from 4D exchange matrix
                  std::vector<double> Jij = internal::biquadratic_exchange_constants.get_exchange_values(imaterial, jmaterial, shell);
                  // set exchange field, normalising to mu_s^i
                  exchange::internal::bq_i_exchange_list[nn].Jij = cs::unit_cell.biquadratic.interaction[i].Jij[0][0] * Jij[0] * imus;
   					// reset interation id to neighbour number - causes segfault if nn out of range
   					exchange::internal::biquadratic_neighbour_interaction_type_array[nn] = nn;
   				}
   			}
   			break;

         case exchange::vectorial: // normalised vectorial exchange
      			// unroll material calculations
      			std::cout << "Using normalised vectorial form of exchange interaction with " << cs::unit_cell.biquadratic.interaction.size() << " total interactions." << std::endl;
      			zlog << zTs() << "Unrolled exchange template requires " << 3.0*double(exchange::internal::biquadratic_neighbour_list_array.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
      			exchange::internal::bq_v_exchange_list.reserve(exchange::internal::biquadratic_neighbour_list_array.size());
      			// loop over all interactions
      			for(int atom = 0; atom < atoms::num_atoms; atom++){
      				const int imaterial = atoms::type_array[atom];
                  const double imus = 1.0 / mp::material[imaterial].mu_s_SI; // get inverse spin moment
      				for(int nn = exchange::internal::biquadratic_neighbour_list_start_index[atom];nn <= exchange::internal::biquadratic_neighbour_list_end_index[atom]; nn++){
      					const int natom = exchange::internal::biquadratic_neighbour_list_array[nn];
      					const int jmaterial = atoms::type_array[natom];
      					exchange::internal::bq_v_exchange_list.push_back(tmp_zvec);
                     // get unit cell interaction id
                     int i = exchange::internal::biquadratic_neighbour_interaction_type_array[nn];
                     // get shell ID for interaction
                     const int shell = cs::unit_cell.bilinear.interaction[i].shell;
                     // get exchange value from 4D exchange matrix
                     std::vector<double> Jij = internal::biquadratic_exchange_constants.get_exchange_values(imaterial, jmaterial, shell);
                     // set exchange field, normalising to mu_s^i
                     if( Jij.size() == 3 ){
                        exchange::internal::bq_v_exchange_list[nn].Jij[0] = cs::unit_cell.biquadratic.interaction[i].Jij[0][0] * Jij[0] * imus;
                        exchange::internal::bq_v_exchange_list[nn].Jij[1] = cs::unit_cell.biquadratic.interaction[i].Jij[1][1] * Jij[1] * imus;
                        exchange::internal::bq_v_exchange_list[nn].Jij[2] = cs::unit_cell.biquadratic.interaction[i].Jij[2][2] * Jij[2] * imus;
                     }
                     else if( Jij.size() == 1 ){
                        exchange::internal::bq_v_exchange_list[nn].Jij[0] = cs::unit_cell.biquadratic.interaction[i].Jij[0][0] * Jij[0] * imus;
                        exchange::internal::bq_v_exchange_list[nn].Jij[1] = cs::unit_cell.biquadratic.interaction[i].Jij[1][1] * Jij[0] * imus;
                        exchange::internal::bq_v_exchange_list[nn].Jij[2] = cs::unit_cell.biquadratic.interaction[i].Jij[2][2] * Jij[0] * imus;
                     }
                     else{
                        std::cerr     << "Programmer error! Biquadratic exchange values size of " << Jij.size() << " must be 1 or 3 values. Exiting" << std::endl;
                        zlog << zTs() << "Programmer error! Biquadratic exchange values size of " << Jij.size() << " must be 1 or 3 values. Exiting" << std::endl;
                        err::vexit();
                     }
      					// reset interation id to neighbour number - causes segfault if nn out of range
      					exchange::internal::biquadratic_neighbour_interaction_type_array[nn] = nn;
      				}
      			}
      			break;

         case exchange::tensorial: // normalised tensorial exchange
         {
   			std::cout << "Using normalised tensorial form of exchange interaction with " << cs::unit_cell.biquadratic.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 9.0*double(exchange::internal::biquadratic_neighbour_list_array.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			// unroll isotopic interactions
   			exchange::internal::bq_t_exchange_list.reserve(exchange::internal::biquadratic_neighbour_list_array.size());

            // loop over all interactions
            for(int atom = 0; atom < atoms::num_atoms; atom++){
               const int imaterial = atoms::type_array[atom];
               const double imus = 1.0 / mp::material[imaterial].mu_s_SI; // get inverse spin moment
               for(int nn = exchange::internal::biquadratic_neighbour_list_start_index[atom]; nn <= exchange::internal::biquadratic_neighbour_list_end_index[atom]; nn++){
                  const int natom = exchange::internal::biquadratic_neighbour_list_array[nn]; // atom id of neighbour atom
                  const int jmaterial = atoms::type_array[natom]; // material of neighbour atom
                  exchange::internal::bq_t_exchange_list.push_back(tmp_zten);
                  // get unit cell interaction id
                  int i = exchange::internal::biquadratic_neighbour_interaction_type_array[nn];
                  // get shell ID for interaction
                  const int shell = cs::unit_cell.bilinear.interaction[i].shell;
                  // get exchange value from 4D exchange matrix
                  std::vector<double> Jij = internal::bilinear_exchange_constants.get_exchange_values(imaterial, jmaterial, shell);
                  // set exchange field, normalising to mu_s^i
                  // future development may allow for generic inclusion of DMI parameter from exchange tensor, but not currently enabled
                  if( Jij.size() == 3 ){
                     exchange::internal::bq_t_exchange_list[nn].Jij[0][0] = cs::unit_cell.biquadratic.interaction[i].Jij[0][0] * Jij[0] * imus;
                     exchange::internal::bq_t_exchange_list[nn].Jij[0][1] = cs::unit_cell.biquadratic.interaction[i].Jij[0][1] * 0.0; //mp::material[imaterial].Dij[jmaterial];
                     exchange::internal::bq_t_exchange_list[nn].Jij[0][2] = cs::unit_cell.biquadratic.interaction[i].Jij[0][2] * 0.0; //mp::material[imaterial].Dij[jmaterial];

                     exchange::internal::bq_t_exchange_list[nn].Jij[1][0] = cs::unit_cell.biquadratic.interaction[i].Jij[1][0] * 0.0; //mp::material[imaterial].Dij[jmaterial];
                     exchange::internal::bq_t_exchange_list[nn].Jij[1][1] = cs::unit_cell.biquadratic.interaction[i].Jij[1][1] * Jij[1] * imus;
                     exchange::internal::bq_t_exchange_list[nn].Jij[1][2] = cs::unit_cell.biquadratic.interaction[i].Jij[1][2] * 0.0; //mp::material[imaterial].Dij[jmaterial];

                     exchange::internal::bq_t_exchange_list[nn].Jij[2][0] = cs::unit_cell.biquadratic.interaction[i].Jij[2][0] * 0.0; //mp::material[imaterial].Dij[jmaterial];
                     exchange::internal::bq_t_exchange_list[nn].Jij[2][1] = cs::unit_cell.biquadratic.interaction[i].Jij[2][1] * 0.0; //mp::material[imaterial].Dij[jmaterial];
                     exchange::internal::bq_t_exchange_list[nn].Jij[2][2] = cs::unit_cell.biquadratic.interaction[i].Jij[2][2] * Jij[2] * imus;
                  }
                  else if( Jij.size() == 1 ){
                     exchange::internal::bq_t_exchange_list[nn].Jij[0][0] = cs::unit_cell.biquadratic.interaction[i].Jij[0][0] * Jij[0] * imus;
                     exchange::internal::bq_t_exchange_list[nn].Jij[0][1] = cs::unit_cell.biquadratic.interaction[i].Jij[0][1] * 0.0; //mp::material[imaterial].Dij[jmaterial];
                     exchange::internal::bq_t_exchange_list[nn].Jij[0][2] = cs::unit_cell.biquadratic.interaction[i].Jij[0][2] * 0.0; //mp::material[imaterial].Dij[jmaterial];

                     exchange::internal::bq_t_exchange_list[nn].Jij[1][0] = cs::unit_cell.biquadratic.interaction[i].Jij[1][0] * 0.0; //mp::material[imaterial].Dij[jmaterial];
                     exchange::internal::bq_t_exchange_list[nn].Jij[1][1] = cs::unit_cell.biquadratic.interaction[i].Jij[1][1] * Jij[0] * imus;
                     exchange::internal::bq_t_exchange_list[nn].Jij[1][2] = cs::unit_cell.biquadratic.interaction[i].Jij[1][2] * 0.0; //mp::material[imaterial].Dij[jmaterial];

                     exchange::internal::bq_t_exchange_list[nn].Jij[2][0] = cs::unit_cell.biquadratic.interaction[i].Jij[2][0] * 0.0; //mp::material[imaterial].Dij[jmaterial];
                     exchange::internal::bq_t_exchange_list[nn].Jij[2][1] = cs::unit_cell.biquadratic.interaction[i].Jij[2][1] * 0.0; //mp::material[imaterial].Dij[jmaterial];
                     exchange::internal::bq_t_exchange_list[nn].Jij[2][2] = cs::unit_cell.biquadratic.interaction[i].Jij[2][2] * Jij[0] * imus;
                  }
                  else{
                     std::cerr     << "Programmer error! Exchange values size of " << Jij.size() << " must be 1 or 3 values. Exiting" << std::endl;
                     zlog << zTs() << "Programmer error! Exchange values size of " << Jij.size() << " must be 1 or 3 values. Exiting" << std::endl;
                     err::vexit();
                  }
                  // reset interation id to neighbour number - causes segfault if nn out of range
                  exchange::internal::biquadratic_neighbour_interaction_type_array[nn] = nn;

               } // end of neighbour loop
   			} // end of atom loop
         }
   		break;

   		default:
   			terminaltextcolor(RED);
   			std::cerr << "Error! - Unknown unit cell biquadratic exchange type " << internal::biquadratic_exchange_type << "; unable to unroll exchenge template. Exiting" << std::endl;
   			terminaltextcolor(WHITE);
   			err::vexit();
   			break;
   	}

      return;

   }

} // end of internal namespace

} // end of exchange namespace
