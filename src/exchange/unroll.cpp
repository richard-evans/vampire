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
   // Function to initialize exchange module
   //----------------------------------------------------------------------------
   void unroll_exchange_interactions(){

   	// condense interaction list
   	atoms::exchange_type=cs::unit_cell.exchange_type;

   	// temporary class variables
   	zval_t tmp_zval;
   	zvec_t tmp_zvec;
   	zten_t tmp_zten;

   	switch(atoms::exchange_type){
   		case -1:
   			// unroll material calculations
   			std::cout << "Using generic form of exchange interaction with " << cs::unit_cell.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 1.0*double(atoms::neighbour_list_array.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			atoms::i_exchange_list.reserve(atoms::neighbour_list_array.size());
   			// loop over all interactions
   			for(int atom=0;atom<atoms::num_atoms;atom++){
   				const int imaterial=atoms::type_array[atom];
   				for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){
   					const int natom = atoms::neighbour_list_array[nn];
   					const int jmaterial=atoms::type_array[natom];
   					atoms::i_exchange_list.push_back(tmp_zval);
                  // get unit cell interaction id
                  int i = atoms::neighbour_interaction_type_array[nn];
                  atoms::i_exchange_list[nn].Jij=cs::unit_cell.interaction[i].Jij[0][0]*mp::material[imaterial].Jij_matrix[jmaterial][0];
   					// reset interation id to neighbour number - causes segfault if nn out of range
   					atoms::neighbour_interaction_type_array[nn]=nn;
   				}
   			}
   			// now set exchange type to normal isotropic case
   			atoms::exchange_type=0;
   			break;
   		case 0:
   			std::cout << "Using isotropic form of exchange interaction with " << cs::unit_cell.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 1.0*double(cs::unit_cell.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			// unroll isotopic interactions
   			atoms::i_exchange_list.reserve(cs::unit_cell.interaction.size());
   			for(unsigned int i=0;i<cs::unit_cell.interaction.size();i++){
   				int iatom = cs::unit_cell.interaction[i].i;
   				int imat = cs::unit_cell.atom[iatom].mat;
   				atoms::i_exchange_list.push_back(tmp_zval);
   				atoms::i_exchange_list[i].Jij=-cs::unit_cell.interaction[i].Jij[0][0]/mp::material[imat].mu_s_SI;
   			}
   			break;
   		case 1:
   			std::cout << "Using vectorial form of exchange interaction with " << cs::unit_cell.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 3.0*double(cs::unit_cell.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			// unroll isotopic interactions
   			atoms::v_exchange_list.reserve(cs::unit_cell.interaction.size());
   			for(unsigned int i=0;i<cs::unit_cell.interaction.size();i++){
   				int iatom = cs::unit_cell.interaction[i].i;
   				int imat = cs::unit_cell.atom[iatom].mat;
   				atoms::v_exchange_list.push_back(tmp_zvec);
   				atoms::v_exchange_list[i].Jij[0]=-cs::unit_cell.interaction[i].Jij[0][0]/mp::material[imat].mu_s_SI;
   				atoms::v_exchange_list[i].Jij[1]=-cs::unit_cell.interaction[i].Jij[1][1]/mp::material[imat].mu_s_SI;
   				atoms::v_exchange_list[i].Jij[2]=-cs::unit_cell.interaction[i].Jij[2][2]/mp::material[imat].mu_s_SI;
   			}
   			break;
   		case 2:
   			std::cout << "Using tensorial form of exchange interaction with " << cs::unit_cell.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 9.0*double(cs::unit_cell.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			// unroll isotopic interactions
   			atoms::t_exchange_list.reserve(cs::unit_cell.interaction.size());
   			for(unsigned int i=0;i<cs::unit_cell.interaction.size();i++){
   				int iatom = cs::unit_cell.interaction[i].i;
   				int imat = cs::unit_cell.atom[iatom].mat;
   				atoms::t_exchange_list.push_back(tmp_zten);

   				atoms::t_exchange_list[i].Jij[0][0]=-cs::unit_cell.interaction[i].Jij[0][0]/mp::material[imat].mu_s_SI;
   				atoms::t_exchange_list[i].Jij[0][1]=-cs::unit_cell.interaction[i].Jij[0][1]/mp::material[imat].mu_s_SI;
   				atoms::t_exchange_list[i].Jij[0][2]=-cs::unit_cell.interaction[i].Jij[0][2]/mp::material[imat].mu_s_SI;

   				atoms::t_exchange_list[i].Jij[1][0]=-cs::unit_cell.interaction[i].Jij[1][0]/mp::material[imat].mu_s_SI;
   				atoms::t_exchange_list[i].Jij[1][1]=-cs::unit_cell.interaction[i].Jij[1][1]/mp::material[imat].mu_s_SI;
   				atoms::t_exchange_list[i].Jij[1][2]=-cs::unit_cell.interaction[i].Jij[1][2]/mp::material[imat].mu_s_SI;

   				atoms::t_exchange_list[i].Jij[2][0]=-cs::unit_cell.interaction[i].Jij[2][0]/mp::material[imat].mu_s_SI;
   				atoms::t_exchange_list[i].Jij[2][1]=-cs::unit_cell.interaction[i].Jij[2][1]/mp::material[imat].mu_s_SI;
   				atoms::t_exchange_list[i].Jij[2][2]=-cs::unit_cell.interaction[i].Jij[2][2]/mp::material[imat].mu_s_SI;
   			}
   			break;

         case 3: // normalised vectorial exchange
      			// unroll material calculations
      			std::cout << "Using vectorial form of exchange interaction with " << cs::unit_cell.interaction.size() << " total interactions." << std::endl;
      			zlog << zTs() << "Unrolled exchange template requires " << 3.0*double(atoms::neighbour_list_array.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
      			atoms::v_exchange_list.reserve(atoms::neighbour_list_array.size());
      			// loop over all interactions
      			for(int atom=0;atom<atoms::num_atoms;atom++){
      				const int imaterial=atoms::type_array[atom];
      				for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){
      					const int natom = atoms::neighbour_list_array[nn];
      					const int jmaterial=atoms::type_array[natom];
      					atoms::v_exchange_list.push_back(tmp_zvec);
                     // get unit cell interaction id
                     int i = atoms::neighbour_interaction_type_array[nn];
                     atoms::v_exchange_list[nn].Jij[0]=cs::unit_cell.interaction[i].Jij[0][0]*mp::material[imaterial].Jij_matrix[jmaterial][0];
                     atoms::v_exchange_list[nn].Jij[1]=cs::unit_cell.interaction[i].Jij[1][1]*mp::material[imaterial].Jij_matrix[jmaterial][1];
                     atoms::v_exchange_list[nn].Jij[2]=cs::unit_cell.interaction[i].Jij[2][2]*mp::material[imaterial].Jij_matrix[jmaterial][2];
      					// reset interation id to neighbour number - causes segfault if nn out of range
      					atoms::neighbour_interaction_type_array[nn]=nn;
      				}
      			}
      			// now set exchange type to normal vectorial case
      			atoms::exchange_type=1;
      			break;

         case 4: // normalised tensorial exchange
   			std::cout << "Using tensorial form of exchange interaction with " << cs::unit_cell.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 9.0*double(atoms::neighbour_list_array.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			// unroll isotopic interactions
   			atoms::t_exchange_list.reserve(atoms::neighbour_list_array.size());
            // loop over all interactions
            for(int atom=0;atom<atoms::num_atoms;atom++){
               const int imaterial=atoms::type_array[atom];
               for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){

                  const int natom = atoms::neighbour_list_array[nn];
                  const int jmaterial=atoms::type_array[natom];

                  atoms::t_exchange_list.push_back(tmp_zten);

                  // get unit cell interaction id
                  int i = atoms::neighbour_interaction_type_array[nn];

                  atoms::t_exchange_list[nn].Jij[0][0]=-cs::unit_cell.interaction[i].Jij[0][0] * mp::material[imaterial].Jij_matrix[jmaterial][0];
                  atoms::t_exchange_list[nn].Jij[0][1]=-cs::unit_cell.interaction[i].Jij[0][1] * 0.0; //mp::material[imaterial].Dij[jmaterial];
                  atoms::t_exchange_list[nn].Jij[0][2]=-cs::unit_cell.interaction[i].Jij[0][2] * 0.0; //mp::material[imaterial].Dij[jmaterial];

                  atoms::t_exchange_list[nn].Jij[1][0]=-cs::unit_cell.interaction[i].Jij[1][0] * 0.0; //mp::material[imaterial].Dij[jmaterial];
                  atoms::t_exchange_list[nn].Jij[1][1]=-cs::unit_cell.interaction[i].Jij[1][1] * mp::material[imaterial].Jij_matrix[jmaterial][1];
                  atoms::t_exchange_list[nn].Jij[1][2]=-cs::unit_cell.interaction[i].Jij[1][2] * 0.0; //mp::material[imaterial].Dij[jmaterial];

                  atoms::t_exchange_list[nn].Jij[2][0]=-cs::unit_cell.interaction[i].Jij[2][0] * 0.0; //mp::material[imaterial].Dij[jmaterial];
                  atoms::t_exchange_list[nn].Jij[2][1]=-cs::unit_cell.interaction[i].Jij[2][1] * 0.0; //mp::material[imaterial].Dij[jmaterial];
                  atoms::t_exchange_list[nn].Jij[2][2]=-cs::unit_cell.interaction[i].Jij[2][2] * mp::material[imaterial].Jij_matrix[jmaterial][2];

                  // reset interation id to neighbour number - causes segfault if nn out of range
                  atoms::neighbour_interaction_type_array[nn]=nn;

               } // end of neighbour loop
   			} // end of atom loop
   			break;

   		default:
   			terminaltextcolor(RED);
   			std::cerr << "Error! - Unknown unit cell exchange type " << atoms::exchange_type << "; unable to unroll exchenge template. Exiting" << std::endl;
   			terminaltextcolor(WHITE);
   			err::vexit();
   			break;
   	}

      return;

   }

} // end of internal namespace

} // end of exchange namespace
