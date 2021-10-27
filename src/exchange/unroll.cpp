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
   // Function to unroll neighbour list into 1D
   //----------------------------------------------------------------------------
   void unroll_exchange_interactions(std::vector<std::vector <neighbours::neighbour_t> >& bilinear){

      // if dmi is enabled then set exchange type to force normalised tensor form of exchange
      if(internal::enable_dmi || internal::enable_kitaev){
         internal::exchange_type = exchange::tensorial;
         internal::use_material_exchange_constants = true;
      }

      // if normalised exchange is to be used (taking values from exchange matrix
      // defined in the material file), then unroll normalised interactions
      if(internal::use_material_exchange_constants){
         unroll_normalised_exchange_interactions(bilinear);
         return;
      }

   	// temporary class variables
   	zval_t tmp_zval;
   	zvec_t tmp_zvec;
   	zten_t tmp_zten;

   	switch(internal::exchange_type){

   		case exchange::isotropic:
   			std::cout << "Using isotropic form of exchange interaction with " << cs::unit_cell.bilinear.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 1.0*double(cs::unit_cell.bilinear.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			// unroll isotopic interactions
   			atoms::i_exchange_list.reserve(cs::unit_cell.bilinear.interaction.size());
   			for(unsigned int i=0;i<cs::unit_cell.bilinear.interaction.size();i++){
   				int iatom = cs::unit_cell.bilinear.interaction[i].i;
   				int imat = cs::unit_cell.atom[iatom].mat;
   				atoms::i_exchange_list.push_back(tmp_zval);
   				atoms::i_exchange_list[i].Jij = cs::unit_cell.bilinear.interaction[i].Jij[0][0]/mp::material[imat].mu_s_SI;
   			}
   			break;

   		case exchange::vectorial:
   			std::cout << "Using vectorial form of exchange interaction with " << cs::unit_cell.bilinear.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 3.0*double(cs::unit_cell.bilinear.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			// unroll isotopic interactions
   			atoms::v_exchange_list.reserve(cs::unit_cell.bilinear.interaction.size());
   			for(unsigned int i=0;i<cs::unit_cell.bilinear.interaction.size();i++){
   				int iatom = cs::unit_cell.bilinear.interaction[i].i;
   				int imat = cs::unit_cell.atom[iatom].mat;
   				atoms::v_exchange_list.push_back(tmp_zvec);
   				atoms::v_exchange_list[i].Jij[0] = cs::unit_cell.bilinear.interaction[i].Jij[0][0]/mp::material[imat].mu_s_SI;
   				atoms::v_exchange_list[i].Jij[1] = cs::unit_cell.bilinear.interaction[i].Jij[1][1]/mp::material[imat].mu_s_SI;
   				atoms::v_exchange_list[i].Jij[2] = cs::unit_cell.bilinear.interaction[i].Jij[2][2]/mp::material[imat].mu_s_SI;
   			}
   			break;

   		case exchange::tensorial:
   			std::cout << "Using tensorial form of exchange interaction with " << cs::unit_cell.bilinear.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 9.0*double(cs::unit_cell.bilinear.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			// unroll isotopic interactions
   			atoms::t_exchange_list.reserve(cs::unit_cell.bilinear.interaction.size());
   			for(unsigned int i=0;i<cs::unit_cell.bilinear.interaction.size();i++){
   				int iatom = cs::unit_cell.bilinear.interaction[i].i;
   				int imat = cs::unit_cell.atom[iatom].mat;
   				atoms::t_exchange_list.push_back(tmp_zten);

   				atoms::t_exchange_list[i].Jij[0][0] = cs::unit_cell.bilinear.interaction[i].Jij[0][0]/mp::material[imat].mu_s_SI;
   				atoms::t_exchange_list[i].Jij[0][1] = cs::unit_cell.bilinear.interaction[i].Jij[0][1]/mp::material[imat].mu_s_SI;
   				atoms::t_exchange_list[i].Jij[0][2] = cs::unit_cell.bilinear.interaction[i].Jij[0][2]/mp::material[imat].mu_s_SI;

   				atoms::t_exchange_list[i].Jij[1][0] = cs::unit_cell.bilinear.interaction[i].Jij[1][0]/mp::material[imat].mu_s_SI;
   				atoms::t_exchange_list[i].Jij[1][1] = cs::unit_cell.bilinear.interaction[i].Jij[1][1]/mp::material[imat].mu_s_SI;
   				atoms::t_exchange_list[i].Jij[1][2] = cs::unit_cell.bilinear.interaction[i].Jij[1][2]/mp::material[imat].mu_s_SI;

   				atoms::t_exchange_list[i].Jij[2][0] = cs::unit_cell.bilinear.interaction[i].Jij[2][0]/mp::material[imat].mu_s_SI;
   				atoms::t_exchange_list[i].Jij[2][1] = cs::unit_cell.bilinear.interaction[i].Jij[2][1]/mp::material[imat].mu_s_SI;
   				atoms::t_exchange_list[i].Jij[2][2] = cs::unit_cell.bilinear.interaction[i].Jij[2][2]/mp::material[imat].mu_s_SI;
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
