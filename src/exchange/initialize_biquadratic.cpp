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
#include "errors.hpp"
#include "exchange.hpp"
#include "unitcell.hpp"
#include "vio.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

namespace internal{

   void initialize_biquadratic_exchange(){

      // if biquadratic exchange is not needed then do nothing
      if(!exchange::biquadratic) return;

      // if normalised exchange is to be used (taking values from exchange matrix
      // defined in the material file), then unroll normalised interactions
      if(internal::use_material_biquadratic_exchange_constants){
         unroll_normalised_biquadratic_exchange_interactions();
         return;
      }

      // temporary class variables
   	value_t  tmp_val;
   	vector_t tmp_vec;
      tensor_t tmp_ten;

   	switch(internal::biquadratic_exchange_type){

   		case exchange::isotropic:
            std::cout << "Initialising biquadratic exchange interactions." << std::endl;
            zlog << zTs() << "Initialising biquadratic exchange interactions." << std::endl;
            zlog << zTs() << "Unrolled biquadratic exchange template requires "
                 << 1.0*double(cs::unit_cell.biquadratic.interaction.size())*double(sizeof(double))*1.0e-6
                 << "MB RAM" << std::endl;

   			// unroll isotopic interactions
   			exchange::internal::bq_i_exchange_list.reserve(cs::unit_cell.biquadratic.interaction.size());
   			for(unsigned int i=0;i<cs::unit_cell.biquadratic.interaction.size();i++){

   				int iatom = cs::unit_cell.biquadratic.interaction[i].i;
   				int imat = cs::unit_cell.atom[iatom].mat;

   				exchange::internal::bq_i_exchange_list.push_back(tmp_val);

   				exchange::internal::bq_i_exchange_list[i].Jij = cs::unit_cell.biquadratic.interaction[i].Jij[0][0]/mp::material[imat].mu_s_SI;

   			}
   			break;

   		case exchange::vectorial:
   			std::cout << "Using vectorial form of biquadratic exchange interaction with " << cs::unit_cell.biquadratic.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled biquadratic exchange template requires " << 3.0*double(cs::unit_cell.biquadratic.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;

   			// unroll vectorial interactions
   			exchange::internal::bq_v_exchange_list.reserve(cs::unit_cell.biquadratic.interaction.size());
   			for(unsigned int i=0;i<cs::unit_cell.biquadratic.interaction.size();i++){

   				int iatom = cs::unit_cell.biquadratic.interaction[i].i;
   				int imat = cs::unit_cell.atom[iatom].mat;

   				exchange::internal::bq_v_exchange_list.push_back(tmp_vec);

   				exchange::internal::bq_v_exchange_list[i].Jij[0] = cs::unit_cell.biquadratic.interaction[i].Jij[0][0]/mp::material[imat].mu_s_SI;
   				exchange::internal::bq_v_exchange_list[i].Jij[1] = cs::unit_cell.biquadratic.interaction[i].Jij[1][1]/mp::material[imat].mu_s_SI;
   				exchange::internal::bq_v_exchange_list[i].Jij[2] = cs::unit_cell.biquadratic.interaction[i].Jij[2][2]/mp::material[imat].mu_s_SI;

   			}
   			break;

   		case exchange::tensorial:

   			std::cout << "Using tensorial form of biquadratic exchange interaction with " << cs::unit_cell.biquadratic.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled biquadratic exchange template requires " << 9.0*double(cs::unit_cell.biquadratic.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;

            // unroll tensorial interactions
   			exchange::internal::bq_t_exchange_list.reserve(cs::unit_cell.biquadratic.interaction.size());
   			for(unsigned int i=0;i<cs::unit_cell.biquadratic.interaction.size();i++){

   				int iatom = cs::unit_cell.biquadratic.interaction[i].i;
   				int imat = cs::unit_cell.atom[iatom].mat;

   				exchange::internal::bq_t_exchange_list.push_back(tmp_ten);

   				exchange::internal::bq_t_exchange_list[i].Jij[0][0] = cs::unit_cell.biquadratic.interaction[i].Jij[0][0]/mp::material[imat].mu_s_SI;
   				exchange::internal::bq_t_exchange_list[i].Jij[0][1] = cs::unit_cell.biquadratic.interaction[i].Jij[0][1]/mp::material[imat].mu_s_SI;
   				exchange::internal::bq_t_exchange_list[i].Jij[0][2] = cs::unit_cell.biquadratic.interaction[i].Jij[0][2]/mp::material[imat].mu_s_SI;

   				exchange::internal::bq_t_exchange_list[i].Jij[1][0] = cs::unit_cell.biquadratic.interaction[i].Jij[1][0]/mp::material[imat].mu_s_SI;
   				exchange::internal::bq_t_exchange_list[i].Jij[1][1] = cs::unit_cell.biquadratic.interaction[i].Jij[1][1]/mp::material[imat].mu_s_SI;
   				exchange::internal::bq_t_exchange_list[i].Jij[1][2] = cs::unit_cell.biquadratic.interaction[i].Jij[1][2]/mp::material[imat].mu_s_SI;

   				exchange::internal::bq_t_exchange_list[i].Jij[2][0] = cs::unit_cell.biquadratic.interaction[i].Jij[2][0]/mp::material[imat].mu_s_SI;
   				exchange::internal::bq_t_exchange_list[i].Jij[2][1] = cs::unit_cell.biquadratic.interaction[i].Jij[2][1]/mp::material[imat].mu_s_SI;
   				exchange::internal::bq_t_exchange_list[i].Jij[2][2] = cs::unit_cell.biquadratic.interaction[i].Jij[2][2]/mp::material[imat].mu_s_SI;

   			}
   			break;

   		default:
   			terminaltextcolor(RED);
   			std::cerr << "Error! - Unknown unit cell biquadratic exchange type " << internal::biquadratic_exchange_type << "; unable to unroll exchange template. Exiting" << std::endl;
   			terminaltextcolor(WHITE);
   			err::vexit();
   			break;

   	}

      return;

   }

} // end of internal namespace

} // end of exchange namespace
