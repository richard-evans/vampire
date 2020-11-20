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
#include "atoms.hpp" // for exchange list type defs
#include "exchange.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

   //-----------------------------------------------------------------------------
   // Function to calculate exchange fields for spins between start and end index
   //-----------------------------------------------------------------------------
   void fields(const int start_index, // first atom for exchange interactions to be calculated
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


   	// Calculate standard (bilinear) exchange fields
      exchange::internal::exchange_fields(start_index, end_index,
                                neighbour_list_start_index, neighbour_list_end_index,
                                type_array, neighbour_list_array, neighbour_interaction_type_array,
                                i_exchange_list, v_exchange_list, t_exchange_list,
                                spin_array_x, spin_array_y, spin_array_z,
                                field_array_x, field_array_y, field_array_z);

      // calculate biquadratic exchange field
      if(exchange::biquadratic){
         exchange::internal::biquadratic_exchange_fields(start_index, end_index,
                                                         exchange::internal::biquadratic_neighbour_list_start_index, exchange::internal::biquadratic_neighbour_list_end_index,
                                                         type_array, exchange::internal::biquadratic_neighbour_list_array, exchange::internal::biquadratic_neighbour_interaction_type_array,
                                                         internal::bq_i_exchange_list, internal::bq_v_exchange_list, internal::bq_t_exchange_list,
                                                         spin_array_x, spin_array_y, spin_array_z,
                                                         field_array_x, field_array_y, field_array_z);
      }

   	return;

   }

} // end of exchange namespace
