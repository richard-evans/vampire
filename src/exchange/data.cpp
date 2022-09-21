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
#include "exchange.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------
   bool biquadratic = false; // flag to enable biquadratic exchange calculation

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside exchange module
      //------------------------------------------------------------------------
      std::vector<internal::mp_t> mp; // array of material properties

      exchange_matrix_4D_t bilinear_exchange_constants; // array of exchange constants
      exchange_matrix_4D_t biquadratic_exchange_constants; // array of biquadratic exchange constants

      bool enable_dmi = false; // flag to enable dmi calculation
      bool enable_kitaev = false; // flag to enable Kitaev calculation

      double dmi_cutoff_range = 2.6; // cutoff range for DMI calculation (Ångstroms)
      double kitaev_cutoff_range = 2.6; // cutoff range for Kitaev calculation (Ångstroms)
      double exchange_factor = 1.0; // scaling factor for exchange constants (usually to correct for ab-initio)

      exchange_t exchange_type = isotropic; // exchange type to use in simulation
      exchange_t biquadratic_exchange_type = isotropic; // biquadratic exchange type to use in simulation

      exchange_t minimum_needed_exchange_type = isotropic; // minimum support required for bilinear exchange type (for vectorial constants and DMI)

      bool use_material_exchange_constants = true; // flag to enable material exchange parameters
      bool use_material_biquadratic_exchange_constants = true; // flag to enable material biquadratic exchange parameters

      std::vector <int> biquadratic_neighbour_list_array; // 1D list of biquadratic neighbours
      std::vector <int> biquadratic_neighbour_interaction_type_array; // 1D list of biquadratic exchange interaction types
      std::vector <int> biquadratic_neighbour_list_start_index; // list of first biquadratic neighbour for atom i
      std::vector <int> biquadratic_neighbour_list_end_index;   // list of last biquadratic neighbour for atom i

      std::vector <exchange::internal::value_t  > bq_i_exchange_list(0); // list of isotropic biquadratic exchange constants
      std::vector <exchange::internal::vector_t > bq_v_exchange_list(0); // list of vectorial biquadratic exchange constants
      std::vector <exchange::internal::tensor_t > bq_t_exchange_list(0); // list of tensorial biquadratic exchange constants

   } // end of internal namespace

} // end of exchange namespace
