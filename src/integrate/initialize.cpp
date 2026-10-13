//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "integrate.hpp"
#include "material.hpp"

// integrate module headers
#include "internal.hpp"

namespace integrate{

   //----------------------------------------------------------------------------
   // Function to initialize integrate module
   //----------------------------------------------------------------------------
   void initialize(){

      const int num_materials = mp::num_materials;

      // Unroll LSF coefficient arrays
      internal::lsf_second_order_coefficient.resize(num_materials,0.0);
      internal::lsf_fourth_order_coefficient.resize(num_materials,0.0);
      internal::lsf_sixth_order_coefficient.resize(num_materials,0.0);

      // loop over materials set by user
      for(unsigned int m=0; m < internal::mp.size(); ++m){

         // Assign Landau coefficients based on material
         if(internal::mp[m].lsf_second_order_coefficient.is_set()) internal::lsf_second_order_coefficient[m] = internal::mp[m].lsf_second_order_coefficient.get();
         if(internal::mp[m].lsf_fourth_order_coefficient.is_set()) internal::lsf_fourth_order_coefficient[m] = internal::mp[m].lsf_fourth_order_coefficient.get();
         if(internal::mp[m].lsf_sixth_order_coefficient.is_set()) internal::lsf_sixth_order_coefficient[m] = internal::mp[m].lsf_sixth_order_coefficient.get();

      }

      // Check and log rotation based LLG integrators
      internal::initialize_rotation_integrators();

      return;

   }

} // end of integrate namespace

