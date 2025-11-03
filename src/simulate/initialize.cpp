//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2020. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "material.hpp"
#include "sim.hpp"
#include "internal.hpp"

namespace sim{

   //-------------------------------------------------------------------------------
   // initialise sim namespace variables
   //-------------------------------------------------------------------------------
   void initialize(int num_materials){

      // unroll slonczewski spin transfer torque arrays
      sim::internal::stt_asm.resize(num_materials,0.0);
      sim::internal::stt_rj.resize(num_materials,0.0);
      sim::internal::stt_pj.resize(num_materials,0.0);

      // unroll spin orbit torque arrays
      sim::internal::sot_asm.resize(num_materials,0.0);
      sim::internal::sot_rj.resize(num_materials,0.0);
      sim::internal::sot_pj.resize(num_materials,0.0);

      sim::internal::vcmak.resize(num_materials, 0.0);

      // Unroll LSF coefficient arrays
      sim::internal::lsf_second_order_coefficient.resize(num_materials,0.0);
      sim::internal::lsf_fourth_order_coefficient.resize(num_materials,0.0);
      sim::internal::lsf_sixth_order_coefficient.resize(num_materials,0.0);

      // loop over materials set by user
      for(unsigned int m=0; m < sim::internal::mp.size(); ++m){
         // copy values set by user to arrays
         if(sim::internal::mp[m].stt_asm.is_set()) sim::internal::stt_asm[m] = sim::internal::mp[m].stt_asm.get();
         if(sim::internal::mp[m].stt_rj.is_set())  sim::internal::stt_rj[m]  = sim::internal::mp[m].stt_rj.get();
         if(sim::internal::mp[m].stt_pj.is_set())  sim::internal::stt_pj[m]  = sim::internal::mp[m].stt_pj.get();

         if(sim::internal::mp[m].sot_asm.is_set()) sim::internal::sot_asm[m] = sim::internal::mp[m].sot_asm.get();
         if(sim::internal::mp[m].sot_rj.is_set())  sim::internal::sot_rj[m]  = sim::internal::mp[m].sot_rj.get();
         if(sim::internal::mp[m].sot_pj.is_set())  sim::internal::sot_pj[m]  = sim::internal::mp[m].sot_pj.get();

         // set vcma coefficients (requires sim::internal::enable_vcma_fields == true) but this should be default
         if(sim::internal::mp[m].vcmak.is_set()){
            const double imu_s = 1.0 / mp::material[m].mu_s_SI; // calculate inverse moment
            sim::internal::vcmak[m] = imu_s * sim::internal::mp[m].vcmak.get();
         }

         // Assign Landau coefficients based on material
         if(sim::internal::mp[m].lsf_second_order_coefficient.is_set()) sim::internal::lsf_second_order_coefficient[m] = sim::internal::mp[m].lsf_second_order_coefficient.get();
         if(sim::internal::mp[m].lsf_fourth_order_coefficient.is_set()) sim::internal::lsf_fourth_order_coefficient[m] = sim::internal::mp[m].lsf_fourth_order_coefficient.get();
         if(sim::internal::mp[m].lsf_sixth_order_coefficient.is_set()) sim::internal::lsf_sixth_order_coefficient[m] = sim::internal::mp[m].lsf_sixth_order_coefficient.get();

      }

      return;
   }

} // end of namespace gpu
