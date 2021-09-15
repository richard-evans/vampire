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
#include "sim.hpp"
#include "internal.hpp"

namespace sim{

   //-------------------------------------------------------------------------------
   // Function to call correct initialization function depending on cuda,opencl etc
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

      // loop over materials set by user
      for(unsigned int m=0; m<sim::internal::mp.size(); ++m){
         // copy values set by user to arrays
         if(sim::internal::mp[m].stt_asm.is_set()) sim::internal::stt_asm[m] = sim::internal::mp[m].stt_asm.get();
         if(sim::internal::mp[m].stt_rj.is_set())  sim::internal::stt_rj[m]  = sim::internal::mp[m].stt_rj.get();
         if(sim::internal::mp[m].stt_pj.is_set())  sim::internal::stt_pj[m]  = sim::internal::mp[m].stt_pj.get();

         if(sim::internal::mp[m].sot_asm.is_set()) sim::internal::sot_asm[m] = sim::internal::mp[m].sot_asm.get();
         if(sim::internal::mp[m].sot_rj.is_set())  sim::internal::sot_rj[m]  = sim::internal::mp[m].sot_rj.get();
         if(sim::internal::mp[m].sot_pj.is_set())  sim::internal::sot_pj[m]  = sim::internal::mp[m].sot_pj.get();
      }

      return;
   }

} // end of namespace gpu
