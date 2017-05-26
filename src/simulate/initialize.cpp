//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
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

      // unroll slonczewski torque arrays
      sim::internal::slonczewski_aj.resize(num_materials,0.0);
      sim::internal::slonczewski_bj.resize(num_materials,0.0);
      // loop over materials set by user
      for(unsigned int m=0; m<sim::internal::mp.size(); ++m){
         // copy values set by user to arrays
         if(sim::internal::mp[m].slonczewski_aj.is_set()) sim::internal::slonczewski_aj[m] = sim::internal::mp[m].slonczewski_aj.get();
         if(sim::internal::mp[m].slonczewski_bj.is_set()) sim::internal::slonczewski_bj[m] = sim::internal::mp[m].slonczewski_bj.get();
      }

      return;
   }

} // end of namespace sim
