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
#include "ltmp.hpp"
#include "vmpi.hpp"

// Local temperature pulse headers
#include "internal.hpp"

namespace ltmp{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Function to calculate the local temperature using the two temperature model
      //-----------------------------------------------------------------------------
      void calculate_local_temperature_gradient(){

         // Precalculate temperature gradient
         const double Tmin = ltmp::internal::minimum_temperature;
         const double Tmax = ltmp::internal::maximum_temperature;

         // Calculate new electron and lattice temperatures with temperature gradient
         for(unsigned int cell=0; cell<ltmp::internal::attenuation_array.size(); ++cell){

            // Determine cell temperature
            const double sqrtT = sqrt(Tmin + Tmax*attenuation_array[cell]);

            // Assume Te = Tp = T and save
            root_temperature_array[2*cell+0] = sqrtT;
            root_temperature_array[2*cell+1] = sqrtT;

         }

         // optionally output cell data
         if(ltmp::internal::output_microcell_data) ltmp::internal::write_cell_temperature_data();

         return;

      }

   } // end of namespace internal
} // end of namespace ltmp
