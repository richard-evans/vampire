//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2014. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "gpu.hpp"
#include "cuda.hpp"
#include "errors.hpp"
#include "program.hpp"

namespace gpu{

   //---------------------------------------------------------------------------
   // Function to determine whether the program may read the spins on the host
   // directly during the simulation
   //---------------------------------------------------------------------------
   bool program_reads_host_spins(){

      // programs which only use statistics of the spins during the simulation
      switch(::program::program){
         case ::program::program_t::time_series:
         case ::program::program_t::hysteresis_loop:
         case ::program::program_t::curie_temperature:
         case ::program::program_t::field_cool:
         case ::program::program_t::laser_pulse:
         case ::program::program_t::hamr_simulation:
         case ::program::program_t::partial_hysteresis_loop:
         case ::program::program_t::fmr:
            return false;
         default:
            return true;
      }

   }

   namespace stats{

      //-------------------------------------------------------------------------------
      // Function to call correct statistics update function on device
      //-------------------------------------------------------------------------------
      void update(){

         #ifdef CUDA
            vcuda::stats::update();
         #endif

         return;
      }

      //-------------------------------------------------------------------------------
      // Function to get statistics module data from device for output
      //
      // Remember to use the set_magnetization in the magnetization_statistic_t instances
      //
      // There is a performance penalty for rapid updates of statistics for output
      // to disk.
      //
      //-------------------------------------------------------------------------------
      void get(){

         #ifdef CUDA
            vcuda::stats::get();
         #endif

         return;
      }

      //-------------------------------------------------------------------------------
      // Function to reset statistics counters on device
      //-------------------------------------------------------------------------------
      void reset(){

         #ifdef CUDA
            vcuda::stats::reset();
         #endif

         return;
      }

   } // end of stats namespace

} // end of namespace gpu
