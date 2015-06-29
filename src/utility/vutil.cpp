//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// System headers
#include <chrono>

// Program headers

//---------------------------------------------------------------------
// Namespace including assorted utility functions
//---------------------------------------------------------------------
namespace vutil{

   // simple class for performing code timing
   class vtimer_t{
      
   private:
      std::chrono::high_resolution_clock::time_point start_time;
      
   public:
      // start the timer
      void start(){
         start_time = std::chrono::high_resolution_clock::now();
      }
      
      // get the elapsed time in seconds
      double elapsed_time(){
         
         // get current time
         std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
         
         // work out elapsed time
         return double(chrono::duration_cast<chrono::seconds>(end_time - start_time).count());
         
      }
   }
   
} // end of namespace vutil

