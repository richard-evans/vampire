#ifndef VUTIL_H_
#define VUTIL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
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
      std::chrono::high_resolution_clock::time_point end_time;

   public:
      // start the timer
      void start(){
         start_time = std::chrono::high_resolution_clock::now();
      }

      // start the timer
      void stop(){
         end_time = std::chrono::high_resolution_clock::now();
      }

      // get the elapsed time in milliseconds
      double elapsed_time(){

         // work out elapsed time
         return 1.e-9*double(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());

      }
   };

} // end of namespace vutil

#endif //VUTIL_H_
