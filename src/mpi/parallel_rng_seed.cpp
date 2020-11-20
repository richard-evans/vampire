//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "vmpi.hpp"
#include "mtrand.hpp"

// Internal vmpi header

namespace vmpi{

   //------------------------------------------------------------------------------
   // Function to generate parallel safe random seed on multiple processors.
   //
   // The source seed generates the same random number on each processor which is
   // then offset by num_processors to ensure each processor has a unique seed
   // (avoids the birthday problem).
   //
   // Implementation uses language defined wraparound behaviour of unsigned int.
   //------------------------------------------------------------------------------
   uint32_t parallel_rng_seed(int seed){

      // instatiate random number generator class
      MTRand_int32 grnd(seed);

      // number to store randomised seed
      uint32_t seed_number;

      // generate a random number through 1000 iterations of seed
      for(int i = 0; i<1000; i++) seed_number = grnd();

      // may want to broadcast seed_number here to force same number on all processors

      // start with long seed to account for high core counts
      uint64_t long_seed = static_cast<uint64_t>(seed_number) + static_cast<uint64_t>(vmpi::my_rank) * static_cast<uint64_t>(vmpi::num_processors);

      // truncate seed to 32 bit integer with wrap around
      uint32_t short_seed = static_cast<uint32_t>(long_seed);

      // return final parallel seed
      return short_seed;

   }

}
