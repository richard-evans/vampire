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
#include "config.hpp"

// config module headers
#include "internal.hpp"

namespace config {

   namespace internal {

      //------------------------------------------------------------------------
      // Function to copy data to temporary buffer for output
      // (needs to be templated for float output)
      //------------------------------------------------------------------------
      void copy_data_to_buffer(const std::vector<double> &x, // vector data
                               const std::vector<double> &y,
                               const std::vector<double> &z,
                               const std::vector<uint64_t> &mask,
                               std::vector<double> &buffer){

         // copy total number of output data to const for compiler
         const unsigned int data_size = mask.size();

         // loop over all atoms to be output
         for (unsigned int id = 0; id < data_size; ++id){

            // determine next datum to be output
            const unsigned int index = mask[id];

            // copy data to be output to main output buffer
            buffer[3 * id + 0] = x[index];
            buffer[3 * id + 1] = y[index];
            buffer[3 * id + 2] = z[index];

         }

         return;

      }

   } // end of internal namespace

} // end of config namespace
