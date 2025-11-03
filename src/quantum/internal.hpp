//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2025. All rights reserved.
//
//   Email: fried-conrad.weber@uni-potsdam.de
//
//------------------------------------------------------------------------------
//

#ifndef QUANTUM_INTERNAL_H_
#define QUANTUM_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the quantum module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "quantum.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-----------------------------------------------------------------------------
      // internal materials class for storing material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

          private:

          public:

             //------------------------------
             // material parameter variables
             //------------------------------
             double test;

             // constructor
             mp_t (const unsigned int max_materials = 100):
                test(0.0) // constructor initialisation of test variable
             {
                // constructor body for initialising more complex data/arrays
             }; // end of constructor

       }; // end of internal::mp class

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      extern bool enabled; // bool to enable module

      extern std::vector<internal::mp_t> mp; // array of material properties

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

   } // end of internal namespace

} // end of quantum namespace

#endif //QUANTUM_INTERNAL_H_
