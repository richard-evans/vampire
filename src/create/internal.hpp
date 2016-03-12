#ifndef CREATE_INTERNAL_H_
#define CREATE_INTERNAL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// create implementation. These functions should
// not be accessed outside of the create module.
//---------------------------------------------------------------------

// transitional arrangement - need old create header for catom_t
// Should eventually move class definition to this file
#include "create.hpp"

// Vampire headers
#include "material.hpp"
#include "mtrand.hpp"

namespace create{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Internal shared variables used for creation
      //-----------------------------------------------------------------------------
      // materials class for storing material parameters
      class mp_t{

      private:

      public:
         // variables
         bool alloy_master; // flag specifying if material is host material
         std::vector<double> alloy_fraction;
         // constructor
         mp_t ():
         	alloy_master(false){
               // resize array of alloy substitution fractions
               alloy_fraction.resize(mp::max_materials,0.0);
            };
      };

      extern std::vector<create::internal::mp_t> mp; // array of material properties
      extern MTRand grnd; // general random number generator for create functions

      // Internal functions for create module
      extern void alloy(std::vector<cs::catom_t> & catom_array);

   } // end of internal namespace
} // end of create namespace

#endif //CREATE_INTERNAL_H_
