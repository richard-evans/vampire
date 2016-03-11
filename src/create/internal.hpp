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

#include "mtrand.hpp"

namespace create{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Internal shared variables used for creation
      //-----------------------------------------------------------------------------
      // simple initialised class for set variables
      class set_double_t{

      private:
         double value; // value
         bool setf; // flag specifiying variable has been set

      public:
         // class functions
         // constructor
         set_double_t() : value(0.0), setf(false) { }

         // setting function
         void set(double in_value){
            value = in_value;
            setf = true;
         };

         // get value function
         double get(){ return value; };
         // check if variable is set
         bool is_set(){ return setf; };

      };

      struct mp_t{
         set_double_t slonczewski_aj; // adiabatic torque
         set_double_t slonczewski_bj; // field torque (non-adiabatic)
      };

      extern std::vector<create::internal::mp_t> mp; // array of material properties
      extern MTRand grnd; // general random number generator for create functions

      // Internal functions for create module
      extern void alloy(std::vector<cs::catom_t> & catom_array);

   } // end of internal namespace
} // end of create namespace

#endif //CREATE_INTERNAL_H_
