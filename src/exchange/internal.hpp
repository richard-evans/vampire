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

#ifndef EXCHANGE_INTERNAL_H_
#define EXCHANGE_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the exchange module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "exchange.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------
      //-----------------------------------------------------------------------------
      // materials class for storing exchange material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

         private:

         public:

            // variables
            std::vector<double> dmi; // Dzyaloshinskii-Moriya interaction constant

            // constructor
            mp_t (const unsigned int max_materials = 100)
            {
               // resize arrays to correct size
               dmi.resize(max_materials, 0.0); // initialise pair anisotropy constants to zero

            }; // end of constructor

      }; // end of exchange::internal::mp class

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
      extern std::vector<internal::mp_t> mp; // array of material properties

      extern bool enable_dmi; // flag to enable dmi calculation

      extern double dmi_cutoff_range; // cutoff range for DMI calculation (Ångstroms)

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      void calculate_dmi(std::vector<std::vector <cs::neighbour_t> >& cneighbourlist);
      void unroll_exchange_interactions();

   } // end of internal namespace

} // end of exchange namespace

#endif //EXCHANGE_INTERNAL_H_
