//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "spininitialize.hpp"

// spininitialize module headers
#include "internal.hpp"

namespace spininitialize{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside spininitialize module
      //------------------------------------------------------------------------

      bool enabled; // bool to enable module

      std::vector<internal::mp_t> mp; // array of material spin initialisation properties

      std::vector<std::string> vector_field_filenames; // cache of loaded vector field file names
      std::vector< std::vector<field_point_t> > vector_field_data; // cache of loaded vector field data points

   } // end of internal namespace

} // end of spininitialize namespace
