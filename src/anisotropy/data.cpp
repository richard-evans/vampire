//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: sw766@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

    //------------------------------------------------------------------------------
    // Externally visible variables
    //------------------------------------------------------------------------------

    bool enabled = false;

    namespace internal{

        //------------------------------------------------------------------------
        // Shared variables inside anisotropy module
        //------------------------------------------------------------------------

        bool initialised = false;

        bool uniaxial = false;

        std::vector<double> field_array;

    } // end of internal namespace

} // end of anisotropy namespace
