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

namespace anisotropy {

    //----------------------------------------------------------------------------
    // function to calculate anisotropy fields
    //----------------------------------------------------------------------------
    int calculate_energies()
    {
        double energy =
            second_order_tensor_anisotropy() +
            third_order_tensor_anisotropy() +
            fourth_order_tensor_anisotropy();

        return EXIT_SUCCESS;
    }

} // end of anisotropy namespace
