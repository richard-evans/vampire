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

namespace anisotropy
{
    //----------------------------------------------------------------------------
    // Function to calculate uniaxial fields
    //----------------------------------------------------------------------------
    void calculate_anisotropy_fields (
        const int num_atoms,                        // atoms in system
        const std::vector<int>& atom_type_array,
        const std::vector<double>& atom_coords_x,
        const std::vector<double>& atom_coords_y,
        const std::vector<double>& atom_coords_z)
        {
            if (anisotropy::uniaxial) calculate_uniaxial_fields(num_atoms);
        }

        return;

    }

} // end of anisotropy namespace
