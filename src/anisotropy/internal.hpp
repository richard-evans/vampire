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

#ifndef ANISOTROPY_INTERNAL_H_
#define ANISOTROPY_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the anisotropy module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "material.hpp"
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

    extern bool enabled; // enable anisotropy calculation

    extern int internally_calculate_anisotropy_fields(const int num_atoms);

    namespace internal{

        //-------------------------------------------------------------------------
        // internal data type definitions
        //-------------------------------------------------------------------------
        extern bool initialised; // check module has been initialised

        extern bool uniaxial;

        //-------------------------------------------------------------------------
        // Internal shared variables
        //-------------------------------------------------------------------------

        extern int num_atoms;
        extern std::vector<int> atom_type_array;
        extern std::vector<zkval_t> materialscalaranisotropyarray;

        extern std::vector<double> spin_array_x;
        extern std::vector<double> spin_array_y;
        extern std::vector<double> spin_array_z;

        extern std::vector<double> field_array;

        //-------------------------------------------------------------------------
        // internal function declarations
        //-------------------------------------------------------------------------

    } // end of internal namespace

} // end of anisotropy namespace

#endif //ANISOTROPY_INTERNAL_H_
