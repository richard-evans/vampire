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

    //-------------------------------------------------------------------------
    // function declarations
    //-------------------------------------------------------------------------
    extern int calculate_energies();
    extern double second_order_tensor_anisotropy();
    extern double third_order_tensor_anisotropy();
    extern double fourth_order_tensor_anisotropy();

    namespace internal{

        //-------------------------------------------------------------------------
        // internal data type definitions
        //-------------------------------------------------------------------------
        extern bool initialised; // check module has been initialised

        /* anisotropy type flags */
        extern bool uniaxial_first_order;
        extern bool neel;

        extern bool uniaxial_second_order;

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

        extern std::vector<std::vector<std::vector<double> > > second_order_tensor;
        extern std::vector<std::vector<std::vector<std::vector<double> > > > third_order_tensor;
        extern std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > fourth_order_tensor;

        //-------------------------------------------------------------------------
        // internal function declarations
        //-------------------------------------------------------------------------

    } // end of internal namespace

} // end of anisotropy namespace

#endif //ANISOTROPY_INTERNAL_H_
