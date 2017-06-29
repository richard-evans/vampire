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

    namespace internal{

        //------------------------------------------------------------------------
        // Shared variables inside anisotropy module
        //------------------------------------------------------------------------

        bool initialised = false;

        bool uniaxial_first_order = false;
        bool neel = false;

        bool uniaxial_second_order = false;

        std::vector<double> field_array;

        int num_atoms;
        std::vector<int> atom_type_array;
        std::vector<zkval_t> materialscalaranisotropyarray;

        std::vector<double> spin_array_x;
        std::vector<double> spin_array_y;
        std::vector<double> spin_array_z;

        std::vector<std::vector<std::vector<double> > > second_order_tensor;
        std::vector<std::vector<std::vector<std::vector<double> > > > third_order_tensor;
        std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > fourth_order_tensor;

    } // end of internal namespace

} // end of anisotropy namespace
