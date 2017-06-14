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
#include <vector>

// Vampire headers
#include "anisotropy.hpp"
#include "atoms.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "material.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

    //----------------------------------------------------------------------------
    // function to initialize anisotropy module
    //----------------------------------------------------------------------------
    void initialise (
        const int num_atoms,
        std::vector<int>& atom_type_array,
        std::vector<zkval_t>& materialscalaranisotropyarray,
        std::vector<double>& atom_coords_x,
        std::vector<double>& atom_coords_y,
        std::vector<double>& atom_coords_z,
        std::vector<double>& spin_array_x,
        std::vector<double>& spin_array_y,
        std::vector<double>& spin_array_z)
        {
            /* output informative message */
            zlog << zTs() << "Initialising data structures for anisotropy calculation." << std::endl;

            /* check for prior initialisation */
            if (internal::initialised)
            {
               zlog << zTs() << "Warning: Anisotropy calculation already initialised. Continuing." << std::endl;
               return;
            }

            /* internalise variables */
            internal::atom_type_array = atom_type_array;
            internal::materialscalaranisotropyarray = materialscalaranisotropyarray;

            internal::field_array.resize(atoms::num_atoms);

            /*
             * initialise tensors
             */

            if (uniaxial)
            {


                /* resize tensors */
                internal::second_order_tensor.resize(9);

            }

            internal::initialised = true;

            return;
        }

    } // end of anisotropy namespace
