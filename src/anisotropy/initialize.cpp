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
#include "errors.hpp"
#include "vio.hpp"

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
            /* check if anisotropy calculation enabled. If not, do nothing */
            if (!anisotropy::enabled) return;

            /* output informative message */
            zlog << zTs() << "Initialising data structures for anisotropy calculation." << std::endl;

            /* check for prior initialisation */
            if (anisotropy::internal::initialised)
            {
               zlog << zTs() << "Warning: Anisotropy calculation already initialised. Continuing." << std::endl;
               return;
            }

            internal::num_atoms = num_atoms;
            // anisotropy::internal::atom_type_array = atom_type_array;
            // anisotropy::internal::materialscalaranisotropyarray = materialscalaranisotropyarray;

            // anisotropy::internal::spin_array_x = spin_array_x;
            // anisotropy::internal::spin_array_y = spin_array_y;
            // anisotropy::internal::spin_array_z = spin_array_z;

            // anisotropy::internal::field_array.resize(num_atoms);

            // anisotropy::calculate_anisotropy_fields(num_atoms);
        }

    } // end of anisotropy namespace
