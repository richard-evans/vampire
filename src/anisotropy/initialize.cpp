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

    //----------------------------------------------------------------------------
    // function to initialize anisotropy module
    //----------------------------------------------------------------------------
    void initialise (
        const int num_atoms,
        const std::vector<int>& atom_type_array,
        const std::vector<double>& materialscalaranisotropyarray,
        const std::vector<double>& atom_coords_x,
        const std::vector<double>& atom_coords_y,
        const std::vector<double>& atom_coords_z,
        const std::vector<double>& spin_array_x,
        const std::vector<double>& spin_array_y,
        const std::vector<double>& spin_array_z)
        {
            /* check if anisotropy calculation enabled. If not, do nothing */
            if (!anisotropy::internal::enabled) return;

            /* output informative message */
            zlog << zTs() << "Initialising data structures for anisotropy calculation." << std::endl;

            /* check for prior initialisation */
            if (anisotropy::internal::initialised)
            {
               zlog << zTs() << "Warning: Anisotropy calculation already initialised. Continuing." << std::endl;
               return;
            }

            anisotropy::internal::num_atoms = num_atoms;
            anisotropy::internal::atom_type_array = atom_type_array;
            anisotropy::internal::materialscalaranisotropyarray = materialscalaranisotropyarray;

            anisotropy::internal::spin_array_x = spin_array_x;
            anisotropy::internal::spin_array_y = spin_array_y;
            anisotropy::internal::spin_array_z = spin_array_z;

            return;
        }

    } // end of anisotropy namespace
