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
    // function to calculate anisotropy fields
    //----------------------------------------------------------------------------
    void calculate_anisotropy_fields(
        const int num_atoms,
        const std::vector<int>& atom_type_array,
        const std::vector<zkval_t>& materialscalaranisotropyarray,
        const std::vector<double>& atom_coords_x,
        const std::vector<double>& atom_coords_y,
        const std::vector<double>& atom_coords_z,
        const std::vector<double>& spin_array_x,
        const std::vector<double>& spin_array_y,
        const std::vector<double>& spin_array_z)
        {
            if (anisotropy::internal::uniaxial) calculate_uniaxial_fields(num_atoms);
            // if (anisotropy::internal::2ndorder) calculate_2ndorder_fields(num_atoms);
        }

        double calculate_uniaxial_fields(num_atoms)
        {
            for (int atom=0; atom<num_atoms; ++atom)
            {
                int mat = atoms_type_array[atom];
                double uniaxial_field = 2.0 * materialscalaranisotropyarray[mat]
                                            * spin_array_z[atom];

                anisotropy::internal::field_array[atom] += uniaxial_field;
            }
        }

        return;

    } // end of anisotropy namespace
