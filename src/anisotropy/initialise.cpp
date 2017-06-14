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

                /* initialise all elements to zero */
                for (int i = 0; i < internal::second_order_tensor.size(); ++i)
                    internal::second_order_tensor.at(i) = 0;

                for (int mat = 0; mat < mp::num_materials; mat++)
                {
                    double ex = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(0);
                    double ey = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(1);
                    double ez = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(2);

                    internal::second_order_tensor.at(0) += mp::material.at(mat).Ku * ex * ex;
                    internal::second_order_tensor.at(1) += mp::material.at(mat).Ku * ex * ey;
                    internal::second_order_tensor.at(2) += mp::material.at(mat).Ku * ex * ez;

                    internal::second_order_tensor.at(3) += mp::material.at(mat).Ku * ey * ex;
                    internal::second_order_tensor.at(4) += mp::material.at(mat).Ku * ey * ey;
                    internal::second_order_tensor.at(5) += mp::material.at(mat).Ku * ey * ez;

                    internal::second_order_tensor.at(6) += mp::material.at(mat).Ku * ez * ex;
                    internal::second_order_tensor.at(7) += mp::material.at(mat).Ku * ez * ey;
                    internal::second_order_tensor.at(8) += mp::material.at(mat).Ku * ez * ez;
                }
            }

            internal::initialised = true;

            return;
        }

    } // end of anisotropy namespace
