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

            /* resize tensors */
            internal::second_order_tensor.resize(atoms::num_atoms);
            for (int i=0; i<atoms::num_atoms; i++)
                internal::second_order_tensor.at(i).resize(9);

            /* initialise all elements to zero */
            for (int i=0; i<atoms::num_atoms; i++)
                for (int j=0; j<9; j++)
                    internal::second_order_tensor.at(i).at(j) = 0;

            if (uniaxial)
            {
                for (int atom=0; atom<atoms::num_atoms; atom++)
                {
                    int mat = atoms::type_array.at(atom);

                    double ex = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(0);
                    double ey = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(1);
                    double ez = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(2);

                    internal::second_order_tensor.at(atom).at(0) += mp::material.at(mat).Ku * ex * ex;
                    internal::second_order_tensor.at(atom).at(1) += mp::material.at(mat).Ku * ex * ey;
                    internal::second_order_tensor.at(atom).at(2) += mp::material.at(mat).Ku * ex * ez;

                    internal::second_order_tensor.at(atom).at(3) += mp::material.at(mat).Ku * ey * ex;
                    internal::second_order_tensor.at(atom).at(4) += mp::material.at(mat).Ku * ey * ey;
                    internal::second_order_tensor.at(atom).at(5) += mp::material.at(mat).Ku * ey * ez;

                    internal::second_order_tensor.at(atom).at(6) += mp::material.at(mat).Ku * ez * ex;
                    internal::second_order_tensor.at(atom).at(7) += mp::material.at(mat).Ku * ez * ey;
                    internal::second_order_tensor.at(atom).at(8) += mp::material.at(mat).Ku * ez * ez;
                }
            }

            if (neel)
            {
                for (int atom = 0; atom < atoms::num_atoms; atom++)
                {
                    int mat = atoms::type_array.at(atom);
                    int start = atoms::nearest_neighbour_list_si.at(atom);
                    int end = atoms::nearest_neighbour_list_ei.at(atom);

                    double Ks = 0.5 * 2.0 * mp::material.at(mat).Ks; // note factor two from differentiation

                    for (int n = start; n < end; n++)
                    {
                        double eijx = atoms::eijx.at(n);
                        double eijy = atoms::eijy.at(n);
                        double eijz = atoms::eijz.at(n);

                        internal::second_order_tensor.at(atom).at(0) += eijx * eijx * Ks;
                        internal::second_order_tensor.at(atom).at(1) += eijx * eijy * Ks;
                        internal::second_order_tensor.at(atom).at(2) += eijx * eijz * Ks;

                        internal::second_order_tensor.at(atom).at(3) += eijy * eijx * Ks;
                        internal::second_order_tensor.at(atom).at(4) += eijy * eijy * Ks;
                        internal::second_order_tensor.at(atom).at(5) += eijy * eijz * Ks;

                        internal::second_order_tensor.at(atom).at(6) += eijz * eijx * Ks;
                        internal::second_order_tensor.at(atom).at(7) += eijz * eijy * Ks;
                        internal::second_order_tensor.at(atom).at(8) += eijz * eijz * Ks;
                    }

                }
            }

            internal::initialised = true;

            return;
        }

    } // end of anisotropy namespace
