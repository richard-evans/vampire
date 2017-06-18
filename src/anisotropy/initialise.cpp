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
            for (int atom = 0; atom < atoms::num_atoms; ++atom)
            {
                internal::second_order_tensor.at(atom).resize(3);
                for (int i = 0; i < 3; ++i)
                    internal::second_order_tensor.at(atom).at(i).resize(3);
            }

            internal::third_order_tensor.resize(atoms::num_atoms);
            for (int atom = 0; atom < atoms::num_atoms; ++atom)
            {
                internal::third_order_tensor.at(atom).resize(3);
                for (int i = 0; i < 3; ++i)
                {
                    internal::third_order_tensor.at(atom).at(i).resize(3);
                    for (int j = 0; j < 3; ++j)
                        internal::third_order_tensor.at(atom).at(i).at(j).resize(3);
                }
            }

            /* initialise all elements to zero */
            for (int atom = 0; atom < atoms::num_atoms; ++atom)
            for (int i = 0; i<3; ++i)
            for (int j=0; j<3; j++)
                internal::second_order_tensor.at(atom).at(i).at(j) = 0;

            if (uniaxial_first_order)
            {
                for (int atom=0; atom<atoms::num_atoms; ++atom)
                {
                    int mat = atoms::type_array.at(atom);
                    double Ku = mp::material.at(mat).Ku;

                    double e[3];
                    e[0] = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(0);
                    e[1] = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(1);
                    e[2] = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(2);

                    for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        internal::second_order_tensor.at(atom).at(i).at(j) += e[i] * e[j] * Ku;
                }
            }

            if (neel)
            {
                for (int atom = 0; atom < atoms::num_atoms; atom++)
                {
                    int mat = atoms::type_array.at(atom);
                    int start = atoms::nearest_neighbour_list_si.at(atom);
                    int end = atoms::nearest_neighbour_list_ei.at(atom);

                    // surface constant: note factor 2 from differentiation
                    double Ks = 0.5 * 2.0 * mp::material.at(mat).Ks;

                    for (int n = start; n < end; ++n)
                    {
                        double eij[3];
                        eij[0] = atoms::eijx.at(n);
                        eij[1] = atoms::eijy.at(n);
                        eij[2] = atoms::eijz.at(n);

                        for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                            internal::second_order_tensor.at(atom).at(i).at(j)
                                += eij[i] * eij[j] * Ks;
                    }

                }
            }

            if (uniaxial_second_order)
            {
                for (int atom = 0; atom < atoms::num_atoms; atom ++)
                {
                    int mat = atoms::type_array.at(atom);
                    double Ku2 = 4.0*mp::material_second_order_anisotropy_constant_array.at(mat);

                    double e[3];
                    e[0] = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(0);
                    e[1] = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(1);
                    e[2] = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(2);

                    for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                    for (int k = 0; k < 3; ++k)
                        internal::third_order_tensor.at(atom).at(i).at(j).at(k)
                            += e[i] * e[j] * e[k] * Ku2;
                }
            }

            internal::initialised = true;

            return;
        }

    } // end of anisotropy namespace
