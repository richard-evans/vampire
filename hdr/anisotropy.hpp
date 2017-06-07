//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland 2016. All rights reserved.
//
//   Email: sw766@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef ANISOTROPY_H_
#define ANISOTROPY_H_

// C++ standard library headers
//#include <vector>
//#include <string>

// Vampire headers
//#include "anisotropy.hpp"
#include "material.hpp"

//--------------------------------------------------------------------------------
// namespace for variables and functions for anisotropy module
//--------------------------------------------------------------------------------
namespace anisotropy
{
    //-----------------------------------------------------------------------------
    // function to initialise anisotropy module
    //-----------------------------------------------------------------------------
    void initialise(
        const int num_atoms,
        std::vector<int>& atom_type_array,
        std::vector<zkval_t>& materialscalaranisotropyarray,
        std::vector<double>& atom_coords_x,
        std::vector<double>& atom_coords_y,
        std::vector<double>& atom_coords_z,
        std::vector<double>& spin_array_x,
        std::vector<double>& spin_array_y,
        std::vector<double>& spin_array_z);

    //---------------------------------------------------------------------------
    // Function to process input file parameters for anisotropy module
    //---------------------------------------------------------------------------
    bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

    //---------------------------------------------------------------------------
    // Function to process material parameters
    //---------------------------------------------------------------------------
    bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of anisotropy namespace

#endif //ANISOTROPY_H_
