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
#include <string>

// Vampire headers
#include "anisotropy.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for anisotropy module
//--------------------------------------------------------------------------------
namespace anisotropy{

    extern bool uniaxial;

    //-----------------------------------------------------------------------------
    // Function to initialise anisotropy module
    //-----------------------------------------------------------------------------
    void initialize();

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
