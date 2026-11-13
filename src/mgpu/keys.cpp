//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) S Jenkins 2026. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <cstdint>
#include <vector>

// Vampire headers
#include "create.hpp"
#include "mgpu.hpp"

namespace mgpu{

// the keys are stored by the create module (create::set_atom_keys)
std::vector<uint64_t>& atom_key = create::atom_keys;

} // end of mgpu namespace
