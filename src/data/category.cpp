//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and Mara strungaru 2023. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Vampire header files
#include "atoms.hpp"
#include "category.hpp"
#include "grains.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

// System header files
#include <cmath>
#include <iostream>

namespace cat{

	// vector of categories for book keeping
	std::vector <category_t> category(0);

	// constructor
	category_t::category_t():
		magm(0.0),
		mean_magm(0.0),
		mx(0.0),
		my(0.0),
		mz(0.0),
		torque(0.0),
		mean_torque(0.0),
		energy(0.0),
		theta(0.0),
		phi(0.0),
		spin_temp(0.0),
		mean_spin_temp(0.0)
	{
		category_atoms.resize(0,0);
	}

} // End of namespace cat
