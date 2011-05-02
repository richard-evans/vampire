///
///	@file
///	@brief Category Class data structure
///
///	Categories identify a collection of atoms in a way similar to, but separate from,  their material. 
///	They support output of categorised magnetisation, eg in a plane of equal height, or radial dependent magneisation. 
///	Also allows constraint of a subset of atoms using CMC/Lagrange. Can be used to track sublattices in a single material.
///
///	@author Richard Evans, richard.evans@york.ac.uk
///
///	@section License
///	Use of this code, either in source or compiled form, is subject to license from the authors.
///	Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2011. All Rights Reserved.
///
///	 @internal
///	Created  30/04/2011
///	Revision  1.0
///	Copyright  Copyright (c) 2011, Richard Evans
///
///=====================================================================================
///

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
		phi(0.0)
	{
		category_atoms.resize(0,0);
		
	}
		
		
} // End of namespace grains
