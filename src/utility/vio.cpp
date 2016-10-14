//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains vin and vout namespaces for file input.output in vampire.
///
/// @details File and screen input and output are controlled via the separate namespaces.
///
/// @section notes Implementation Notes
/// This is a list of other notes, not related to functionality but rather to implementation.
/// Also include references, formulae and other notes here.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    15/01/2010
/// @internal
///	Created:		15/01/2010
///	Revision:	  ---
///=====================================================================================
///

// Headers
#include "atoms.hpp"
#include "cells.hpp"
#include "create.hpp"
#include "config.hpp"
#include "demag.hpp"
#include "errors.hpp"
#include "grains.hpp"
#include "ltmp.hpp"
#include "voronoi.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "units.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "../vio/internal.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#ifdef WIN_COMPILE
#include <direct.h>
#endif


// Unused functions 
//
//namespace vout{
//
//	#ifdef MPICF
//	null_streambuf nullbuf;
//
//	void redirect(std::ostream& strm, std::string filename) {
//		errfile.open(filename.c_str());
//		// redirect ouput into the file
//		strm.rdbuf (errfile.rdbuf());
//	}
//
//	void nullify(std::ostream& strm){
//		strm.rdbuf(&nullbuf);
//	}
//	#endif
//} // end of namespace vout
//