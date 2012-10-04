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
/// @brief Contains Vampire error functions and variables

/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    10/03/2011
/// @internal
///	Created:		10/03/2011
///	Revision:	  ---
///=====================================================================================
///

// Standard Libraries
#include <cstdlib>
#include <iostream>

// zspin header files
#include "errors.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

namespace err
{
	bool check=false;

	// Default vexit
	void vexit(){
		
		// check calling of routine if error checking is activated
		if(err::check==true){std::cout << "err::vexit has been called" << std::endl;}

		// Abort MPI processes for parallel execution
		#ifdef MPICF
		std::cerr << "Fatal error on rank: " << vmpi::my_rank << ": Aborting program." << std::endl;
		zlog << zTs() << "Fatal error on rank: " << vmpi::my_rank << ": Aborting program." << std::endl;
		// concatenate log and sort
		system("ls zlog.* | xargs cat | sort -n > zlog");
		MPI::COMM_WORLD.Abort(EXIT_FAILURE);
		// MPI program dies ungracefully here
		#endif
		
		// Print error message to screen and log
		zlog << zTs() << "Fatal error: Aborting program." << std::endl;
		std::cout << "Fatal error: Aborting program. See zlog file for details." << std::endl;

		// Now exit program disgracefully
		exit(EXIT_FAILURE);
	}

	// zexit with error message
	void zexit(std::string message){
		
		// check calling of routine if error checking is activated
		if(err::check==true){std::cout << "err::zexit(std::string) has been called" << std::endl;}
		
		// Abort MPI processes for parallel execution
		#ifdef MPICF
		std::cerr << "Fatal error on rank " << vmpi::my_rank << ": " << message << std::endl;
		std::cerr << "Aborting program." << std::endl;
		zlog << zTs() << "Fatal error on rank " << vmpi::my_rank << ": " << message << std::endl;
		zlog << zTs() << "Aborting program." << std::endl;

		// concatenate log and sort
		system("ls zlog.* | xargs cat | sort -n > zlog");
		MPI::COMM_WORLD.Abort(EXIT_FAILURE);
		// MPI program dies ungracefully here
		#else
			// Print Error message to screen and log
			std::cerr << "Fatal error: " << message << std::endl;
			zlog << zTs() << "Fatal error: " << message << std::endl;

		#endif
		
		// Print error message to screen and log
		zlog << zTs() << "Aborting program." << std::endl;
		std::cout << "Aborting program. See zlog file for details." << std::endl;

		// Now exit program disgracefully
		exit(EXIT_FAILURE);
	}
	
} // end of namespace err

