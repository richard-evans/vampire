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

// Vampire Header files
#include "errors.hpp"
#include "vmpi.hpp"

namespace err
{
	bool check=false;

	void vexit(){
		
		// check calling of routine if error checking is activated
		if(err::check==true){std::cout << "err::vexit has been called" << std::endl;}

		// Abort MPI processes for parallel execution
		#ifdef MPICF
		std::cerr << "Aborting MPI Processes initiated on rank: " << vmpi::my_rank << " ...";
		MPI::COMM_WORLD.Abort(EXIT_FAILURE);
		std::cerr << "Done." << std::endl;
		#endif
		
		// Now exit program disgracefully
		exit(EXIT_FAILURE);
}

} // end of namespace err

