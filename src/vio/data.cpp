//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Rory Pond 2016. All rights reserved.
//
//   Email: rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
// Headers
#include "vio.hpp"

// vio module headers
#include "internal.hpp"

// Global output filestreams
std::ofstream zinfo;
std::ofstream zlog;
std::ofstream zmag;
std::ofstream zgrain;

namespace vout{

	// Namespace variable declarations
	std::vector<unsigned int> file_output_list(0);
	std::vector<unsigned int> screen_output_list(0);
	std::vector<unsigned int> grain_output_list(0);

	// Variables to control rate of data output to screen, output file and grain file
	int output_rate=1;
	int output_grain_rate=1;
	//int output_screen_rate=1; needs to be implemented

	bool gnuplot_array_format=false;

}
