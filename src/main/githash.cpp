//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// magic to turn pre-processor variables to strings
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define GITHASH TOSTRING(GHASH)

// C++ standard library headers
#include <string>

namespace vinfo{

// function to return git commit of unique code version at HEAD of repository
std::string githash(){
	return GITHASH;
}

} // end of vinfo namespace
