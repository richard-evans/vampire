//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans and Joe Barker 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// vio module headers
#include "internal.hpp"

// C++ headers
#include <iostream>
#include <fstream>

#ifdef MPICF

struct null_streambuf
: public std::streambuf
{
  void overflow(char c)
  {
  }
} nullbuf;

namespace vout{

   // global scope error file for redirected stream
   std::ofstream errfile;

   //------------------------------------------------------------------------------
   // Function to redirect output streams to file
   //------------------------------------------------------------------------------
	void redirect(std::ostream& strm, std::string filename) {
		errfile.open(filename.c_str());
		// redirect ouput into the file
		strm.rdbuf (vout::errfile.rdbuf());
	}

   //------------------------------------------------------------------------------
   // Function to nullify output streams for parallel version
   //------------------------------------------------------------------------------
	void nullify(std::ostream& strm){
		strm.rdbuf(&nullbuf);
	}

} // end of namespace vout

#endif
