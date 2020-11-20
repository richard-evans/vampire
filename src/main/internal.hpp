//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef VMAIN_INTERNAL_H_
#define VMAIN_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the main module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <cstdint>
#include <string>

namespace vmain{

namespace internal{

   //-------------------------------------------------------------------------
   // Internal shared variables
   //-------------------------------------------------------------------------
   extern std::string input_file_name;

   //-------------------------------------------------------------------------
   // Internal function declarations
   //-------------------------------------------------------------------------
   void command_line_args(int argc, char* argv[]); // process command line arguments


} // end of internal namespace

} // end of main namespace

#endif //VMAIN_INTERNAL_H_
