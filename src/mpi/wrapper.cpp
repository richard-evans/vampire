//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "vmpi.hpp"

// Internal vmpi header

namespace vmpi{

//------------------------------------------------------------------------------
// Wrapper function for MPI barrier
//------------------------------------------------------------------------------
void barrier(){

   // Wait for all processors just in case anyone else times out
   #ifdef MPICF
      MPI_Barrier(MPI_COMM_WORLD);
   #endif

   return;

}

}
