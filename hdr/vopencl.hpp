//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//------------------------------------------------------------------------------
//

#ifndef VOPENCL_H_
#define VOPENCL_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "vopencl.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for vopencl module
//--------------------------------------------------------------------------------
namespace vopencl
{

   //-----------------------------------------------------------------------------
   // Function for OpenCL acceleration
   //-----------------------------------------------------------------------------
   bool initialize();
   void llg_heun();
   void stats_update();
   void finalize();

   namespace config
   {
      void synchronise();
   }

   namespace stats
   {
      void update();
      void get();
      void reset();
   }

} // end of vopencl namespace

#endif //VOPENCL_H_
