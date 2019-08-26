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

// Vampire headers

//--------------------------------------------------------------------------------
// Namespace for variables and functions for vopencl module
//--------------------------------------------------------------------------------
namespace vopencl
{

   //-----------------------------------------------------------------------------
   // Function for OpenCL acceleration
   //-----------------------------------------------------------------------------
   bool initialize(bool cpu_stats);
   void llg_heun(void);
   void stats_update(void);
   void finalize(void);

   namespace config
   {
      void synchronise(void);
   }

   namespace stats
   {
      void update(void);
      void get(void);
      void reset(void);
   }

} // end of vopencl namespace

#endif //VOPENCL_H_
