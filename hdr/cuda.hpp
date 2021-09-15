#ifndef CUDA_H_
#define CUDA_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

namespace vcuda{

#ifdef CUDA

   //-----------------------------------------------------------------------------
   // Variables used for cuda GPU acceleration
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Functions for cuda GPU acceleration
   //-----------------------------------------------------------------------------
   extern bool initialize(bool cpu_stats);
   extern bool initialize_dipole();
   extern void llg_heun();
   extern void stats_update();
   extern void finalize();
   extern void transfer_spin_positions_from_gpu_to_cpu();
   extern void transfer_dipole_fields_from_cpu_to_gpu();
   extern void transfer_dipole_cells_fields_from_gpu_to_cpu();
   extern void update_dipolar_fields();
   extern void mc_step();

   namespace config{
      extern void synchronise();
   }

   namespace stats
   {
      void update ();
      void get ();
      void reset ();
   } /* stats */

#endif

} // end of vcuda namespace

#endif //CUDA_H_
