#ifndef GPU_H_
#define GPU_H_
//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2015, 2016. All rights reserved.
//
//------------------------------------------------------------------------------
#include <string>

namespace gpu{

   //-----------------------------------------------------------------------------
   // Variables used for GPU acceleration
   //-----------------------------------------------------------------------------
   extern bool acceleration; // flag to enable gpu_acceleration
   extern bool cpu_stats; // flag to calculate stats using cpu
   extern bool initialized; // flag to check for successful initialization

   extern int device; // int specifying gpu device to use for simulation


   //-----------------------------------------------------------------------------
   // Functions for GPU acceleration
   //-----------------------------------------------------------------------------
   extern bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);
   extern void initialize();
   extern void initialize_dipole();
   extern void llg_heun();
   extern void finalize();
   extern void transfer_spin_positions_from_gpu_to_cpu();
   extern void transfer_dipole_fields_from_cpu_to_gpu();
   extern void transfer_dipole_cells_fields_from_gpu_to_cpu();
   extern void update_dipolar_fields();
   extern void mc_step();

   //-----------------------------------------------------------------------------
   // Functions for GPU configuration output
   //-----------------------------------------------------------------------------
   namespace config{
      extern void synchronise();
   }

   //-----------------------------------------------------------------------------
   // Functions for GPU statistics calculation
   //-----------------------------------------------------------------------------
   namespace stats{

      extern void update();
      extern void reset();
      extern void get();

   }

} // end of gpu namespace

#endif //GPU_H_
