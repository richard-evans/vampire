//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) S Jenkins 2026. All rights reserved.
//
//------------------------------------------------------------------------------
#ifndef MGPU_H_
#define MGPU_H_

// C++ standard library headers
#include <cstdint>
#include <vector>

//-----------------------------------------------------------------------------
// Multi-GPU module (make cuda-multigpu): one MPI rank per GPU. The system is
// decomposed geometrically as in the MPI version; each rank integrates its
// core and boundary atoms on its GPU with the kernels of the CUDA module and
// exchanges the spins of halo atoms with its neighbours before every field
// evaluation.
//-----------------------------------------------------------------------------
namespace mgpu{

   // Selects the GPU of this rank and checks that the simulation is supported.
   // Called before the CUDA module is initialised.
   void initialize();

   // Sets up halo exchange on the device and, if requested, orders the exchange
   // neighbours of each atom by atom key. Called after the CUDA module is initialised.
   void initialize_device();

   // Key of each atom (create::atom_keys)
   extern std::vector<uint64_t>& atom_key;

   // Sum the exchange contributions of each atom in atom key order
   // (gpu:decomposition-independent-order, default off). Needed for results that
   // are bit for bit identical on any number of GPUs, e.g. for testing; costs a
   // sort of the neighbour list at start-up.
   extern bool decomposition_independent_order;

   // Performs a single LLG Heun step on the device
   void llg_heun();

   // Frees device and host buffers of the module
   void finalize();

} // end of mgpu namespace

#endif //MGPU_H_
