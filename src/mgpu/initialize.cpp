//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) S Jenkins 2026. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>
#include <string>
#include <vector>

// Vampire headers
#include "create.hpp"
#include "errors.hpp"
#include "gpu.hpp"
#include "mgpu.hpp"
#include "program.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// mgpu module headers
#include "internal.hpp"

namespace mgpu{

namespace internal{

   //---------------------------------------------------------------------------
   // Function to stop with an error message for unsupported simulations
   //---------------------------------------------------------------------------
   void unsupported(const std::string& feature){
      err::fatal("Error: " + feature + " is not supported by the multi-GPU version (vampire-cuda-multigpu)");
   }

} // end of internal namespace

//------------------------------------------------------------------------------
// Function to select the GPU of this rank and to check that the simulation is
// supported. Ranks on the same node use the node's GPUs in turn.
//------------------------------------------------------------------------------
void initialize(){

   if(err::check==true) std::cout << "mgpu::initialize has been called" << std::endl;

   // simulation features not yet distributed over several GPUs (those the GPU
   // version does not support at all are checked in gpu::initialize)
   if(vmpi::mpi_mode != 0) internal::unsupported("MPI modes other than geometric decomposition");

   // every local atom needs its own key (create::set_atom_keys)
   if(!create::atom_keys_valid){
      err::fatal("Error: atoms on rank " + std::to_string(vmpi::my_rank) + " cannot be identified by their unit cell (at most " +
                 std::to_string(1 << 16) + " unit cells along each axis and atoms per unit cell), which the multi-GPU version requires");
   }

   // on several GPUs the host has no up to date spins while the simulation runs, so only
   // programs that use just statistics of the spins are supported
   if(vmpi::num_processors > 1 && gpu::program_reads_host_spins()){
      internal::unsupported("Program " + std::to_string(static_cast<int>(program::program)) + " (sim:program) on more than one GPU");
   }

   // rank within this node
   MPI_Comm node;
   MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node);
   int node_rank, node_size;
   MPI_Comm_rank(node, &node_rank);
   MPI_Comm_size(node, &node_size);

   const int num_devices = internal::device_count();
   if(num_devices < 1){
      err::fatal("Error: no CUDA device found on rank " + std::to_string(vmpi::my_rank));
   }

   // the CUDA module sets the device chosen here (unless set in the input file)
   if(gpu::device == -1) gpu::device = node_rank % num_devices;

   // ranks may see different devices (e.g. one each), so compare PCI addresses
   const int id_size = 32;
   std::vector<char> ids(id_size * node_size, 0);
   internal::device_pci_id(gpu::device, &ids[id_size * node_rank], id_size);
   MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ids.data(), id_size, MPI_CHAR, node);
   MPI_Comm_free(&node);
   const std::string my_id(&ids[id_size * node_rank]);
   int sharing = 0;
   for(int r = 0; r < node_size; r++) if(my_id == std::string(&ids[id_size * r])) sharing++;

   zlog << zTs() << "Rank " << vmpi::my_rank << " is rank " << node_rank << " of " << node_size << " on this node and uses GPU "
        << gpu::device << " of " << num_devices << " (PCI " << my_id << ")" << std::endl;
   if(sharing > 1){
      if(vmpi::my_rank == 0) std::cout << "Warning: rank 0 shares its GPU with " << sharing - 1 << " other rank(s)" << std::endl;
      zlog << zTs() << "Warning: this rank shares its GPU with " << sharing - 1 << " other rank(s)" << std::endl;
   }

}

} // end of mgpu namespace
