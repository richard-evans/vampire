//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2019. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <sstream>

// Vampire headers
#include "dipole.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"

// dipole module headers
#include "internal.hpp"

namespace dipole{
namespace internal{


//-----------------------------------------------------------------
// Function to output solver memory information
//-----------------------------------------------------------------
void output_dipole_solver_mem_info(int num_cells, int num_local_cells){

   // Check memory requirements and print to screen (6 tensor components and 8 bytes per number)
   zlog << zTs() << "\tDipole field calculation requires " << double(num_cells)*double(num_local_cells * 6.0) * 8.0 / 1.0e6 << " MB of RAM" << std::endl;
   std::cout     << "Dipole field calculation requires "   << double(num_cells)*double(num_local_cells * 6.0) * 8.0 / 1.0e6 << " MB of RAM" << std::endl;

   // output parallel information for MPI code ( num cells ^ 2 )
   #ifdef MPICF
      zlog << zTs() << "\tTotal memory for dipole calculation (all CPUs): " << double(num_cells)*double(num_cells * 6.0) * 8.0 / 1.0e6 << " MB of RAM" << std::endl;
      std::cout     << "Total memory for dipole calculation (all CPUs): "   << double(num_cells)*double(num_cells * 6.0) * 8.0 / 1.0e6 << " MB of RAM" << std::endl;
      zlog << zTs() << "\tNumber of local cells for dipole calculation = " << num_local_cells << std::endl;
      zlog << zTs() << "\tNumber of total cells for dipole calculation = " << num_cells << std::endl;
   #endif

   return;

}

//-----------------------------------------------------------------
// Function to stop with a message if the stored tensors of the tensor
// solver need more memory than the node has available
//-----------------------------------------------------------------
void check_tensor_solver_memory(int num_cells, int num_local_cells){

   // tensors of all processes on this node (6 components, 8 bytes each), in kB
   double required_kb = double(num_cells) * double(num_local_cells) * 6.0 * 8.0 / 1024.0;
   #ifdef MPICF
      MPI_Comm node_comm;
      MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, vmpi::my_rank, MPI_INFO_NULL, &node_comm);
      MPI_Allreduce(MPI_IN_PLACE, &required_kb, 1, MPI_DOUBLE, MPI_SUM, node_comm);
      MPI_Comm_free(&node_comm);
   #endif

   const double available_kb = vutil::available_memory_kb();
   int too_large = (available_kb >= 0.0 && required_kb > available_kb);
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &too_large, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   #endif
   if(!too_large) return;

   std::ostringstream message;
   message << "Error - the dipole tensors need " << required_kb * 1024.0 / 1.0e9 << " GB of memory";
   if(available_kb >= 0.0 && required_kb > available_kb) message << " on this node, but only " << available_kb * 1024.0 / 1.0e9 << " GB are available";
   else message << " on one of the nodes, more than is available there";
   message << ". Set dipole:on-the-fly-tensors to compute the tensors at each field update instead of storing them"
           << " (the same fields with the memory of one row of tensors, but each field update then takes about as long"
           << " as the tensor set-up), or run the parallel version on more nodes.";
   err::fatal(message.str());

}

} // end of namespace internal
} // end of namespace dipole
