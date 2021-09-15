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

// Vampire headers
#include "dipole.hpp"
#include "vio.hpp"

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

} // end of namespace internal
} // end of namespace dipole
