//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "dipole.hpp"
#include "sim.hpp"
#include "spintorque.hpp"
#include "spintransport.hpp"
#include "unitcell.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// sim namespace
namespace sim{

// sim::internal namespace
namespace internal{

//------------------------------------------------------------------------------
// Function to increment time counter and associated variables
//------------------------------------------------------------------------------
void increment_time(){

   // set flag checkpoint_loaded_flag to false since first step of simulations was performed
   sim::checkpoint_loaded_flag=false;

	sim::time++;
	// sim::head_position[0]+=sim::head_speed*mp::dt_SI*1.0e10;

   // Update dipole fields
   dipole::calculate_field(sim::time, atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, atoms::m_spin_array, atoms::magnetic);

	if(sim::lagrange_multiplier) update_lagrange_lambda();
   st::update_spin_torque_fields(atoms::x_spin_array,
                               atoms::y_spin_array,
                               atoms::z_spin_array,
                               atoms::type_array,
                               mp::mu_s_array);

   //---------------------------------------------------------------------
   // update spin transport solver resistance and current and STT fields
   //---------------------------------------------------------------------
   // Determine number of local atoms
   #ifdef MPICF
      const unsigned int num_local_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
   #else
      const unsigned int num_local_atoms = atoms::num_atoms;
   #endif

   spin_transport::update(num_local_atoms,
                          atoms::x_spin_array,
                          atoms::y_spin_array,
                          atoms::z_spin_array,
                          atoms::m_spin_array);

   return;

}

} // end of internal namespace
} // end of sim namespace
