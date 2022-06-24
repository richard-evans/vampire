//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Standard Headers
#include <iostream>
#include <fstream>

// Vampire headers
#include "errors.hpp"
#include "atoms.hpp"
#include "cells.hpp"
#include "create.hpp"
#include "dipole.hpp"
#include "grains.hpp"
#include "ltmp.hpp"
#include "material.hpp"
#include "neighbours.hpp"
#include "sim.hpp"
#include "spintorque.hpp"
#include "unitcell.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

// Internal create header
#include "internal.hpp"

namespace cs{

	// System Dimensions
	double system_dimensions[3]={77.0,77.0,77.0};	/// Size of system (A)
	bool pbc[3]={false,false,false};						/// Periodic boundary conditions

	unsigned int total_num_unit_cells[3]={0,0,0};	/// Unit cells for entire system (x,y,z)
	unsigned int local_num_unit_cells[3]={0,0,0};	/// Unit cells on local processor (x,y,z)

	// System Parameters
	int particle_creation_parity=0; /// Offset of particle centre (odd/even)
	double particle_scale=50.0;     /// Diameter of particles/grains (A)
	double particle_spacing=10.0;   /// Spacing Between particles (A)
	double particle_array_offset_x=0.0; /// Offset particle array along x-direction;
	double particle_array_offset_y=0.0; /// Offset particle array along y-direction;
   double particle_shape_factor_x=1.0; /// Normalised particle shape
   double particle_shape_factor_y=1.0; /// Normalised particle shape
   double particle_shape_factor_z=1.0; /// Normalised particle shape

	// Other directives and flags
	bool single_spin=false;
	int system_creation_flags[10]={0,0,0,0,0,0,0,0,0,0};

	bool fill_core_shell=true;
   bool core_shell_particles = false;

   // Variables for multilayer system
   bool multilayers = false;
   bool multilayer_height_category = false; // enable height categorization by multilayer number
   int num_multilayers = 1;

	// Variables for interfacial roughness control
	bool interfacial_roughness=false;
	bool interfacial_roughness_local_height_field=false;
	int interfacial_roughness_type=0; /// Sets peaks (1), troughs (-1) or both (0)
	unsigned int interfacial_roughness_random_seed=23456;
	unsigned int interfacial_roughness_seed_count=20; /// Number of seeds
	double interfacial_roughness_height_field_resolution=3.5; /// Angstroms
	double interfacial_roughness_mean_seed_radius=30.0; /// Angstroms
	double interfacial_roughness_seed_radius_variance=0.5; /// Variance as fraction of mean radius
	double interfacial_roughness_mean_seed_height=3.0; /// Angstroms
	double interfacial_roughness_seed_height_max=1.8; /// Angstroms

	// unit cell container
	unitcell::unit_cell_t unit_cell;

  // Array for storing non-magnetic atoms
  std::vector<nm_atom_t> non_magnetic_atoms_array;

int create(){
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::create has been called" << std::endl;}

	if(vmpi::my_rank==0){
		std::cout << "Creating system" << std::endl;
	}
	//=============================================================
	//      System creation variables
	//=============================================================

   // initialise create module parameters
   create::initialize();

	// Atom creation array
	std::vector<cs::catom_t> catom_array;

	// initialise unit cell for system
	uc::initialise(cs::unit_cell);

   // Instantiate some constants for improved readability
   const double ucx = unit_cell.dimensions[0];
   const double ucy = unit_cell.dimensions[1];
   const double ucz = unit_cell.dimensions[2];
   const unsigned int na = unit_cell.atom.size();

   // Calculate number of global and local unit cells required (rounding up)
   // Must be set before rounding up system dimensions for periodic boundary conditions
   cs::total_num_unit_cells[0]=int(vmath::iceil(cs::system_dimensions[0]/unit_cell.dimensions[0]));
   cs::total_num_unit_cells[1]=int(vmath::iceil(cs::system_dimensions[1]/unit_cell.dimensions[1]));
   cs::total_num_unit_cells[2]=int(vmath::iceil(cs::system_dimensions[2]/unit_cell.dimensions[2]));

	// check for pbc and if so round up system dimensions
	if(cs::pbc[0]==true) cs::system_dimensions[0]=unit_cell.dimensions[0]*(int(vmath::iceil(cs::system_dimensions[0]/unit_cell.dimensions[0])));
	if(cs::pbc[1]==true) cs::system_dimensions[1]=unit_cell.dimensions[1]*(int(vmath::iceil(cs::system_dimensions[1]/unit_cell.dimensions[1])));
	if(cs::pbc[2]==true) cs::system_dimensions[2]=unit_cell.dimensions[2]*(int(vmath::iceil(cs::system_dimensions[2]/unit_cell.dimensions[2])));

	// Set up Parallel Decomposition if required
	#ifdef MPICF
		if(vmpi::mpi_mode==0) vmpi::geometric_decomposition(vmpi::num_processors,cs::system_dimensions);
	#endif

	// Create block of crystal of desired size
	cs::create_crystal_structure(catom_array);

	// Cut system to the correct type, species etc
	create::create_system_type(catom_array);

	// Copy atoms for interprocessor communications
	#ifdef MPICF
	if(vmpi::mpi_mode==0){
		create::internal::copy_halo_atoms(catom_array);
   }
	#endif

   //---------------------------------------------
	// Create Neighbour lists for system
   //---------------------------------------------
   neighbours::list_t bilinear; // bilinear exchange list
   neighbours::list_t biquadratic; // biquadratic exchange list

   // generate bilinear exchange list
   bilinear.generate(catom_array, cs::unit_cell.bilinear, na, ucx, ucy, ucz);

   // optionally create a biquadratic neighbour list
   if(exchange::biquadratic){
      biquadratic.generate(catom_array, cs::unit_cell.biquadratic, na, ucx, ucy, ucz);
   }

	#ifdef MPICF
		create::internal::identify_mpi_boundary_atoms(catom_array,bilinear);
      if(exchange::biquadratic) create::internal::identify_mpi_boundary_atoms(catom_array,biquadratic);
      create::internal::mark_non_interacting_halo(catom_array);
      // Sort Arrays by MPI Type
      create::internal::sort_atoms_by_mpi_type(catom_array, bilinear, biquadratic);
	#endif

	#ifdef MPICF
      // ** Must be done in parallel **
		create::internal::init_mpi_comms(catom_array);
		vmpi::barrier();
	#endif

	// Print informative message
	std::cout << "Copying system data to optimised data structures." << std::endl;
	zlog << zTs() << "Copying system data to optimised data structures." << std::endl;

	create::internal::set_atom_vars(catom_array, bilinear, biquadratic);

   // Determine number of local atoms
   #ifdef MPICF
   #else
      // set number of core atoms for serial code (to allow wraper functions to work seamlessly)
      vmpi::num_core_atoms = atoms::num_atoms;
		// set the number of local atoms on process (all atoms in serial)
	   vmpi::num_local_atoms = atoms::num_atoms;
   #endif

	// Set grain and cell variables for simulation
	grains::set_properties();

	#ifdef MPICF
		//std::cout << "Outputting coordinate data" << std::endl;
		//vmpi::crystal_xyz(catom_array);
	int my_num_atoms=vmpi::num_core_atoms+vmpi::num_bdry_atoms;
   //std::cout << "my_num_atoms == " << my_num_atoms << std::endl;
	int total_num_atoms=0;
	MPI_Reduce(&my_num_atoms,&total_num_atoms, 1,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	std::cout << "Total number of atoms (all CPUs): " << total_num_atoms << std::endl;
   zlog << zTs() << "Total number of atoms (all CPUs): " << total_num_atoms << std::endl;
	#else
	std::cout << "Number of atoms generated: " << atoms::num_atoms << std::endl;
   zlog << zTs() << "Number of atoms generated: " << atoms::num_atoms << std::endl;

	#endif

	return EXIT_SUCCESS;
}

} // end of create namespace
