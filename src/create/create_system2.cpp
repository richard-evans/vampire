//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
/// @file
/// @brief Contains master function for system creation and cs namespace.
///
/// @details This is the detailed description of the funtion of this file
///
/// @section notes Implementation Notes
/// Creation routines are re-written for double precision corrdinates, neighbourlists and
/// mpi decomposition.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    25/02/2010
/// @internal
///	Created:		25/02/2010
///	Revision:	  ---
///=====================================================================================
///
// Standard Headers
#include <iostream>
#include <fstream>

// Vampire headers
#include "errors.hpp"
#include "atoms.hpp"
#include "cells.hpp"
#include "demag.hpp"
#include "grains.hpp"
#include "ltmp.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"
#include "create.hpp"



/// @namespace ns
/// @brief Create System Namespace - includes variables and functions for system creation.
///
/// @internal
///=====================================================================================
///
namespace cs{

	// System Dimensions
	double system_dimensions[3]={77.0,77.0,77.0};	/// Size of system (A)
	double unit_cell_size[3]={3.54,3.54,3.54};		/// Unit Cell Size (A) [Will eventually be local to unit cells]
	bool pbc[3]={false,false,false};						/// Periodic boundary conditions
	bool SelectMaterialByZHeight=false;					/// Toggle overwriting of material id by z-height
	bool SelectMaterialByGeometry=false;					/// Toggle override of input material type by geometry
	unsigned int total_num_unit_cells[3]={0,0,0};	/// Unit cells for entire system (x,y,z)
	unsigned int local_num_unit_cells[3]={0,0,0};	/// Unit cells on local processor (x,y,z)
	std::string crystal_structure="sc";

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
	std::string unit_cell_file="";
	bool fill_core_shell=true;

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
	cs::unit_cell_t unit_cell;

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

	int num_atoms=0;							///> num_atoms for generation routines
	//int supercell_dim[3];				///> integer dimesions of system supercells
	//double supercell_size[3];				///> real supercell size (Angstroms)

	// Atom creation array
	std::vector<cs::catom_t> catom_array;
	std::vector<std::vector<neighbour_t> > cneighbourlist;

	// initialise unit cell for system
	unit_cell_set(cs::unit_cell);

	cs::unit_cell_size[0]=unit_cell.dimensions[0];
	cs::unit_cell_size[1]=unit_cell.dimensions[1];
	cs::unit_cell_size[2]=unit_cell.dimensions[2];

   // Calculate number of global and local unit cells required (rounding up)
   // Must be set before rounding up system dimensions for periodic boundary conditions
   cs::total_num_unit_cells[0]=int(vmath::iceil(cs::system_dimensions[0]/unit_cell.dimensions[0]));
   cs::total_num_unit_cells[1]=int(vmath::iceil(cs::system_dimensions[1]/unit_cell.dimensions[1]));
   cs::total_num_unit_cells[2]=int(vmath::iceil(cs::system_dimensions[2]/unit_cell.dimensions[2]));

	// check for pbc and if so round up system dimensions
	if(cs::pbc[0]==true) cs::system_dimensions[0]=cs::unit_cell_size[0]*(int(vmath::iceil(cs::system_dimensions[0]/cs::unit_cell_size[0])));
	if(cs::pbc[1]==true) cs::system_dimensions[1]=cs::unit_cell_size[1]*(int(vmath::iceil(cs::system_dimensions[1]/cs::unit_cell_size[1])));
	if(cs::pbc[2]==true) cs::system_dimensions[2]=cs::unit_cell_size[2]*(int(vmath::iceil(cs::system_dimensions[2]/cs::unit_cell_size[2])));

	// Set up Parallel Decomposition if required
	#ifdef MPICF
		if(vmpi::mpi_mode==0) vmpi::geometric_decomposition(vmpi::num_processors,cs::system_dimensions);
	#endif

	//      Initialise variables for system creation
	if(cs::system_creation_flags[0]==1){
		// read_coord_file();
	}

	#ifdef MPICF
	// check for staged replicated data generation
	if(vmpi::replicated_data_staged==true && vmpi::mpi_mode==1){

		// copy ppn to constant temporary
		const int ppn=vmpi::ppn;

		for(int n=0;n<ppn;n++){
			bool i_am_it=false;
			if(vmpi::my_rank%ppn==n) i_am_it=true;

			// Only generate system if I am it
			if(i_am_it){

				zlog << zTs() << "Generating staged system on rank " << vmpi::my_rank << "..." << std::endl;
				//std::cerr << zTs() << "Generating staged system on rank " << vmpi::my_rank << "..." << std::endl;

				// Create block of crystal of desired size
				cs::create_crystal_structure(catom_array);

				// Cut system to the correct type, species etc
				cs::create_system_type(catom_array);

				vmpi::set_replicated_data(catom_array);

				// Create Neighbour list for system
				cs::create_neighbourlist(catom_array,cneighbourlist);

				// Identify needed atoms and destroy the rest
				vmpi::identify_boundary_atoms(catom_array,cneighbourlist);

				zlog << zTs() << "Staged system generation on rank " << vmpi::my_rank << " completed." << std::endl;
				//std::cerr << zTs() << "Staged system generation on rank " << vmpi::my_rank << " completed." << std::endl;
			}

			// Wait for process who is it
			MPI::COMM_WORLD.Barrier();

		} // end of loop over processes
	}
	else{
	#endif

	// Create block of crystal of desired size
	cs::create_crystal_structure(catom_array);

	// Cut system to the correct type, species etc
	cs::create_system_type(catom_array);

	// Copy atoms for interprocessor communications
	#ifdef MPICF
	if(vmpi::mpi_mode==0){
		MPI::COMM_WORLD.Barrier(); // wait for everyone
		vmpi::copy_halo_atoms(catom_array);
		MPI::COMM_WORLD.Barrier(); // sync after halo atoms copied
	}
	else if(vmpi::mpi_mode==1){
		vmpi::set_replicated_data(catom_array);
	}
	#else
		//cs::copy_periodic_boundaries(catom_array);
	#endif

	// Create Neighbour list for system
	cs::create_neighbourlist(catom_array,cneighbourlist);

	#ifdef MPICF
		vmpi::identify_boundary_atoms(catom_array,cneighbourlist);
	#endif


	#ifdef MPICF
	} // stop if for staged generation here
	// ** Must be done in parallel **
		vmpi::init_mpi_comms(catom_array);
		MPI::COMM_WORLD.Barrier();
	#endif

	// Set atom variables for simulation
	#ifdef MPICF
	// check for staged replicated data generation
	if(vmpi::replicated_data_staged==true && vmpi::mpi_mode==1){

		// copy ppn to constant temporary
		const int ppn=vmpi::ppn;

		for(int n=0;n<ppn;n++){
			bool i_am_it=false;
			if(vmpi::my_rank%ppn==n) i_am_it=true;

			// Only generate system if I am it
			if(i_am_it){

				zlog << zTs() << "Copying system data to optimised data structures on rank " << vmpi::my_rank << "..." << std::endl;
				//std::cerr << zTs() << "Copying system data to optimised data structures on rank " << vmpi::my_rank << "..." << std::endl;
				cs::set_atom_vars(catom_array,cneighbourlist);
				zlog << zTs() << "Copying on rank " << vmpi::my_rank << " completed." << std::endl;
				//std::cerr << zTs() << "Copying on rank " << vmpi::my_rank << " completed." << std::endl;
			}

			// Wait for process who is it
			MPI::COMM_WORLD.Barrier();

		} // end of loop over processes
	}
	else{
	#endif

	// Print informative message
	std::cout << "Copying system data to optimised data structures." << std::endl;
	zlog << zTs() << "Copying system data to optimised data structures." << std::endl;

	cs::set_atom_vars(catom_array,cneighbourlist);

	#ifdef MPICF
	} // stop if for staged generation here
	#endif

	// Set grain and cell variables for simulation
	grains::set_properties();
	cells::initialise();
	if(sim::hamiltonian_simulation_flags[4]==1) demag::init();

   // Determine number of local atoms
   #ifdef MPICF
      int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
   #else
      int num_local_atoms = atoms::num_atoms;
   #endif
   //----------------------------------------
   // Initialise local temperature data
   //----------------------------------------
   ltmp::initialise(cs::system_dimensions[0],
                  cs::system_dimensions[1],
                  cs::system_dimensions[2],
                  atoms::x_coord_array,
                  atoms::y_coord_array,
                  atoms::z_coord_array,
                  atoms::type_array,
                  num_local_atoms,
                  sim::Teq,
                  sim::pump_power,
                  sim::pump_time,
                  sim::TTG,
                  sim::TTCe,
                  sim::TTCl,
                  mp::dt_SI,
					   sim::Tmin,
					   sim::Tmax);

	//std::cout << num_atoms << std::endl;
	#ifdef MPICF
		//std::cout << "Outputting coordinate data" << std::endl;
		//vmpi::crystal_xyz(catom_array);
	int my_num_atoms=vmpi::num_core_atoms+vmpi::num_bdry_atoms;
	int total_num_atoms=0;
	MPI::COMM_WORLD.Reduce(&my_num_atoms,&total_num_atoms, 1,MPI_INT, MPI_SUM, 0 );
	std::cout << "Total number of atoms (all CPUs): " << total_num_atoms << std::endl;
   zlog << zTs() << "Total number of atoms (all CPUs): " << total_num_atoms << std::endl;
	#else
	std::cout << "Number of atoms generated: " << atoms::num_atoms << std::endl;
   zlog << zTs() << "Number of atoms generated: " << atoms::num_atoms << std::endl;

	#endif

	return EXIT_SUCCESS;
}

}
