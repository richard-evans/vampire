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
#include "material.hpp"
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
	double system_dimensions[3]={77.0,77.0,77.0};	// Size of system (A)
	double unit_cell_size[3]={3.54,3.54,3.54};		// Unit Cell Size (A) [Will eventually be local to unit cells]
	bool pbc[3]={false,false,false};						// Periodic boundary conditions
	unsigned int total_num_unit_cells[3]={0,0,0};	// Unit cells for entire system (x,y,z)
	unsigned int local_num_unit_cells[3]={0,0,0};	// Unit cells on local processor (x,y,z)
	std::string crystal_structure="sc";
	
	// System Parameters
	int particle_creation_parity=0; // Offset of particle centre (odd/even)
	double particle_scale=50.0;     // Diameter of particles/grains (A)
	double particle_spacing=10.0;   // Spacing Between particles (A)


	// Other directives and flags
	bool single_spin=false;
	int system_creation_flags[10]={0,0,0,0,0,0,0,0,0,0};
	std::string unit_cell_file="";
	
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
	demag::init();
	
	//std::cout << num_atoms << std::endl;
	#ifdef MPICF
		//std::cout << "Outputting coordinate data" << std::endl;
		//vmpi::crystal_xyz(catom_array);
	int my_num_atoms=vmpi::num_core_atoms+vmpi::num_bdry_atoms;
	int total_num_atoms=0;
	MPI::COMM_WORLD.Reduce(&my_num_atoms,&total_num_atoms, 1,MPI_INT, MPI_SUM, 0 );
	if(vmpi::my_rank==0){
	  std::cout << "Total number of atoms (all CPUs): " << total_num_atoms << std::endl;
	}

	#endif

	return EXIT_SUCCESS;
}

}
