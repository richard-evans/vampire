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
	std::vector<std::vector<int> > cneighbourlist; 

	// Set up Parallel Decomposition if required 
	#ifdef MPICF
		if(vmpi::mpi_mode==0) vmpi::geometric_decomposition(vmpi::num_processors,cs::system_dimensions);
	#endif

	//      Initialise variables for system creation	
	if(cs::system_creation_flags[0]==1){
		// read_coord_file();
	}
	// Create block of crystal of desired size
	cs::create_crystal_structure(catom_array);
	
	// Cut system to the correct type, species etc
	cs::create_system_type(catom_array);
	
	// Copy atoms for interprocessor communications
	#ifdef MPICF
	if(vmpi::mpi_mode==0){
		MPI::COMM_WORLD.Barrier();
		vmpi::copy_halo_atoms(catom_array);
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
		MPI::COMM_WORLD.Barrier();
		vmpi::identify_boundary_atoms(catom_array,cneighbourlist);
		MPI::COMM_WORLD.Barrier();
		vmpi::init_mpi_comms(catom_array);
		MPI::COMM_WORLD.Barrier();
	#endif

	//      Set atom variables for simulation
	cs::set_atom_vars(catom_array,cneighbourlist);

	//      Set grain and cell variables for simulation
	grains::set_properties();
	cells::initialise();
	demag::init();
	
	//      Generate system files for storage
	num_atoms=catom_array.size();
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

	#else
		if(1==0){
		std::ofstream xyz_file;
		xyz_file.open ("crystal.xyz");
		xyz_file << num_atoms+80 << std::endl;
		xyz_file << "" << std::endl;
	  	
	  	for(int atom=0; atom<num_atoms; atom++){
	  		xyz_file << material_parameters::material[catom_array[atom].material].element << "\t" << 
	  					catom_array[atom].x << "\t" << 
	  					catom_array[atom].y << "\t" << 
	  					catom_array[atom].z << "\t" << std::endl;
	  	}
	
		// Output axes
		for (int i=0;i<100;i+=5){
			xyz_file << "O\t" << float(i) << "\t" << 0.0 << "\t" << 0.0 << std::endl;
			xyz_file << "O\t" << 0.0 << "\t" << float(i) << "\t" << 0.0 << std::endl;
	
			xyz_file << "O\t" << cs::system_dimensions[0] << "\t" << cs::system_dimensions[1]-float(i) << "\t" << 0.0 << std::endl;
			xyz_file << "O\t" << cs::system_dimensions[0]-float(i) << "\t" << 	cs::system_dimensions[1] << "\t" << 0.0 << std::endl;
		}
		
		xyz_file.close();
		}
		#endif

	return EXIT_SUCCESS;
}

}
