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
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
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
#include <iostream>
#include <fstream>
#include "public.hpp"
#include "atoms.hpp"
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

int create(){
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(error_checking::error_check==true){std::cout << "cs::create has been called" << std::endl;}
	
	//=============================================================
	//      System creation variables
	//=============================================================
	
	int num_atoms=0;							///> num_atoms for generation routines
	//int supercell_dim[3];				///> integer dimesions of system supercells
	//double supercell_size[3];				///> real supercell size (Angstroms)
	
	// Atom creation array
	std::vector<cs::catom_t> catom_array; 
	std::vector<std::vector<int> > cneighbourlist; 

	//std::vector<std::vector<int> > vecvec;
	//vecvec.push_back(std::vector<int>());
	//vecvec[0].push_back(3);
	//int i = vecvec[0][0];
	//OK, vecvec.size() == 1; vecvec[0].size() == 1, vecvec[0][0] == 3

	
	//=============================================================
	//      Set up Parallel Decomposition if required 
	//=============================================================
	#ifdef MPICF
		if(mpi_generic::mpi_mode==2) vmpi::geometric_decomposition(vmpi::num_processors,mp::system_dimensions);
	#endif

	//=============================================================
	//      Initialise variables for system creation
	//=============================================================
	
	if(mp::system_creation_flags[0]==1){
		// read_coord_file();
	}
	else{
		//=============================================================
		//      Create block of crystal of desired size
		//=============================================================

		cs::create_crystal_structure(catom_array);
     
		//=============================================================
		//      Cut system to the correct type, species etc
		//=============================================================

		cs::create_system_type(catom_array);
	}
	//=============================================================
	//      Create Neighbour list for system
	//=============================================================
	// Copy atoms for interprocessor communications
	#ifdef MPICF
		vmpi::copy_halo_atoms(catom_array);
	#else
		//cs::copy_periodic_boundaries(catom_array);
	#endif
	
	cs::create_neighbourlist(catom_array,cneighbourlist);
	
	#ifdef MPICF
		vmpi::identify_boundary_atoms(catom_array,cneighbourlist);
		vmpi::init_mpi_comms(catom_array);
	#endif

 	//=============================================================
	//      Set atom variables for simulation
	//=============================================================
	
	cs::set_atom_vars(catom_array,cneighbourlist);

	//=============================================================
	//      Generate system files for storage
	//=============================================================
	num_atoms=catom_array.size();
	//std::cout << num_atoms << std::endl;
	#ifdef MPICF
		//std::cout << "Outputting coordinate data" << std::endl;
		vmpi::crystal_xyz(catom_array);
	#else
		if(1==1){
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
	
			xyz_file << "O\t" << material_parameters::system_dimensions[0] << "\t" << material_parameters::system_dimensions[1]-float(i) << "\t" << 0.0 << std::endl;
			xyz_file << "O\t" << material_parameters::system_dimensions[0]-float(i) << "\t" << 	material_parameters::system_dimensions[1] << "\t" << 0.0 << std::endl;
		}
		
		xyz_file.close();
		}
		#endif
		//exit(0);
		//exit(1);
	/*if(1==2){
	std::ofstream exc_file;
  	exc_file.open ("exchange.txt");
  	exc_file << cs_num_atoms << std::endl;
  	exc_file << "" << std::endl;
  	
  	for(atom=0; atom<cs_num_atoms; atom++){
  		exc_file << atom << "\t";
  		for(int nn=0;nn<material_parameters::hamiltonian_num_neighbours;nn++){ 
  					exc_file << cs_neighbourlist_array[atom][nn]<< "\t";
  		}
  		exc_file << std::endl;
  	}
  	exc_file.close();
	}*/

	//#endif
  	//=============================================================
	//      Deallocate allocated arrays close output files etc
	//=============================================================

  	//--------------------------
	// deallocate cs_coord_array
	//--------------------------
	/*
  	try{for(int i=0; i<init_cs_num_atoms ; i++) delete [] cs_coord_array[i];
    	delete [] cs_coord_array;
    	cs_coord_array=NULL;
    	}
  	catch(...){std::cout << "error deallocating cs_coord_array" << std::endl;exit(1);}
  	
	//-------------------------------
	// deallocate cs_atom_type_array
	//-------------------------------
	
  	try{delete [] cs_atom_type_array;
  		cs_atom_type_array=NULL;
  		}
  	catch(...){std::cout << "error deallocating cs_atom_type_array" << std::endl;exit(1);}

	//-----------------------------------
	// deallocate cs_neighbourlist_array
	//-----------------------------------
	
  	try{for(int i=0; i<cs_num_atoms ; i++) delete [] cs_neighbourlist_array[i];
    	delete [] cs_neighbourlist_array;
    	cs_neighbourlist_array=NULL;
    	}
  	catch(...){std::cout << "error deallocating cs_neighbourlist_array" << std::endl;exit(1);}	
  	
	
	if(final_num_atoms>0){
		if(mpi_generic::my_rank==0){
			std::cout << "System Generated Successfully" << std::endl;
		}
	}
	else{
		std::cerr << "No atoms generated, exiting!" << std::endl;
		exit(1);
	}*/
	return EXIT_SUCCESS;
}

}
