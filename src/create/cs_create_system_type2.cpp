//======================================================================
//                         create_system_type
//   Subroutine to set system size and create desired crystal structure
//
//======================================================================

#include "errors.hpp"
#include "create.hpp"
#include "grains.hpp"
#include "material.hpp"
#include "random.hpp"
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <list>

namespace cs{
	
	
	
	//======================================================================
	//                         create_system_type
	//   Subroutine to set system size and create desired crystal structure
	//
	//======================================================================	
int create_system_type(std::vector<cs::catom_t> & catom_array){

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::create_system_type has been called" << std::endl;}
	
	//----------------------------------------
	// function prototypes
	//----------------------------------------
	
	int particle(std::vector<cs::catom_t> &);
	int particle_array(std::vector<cs::catom_t> &);
	
	int alloy(std::vector<cs::catom_t> &);
	//int particle_array(int,int**,int*);
	//int hex_particle_array(int,int**,int*);
	//int voronoi_film(int**,double*,int*);
	//int pop_template_2D(int**,int,int**,int*,int*);
	//int cs_calc_grain_vol(int,int*);
	//int multilayer(int,int**,int*,int*);
	//int core_shell(int,int**,int*,int*);
	//int mpi_create_system_type(int,int**,int*,int*);
	//----------------------------------------
	// Local variables
	//----------------------------------------
	
	//int* particle_include_array;
	//int** template_array_2D;
	//int atom;
	//int new_num_atoms;


	
	//----------------------------------------------------------------------------------
	// Choose which system type to create
	//----------------------------------------------------------------------------------
	switch(material_parameters::system_creation_flags[2]){
		case 0: // Isolated particle
			particle(catom_array);
			break;
		
		case 1: // Cubic Particle Array
			particle_array(catom_array);
			break;
		
		case 2: // Hexagonal Particle Array
			//hex_particle_array(catom_array);
			break;
			
		case 3: // Voronoi Granular Film
			voronoi_film(catom_array);
			break;
			
		case 4: // Grain Growth Method
			//grain_growth(cs_num_atoms,cs_coord_array,particle_include_array,cs_atom_type_array);
			std::cerr << "Grain growth not yet implemented, exiting" << std::endl;
			exit(0);
			break;
			
		default:{
			std::cerr << "Unknown system type requested, exiting" << std::endl;
			exit (1);
			}
		}

		alloy(catom_array);



		// Delete unneeded atoms
		clear_atoms(catom_array);
		
	//----------------------------------------------------------------------------------
	// Set post construction atom identifiers for multispecies system
	//----------------------------------------------------------------------------------
	//if(material_parameters::system_creation_flags[4]==1){
	//	multilayer(cs_num_atoms,cs_coord_array,cs_atom_type_array,particle_include_array);
	//}
	//else if(material_parameters::system_creation_flags[4]==2){
	//	intermixing(cs_num_atoms,cs_coord_array,cs_atom_type_array,particle_include_array);
	//}

	//-------------------------------------------------------------
	//      Clean Up Parallel Decomposition if required 
	//-------------------------------------------------------------
	//#ifdef MPICF
	//	if(mpi_generic::mpi_mode==2){
	//		mpi_create_system_type(cs_num_atoms,cs_coord_array,cs_atom_type_array,particle_include_array);
	//	}
	//#endif

	//----------------------------------------------------------------------------------
	// Calculate new number of atoms and re-sort cs_* arrays according to included list
	//----------------------------------------------------------------------------------
	/*new_num_atoms = 0;
	std::vector <int> tmp_grain_no_array(cs_num_atoms);
	
	for(atom=0;atom<cs_num_atoms;atom++){
		if(particle_include_array[atom]!=-1){
			tmp_grain_no_array[new_num_atoms]=particle_include_array[atom];
			cs_coord_array[new_num_atoms][0] = cs_coord_array[atom][0];
			cs_coord_array[new_num_atoms][1] = cs_coord_array[atom][1];			
			cs_coord_array[new_num_atoms][2] = cs_coord_array[atom][2];
			cs_atom_type_array[new_num_atoms] = cs_atom_type_array[atom];
			#ifdef MPICF
				using namespace mpi_create_variables;
				mpi_atom_comm_class_array[new_num_atoms]=mpi_atom_comm_class_array[atom];
				mpi_atom_location_array[new_num_atoms]=mpi_atom_location_array[atom];
				mpi_atom_global_coord_array[3*new_num_atoms+0]=mpi_atom_global_coord_array[3*atom+0];
				mpi_atom_global_coord_array[3*new_num_atoms+1]=mpi_atom_global_coord_array[3*atom+1];
				mpi_atom_global_coord_array[3*new_num_atoms+2]=mpi_atom_global_coord_array[3*atom+2];
			#endif
			new_num_atoms++;
		}	
	}*/
/*
	// Save grain numbers
	//std::cout << "Saving grain numbers for " << grains::num_grains << std::endl;
	atoms::grain_array.resize(new_num_atoms,0);
	for(atom=0;atom<new_num_atoms;atom++){
		int grain = tmp_grain_no_array[atom];
		if((grain >= 0) && (grain < (grains::num_grains + mp::num_materials))){
			atoms::grain_array[atom]=grain;
		}
		else {
			std::cerr << "Error - Unknown grain number " << grain << " for atom " << atom << std::endl;
			exit(1);
		}
	}*/
	//---------------------
	// Reset cs_num_atoms
	//---------------------
	//cs_num_atoms = new_num_atoms;
	
	//std::cout << "num atoms\t" << cs_num_atoms << std::endl; 

	return 0;
}

int particle(std::vector<cs::catom_t> & catom_array){
	//====================================================================================
	//
	//										particle
	//
	//					Subroutine to cut single particle from lattice
	//
	//							Version 1.0 R Evans 22/09/2008
	//
	//====================================================================================

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::particle has been called" << std::endl;}	

	//---------------------------------------------------
	// Set particle origin to atom at centre of lattice
	//---------------------------------------------------

	double particle_origin[3];
	// find centre unit cell
	particle_origin[0] = double(round(mp::system_dimensions[0]/(2.0*mp::lattice_constant[0])))*mp::lattice_constant[0];
	particle_origin[1] = double(round(mp::system_dimensions[1]/(2.0*mp::lattice_constant[1])))*mp::lattice_constant[1];
	particle_origin[2] = double(round(mp::system_dimensions[2]/(2.0*mp::lattice_constant[2])))*mp::lattice_constant[2];

	if(mp::particle_creation_parity==1){
		particle_origin[0]+=mp::lattice_constant[0]*0.5;
		particle_origin[1]+=mp::lattice_constant[1]*0.5;
		particle_origin[2]+=mp::lattice_constant[2]*0.5;
	}
	
	// Use particle type flags to determine which particle shape to cut
	switch(material_parameters::system_creation_flags[1]){
		case 0: // Bulk
			bulk(catom_array,0);
			break;
		case 2: // Cylinder
			cylinder(particle_origin,catom_array,0);
			break;
		case 4: // Sphere
			sphere(particle_origin,catom_array,0);
			break;
		case 5: // Truncated Octahedron
			truncated_octahedron(particle_origin,catom_array,0);
			break;
		default:
			std::cout << "Unknown particle type requested for single particle system" << std::endl;
			exit(1);
	}

	return EXIT_SUCCESS;	
}

int particle_array(std::vector<cs::catom_t> & catom_array){
	//====================================================================================
	//
	//									particle_array
	//
	//					Subroutine to cut many particles from lattice
	//					in a cubic configuration
	//
	//							Version 1.0 R Evans 23/09/2008
	//
	//					Note: particle counter doesn't do anything (yet)
	//
	//====================================================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::particle_array has been called" << std::endl;}	

	// Set number of particles in x and y directions
	const double repeat_size = mp::particle_scale+mp::particle_spacing;
	int num_x_particle = round(mp::system_dimensions[0]/repeat_size);
	int num_y_particle = round(mp::system_dimensions[1]/repeat_size);
	

	// Loop to generate cubic lattice points
	int particle_number=0;
	
	for (int x_particle=0;x_particle < num_x_particle;x_particle++){
		for (int y_particle=0;y_particle < num_y_particle;y_particle++){

			double particle_origin[3];
			// find centre unit cell
			//particle_origin[0] = double(round(mp::system_dimensions[0]/(2.0*mp::lattice_constant[0])))*mp::lattice_constant[0];
			//particle_origin[1] = double(round(mp::system_dimensions[1]/(2.0*mp::lattice_constant[1])))*mp::lattice_constant[1];
			//particle_origin[2] = double(round(mp::system_dimensions[2]/(2.0*mp::lattice_constant[2])))*mp::lattice_constant[2];
			// Determine particle origin
			particle_origin[0] = double(x_particle)*repeat_size + repeat_size;
			particle_origin[1] = double(y_particle)*repeat_size + repeat_size;
			particle_origin[2] = double(round(mp::system_dimensions[2]/(2.0*mp::lattice_constant[2])))*mp::lattice_constant[2];

			if(mp::particle_creation_parity==1){
				particle_origin[0]+=mp::lattice_constant[0]*0.5;
				particle_origin[1]+=mp::lattice_constant[1]*0.5;
				particle_origin[2]+=mp::lattice_constant[2]*0.5;
			}
			// Check to see if a complete particle fits within the system bounds
			if((particle_origin[0]<(mp::system_dimensions[0]-mp::particle_scale)) &&
				(particle_origin[1]<(mp::system_dimensions[1]-mp::particle_scale))){

				// Use particle type flags to determine which particle shape to cut
				switch(material_parameters::system_creation_flags[1]){
					case 0: // Bulk
						bulk(catom_array,particle_number);
						break;
					case 2: // Cylinder
						cylinder(particle_origin,catom_array,particle_number);
						break;
					case 4: // Sphere
						sphere(particle_origin,catom_array,particle_number);
						break;
					case 5: // Truncated Octahedron
						truncated_octahedron(particle_origin,catom_array,particle_number);
						break;
					default:
						std::cout << "Unknown particle type requested for single particle system" << std::endl;
						exit(1);
				}
				// Increment Particle Number Counter
				particle_number++;
			}
		}
	}
	grains::num_grains = particle_number;
	
	// Clear unneeded atoms
	//clear_atoms(catom_array);
	// Re-order atoms by particle number
	sort_atoms_by_grain(catom_array);
	
	return EXIT_SUCCESS;	
}
/*
int hex_particle_array(int cs_num_atoms,int** cs_coord_array,int* particle_include_array){
	//====================================================================================
	//
	//									particle_array
	//
	//					Subroutine to cut many particles from lattice
	//					in a cubic configuration
	//
	//							Version 1.0 R Evans 28/10/2008
	//
	//					Note: particle counter doesn't do anything (yet)
	//
	//====================================================================================

	//---------------------------------------------------
	// Local variables
	//---------------------------------------------------
	
	int particle_origin[3];
	int atom;
	int particle_number, x_particle, y_particle, num_x_particle, num_y_particle;
	int particle_parity;
	//int distance_squared, min_distance_squared, closest_atom;
	
	int particle_coords[3], int_particle_scale,int_particle_spacing;
	int int_system_dimensions[3];
	
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "particle_array has been called" << std::endl;}	

	//----------------------------------------------------------
	// set initial particle number
	//----------------------------------------------------------
	
	particle_number = 0;
	
	//----------------------------------------------------------
	// Set number of particles in x and y directions
	//----------------------------------------------------------
	
	num_x_particle = 2+round(material_parameters::system_dimensions[0]/material_parameters::particle_scale);
	num_y_particle = 2+round(material_parameters::system_dimensions[1]/material_parameters::particle_scale);
	
	//----------------------------------------------------------
	// Set integer equivalents of system/particle dimensions
	//----------------------------------------------------------
	
	//calculate int particle radius
	int_particle_scale = round(material_parameters::particle_scale/(material_parameters::lattice_constant[1]));
	//int_particle_spacing = round(material_parameters::particle_spacing/(2.0*material_parameters::lattice_constant[1]));
	
	//double delta_particle_x_const = material_parameters::particle_scale/(material_parameters::lattice_constant[1]);
	//double delta_particle_y_const = material_parameters::particle_scale/(material_parameters::lattice_constant[2])*sqrt(3.0);
	double delta_particle_x = (material_parameters::particle_scale + material_parameters::particle_spacing)/(material_parameters::lattice_constant[1]);
	double delta_particle_y = ((material_parameters::particle_scale + material_parameters::particle_spacing)/(material_parameters::lattice_constant[2]))*sqrt(3.0);
	double delta_particle_x_parity = delta_particle_x*0.5;
	double delta_particle_y_parity = delta_particle_y*0.5;
	
	int_system_dimensions[0] = int(material_parameters::system_dimensions[0]/material_parameters::lattice_constant[0]);
	int_system_dimensions[1] = int(material_parameters::system_dimensions[1]/material_parameters::lattice_constant[1]);
	int_system_dimensions[2] = int(material_parameters::system_dimensions[2]/material_parameters::lattice_constant[2]);
	
	//---------------------------------------------------------
	// Loop to generate cubic hexagonal lattice points
	//---------------------------------------------------------

	for (x_particle=0;x_particle < num_x_particle;x_particle++){
		for (y_particle=0;y_particle < num_y_particle;y_particle++){
			for (particle_parity=0;particle_parity<2;particle_parity++){

				//---------------------------------------------------
				// Determine particle coordinates
				//---------------------------------------------------
				
				//particle_coords[0] = 2*x_particle*(int_particle_scale + int_particle_spacing) + 2*int_particle_scale + (particle_parity*(int_particle_scale + int_particle_spacing));
				//particle_coords[1] = 2*y_particle*(int_particle_scale + int_particle_spacing) + 2*int_particle_scale + (particle_parity*(int_particle_scale + int_particle_spacing));
				//particle_coords[2] = int_system_dimensions[2]/2;
				particle_coords[0] = round((1.0+particle_parity)*delta_particle_x_parity + delta_particle_x*x_particle);
				particle_coords[1] = round((1.0+particle_parity)*delta_particle_y_parity + delta_particle_y*y_particle);
				particle_coords[2] = int_system_dimensions[2]/2;

		//---------------------------------------------------
		// Set particle origin to particle coords
		//---------------------------------------------------
	
		particle_origin[0] = 2*particle_coords[0];
		particle_origin[1] = 6*particle_coords[1];
		particle_origin[2] = 2*particle_coords[2];
		
		if(mp::particle_creation_parity==1){
			std::string cs = "sc";
			if(mp::crystal_structure==cs){
				particle_origin[0]+=1;
				particle_origin[1]+=3;
				particle_origin[2]+=1;
			}
			cs = "fcc";
			if(mp::crystal_structure==cs){
				particle_origin[0]+=1;
				particle_origin[1]+=3;
				particle_origin[2]+=1;
			}
		}
		//-------------------------------------------------------------------
		// Check to see if a complete particle fits within the system bounds
		//-------------------------------------------------------------------
		
		//std::cout << particle_origin[0] << "\t" << 2*(int_system_dimensions[0]-int_particle_scale) << "\t"
		//	 << particle_origin[1] << "\t" << 6*(int_system_dimensions[1]-int_particle_scale) << std::endl;
		
		
		if((particle_origin[0]<2*(int_system_dimensions[0]-int_particle_scale)) &&
		   (particle_origin[1]<6*(int_system_dimensions[1]-int_particle_scale))){
			//------------------------------------------------------------------
			// Use particle type flags to determine which particle shape to cut
			//------------------------------------------------------------------
			switch(material_parameters::system_creation_flags[1]){
				case 0: // Bulk
					for(atom=0;atom<cs_num_atoms;atom++) particle_include_array[atom]=1;
					std::cout << "Warning - bulk particle type requested for particle array system" << std::endl;
					return 0;
					break;
				case 2: // Cylinder
					cs_cylinder(cs_num_atoms,cs_coord_array,particle_include_array,particle_origin);
					break;
				case 4: // Sphere
					cs_sphere(cs_num_atoms,cs_coord_array,particle_include_array,particle_origin);
					break;			
				case 5: // Truncated Octahedron
					cs_truncated_octahedron(cs_num_atoms,cs_coord_array,particle_include_array,particle_origin);
					break;
					
				default:
					std::cout << "Unknown particle type requested for single particle system" << std::endl;
					exit(1);
				}
				
			//------------------------------------------------------------------
			// Increment Particle Number Counter
			//------------------------------------------------------------------		
			particle_number++;
		   }
		   }
		}
	}
		
	return 0;	
}
*/
/*

int pop_template_2D(int** template_array_2D,int cs_num_atoms,int** cs_coord_array,int* particle_include_array,int* template_bounds){

	//====================================================================================
	//
	//											pop_template_2D
	//
	//					Subroutine to map 2D template onto particles include array
	//
	//									Version 0.2 R Evans 16/09/2008
	//
	//====================================================================================

  	for(int atom=0;atom<cs_num_atoms;atom++){
		#ifdef MPICF
		int cx=cs_coord_array[atom][0]+mpi_create_variables::int_mpi_offset[0];
		int cy=cs_coord_array[atom][1]+mpi_create_variables::int_mpi_offset[1];
		#else
		int cx=cs_coord_array[atom][0];
		int cy=cs_coord_array[atom][1];
		#endif

		// Check for coords in template bounds
		if((cx>0) && (cx<template_bounds[0]) && (cy>0) && (cy<template_bounds[1])){
			if(template_array_2D[cx][cy]!=-1){
				particle_include_array[atom]=template_array_2D[cx][cy];
			}
		}
	}

	return 0;

}
*/

int clear_atoms(std::vector<cs::catom_t> & catom_array){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::clear_atoms has been called" << std::endl;}	
	
	// Get original and new number of atoms
	const int num_atoms=catom_array.size();
	int num_included=0;
	for(int a=0;a<num_atoms;a++){
		if(catom_array[a].include==true){
			num_included++;
		}
	}
	
	// check for unneeded
	if(num_atoms!=num_included){
		std::vector<cs::catom_t> tmp_catom_array(num_atoms);
		tmp_catom_array=catom_array;
		catom_array.resize(num_included);
		int atom=0;
		for(int a=0;a<num_atoms;a++){
			if(catom_array[a].include==true){
				catom_array[atom]=tmp_catom_array[a];
				atom++;
			}
		}
		tmp_catom_array.resize(0);
	}
	
	return EXIT_SUCCESS;
}

// comparison function
bool compare(cs::catom_t first,cs::catom_t second){
	if(first.grain<second.grain) return true;
	else return false;
}
	
int sort_atoms_by_grain(std::vector<cs::catom_t> & catom_array){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::sort_atoms_by_grain has been called" << std::endl;}	
	
	// Get number of atoms
	const int num_atoms=catom_array.size();
	
	// Create list object
	std::list <cs::catom_t> catom_list(num_atoms);

	// copy data to list
	copy(catom_array.begin(), catom_array.end(), catom_list.begin());

	// sort date in list
	catom_list.sort(compare);
	
	// copy list to data
	copy(catom_list.begin(), catom_list.end(), catom_array.begin());
	
	return EXIT_SUCCESS;
}

int alloy(std::vector<cs::catom_t> & catom_array){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::alloy has been called" << std::endl;}	


	// loop over all atoms
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		// if atom material is alloy master
		int local_material=catom_array[atom].material;
		if(mp::material[local_material].alloy_master==true){
			//loop over all potential alloy materials
			for(int mat=0;mat<mp::num_materials;mat++){
				double probability = mp::material[local_material].alloy[mat];
				if(mtrandom::grnd() < probability){
					catom_array[atom].material=mat;
				}
			}
		}
	}
				
	return EXIT_SUCCESS;	
}





}
