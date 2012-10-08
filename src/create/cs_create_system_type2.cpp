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
#include "vio.hpp"
#include "vmath.hpp"

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
	int intermixing(std::vector<cs::catom_t> &);
	void dilute(std::vector<cs::catom_t> &);
	void geometry(std::vector<cs::catom_t> &);

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
	switch(cs::system_creation_flags[2]){
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
			err::vexit();
			break;
			
		default:{
			std::cerr << "Unknown system type requested, exiting" << std::endl;
			err::vexit();
			}
		}

		// call geometry function
		geometry(catom_array);
		
		// call intermixing function - must be before alloy function
		intermixing(catom_array);

		// call alloy function
		alloy(catom_array);

		// call dilution function
		dilute(catom_array);
		
		// Delete unneeded atoms
		clear_atoms(catom_array);
		
		// Check for zero atoms generated
		if(catom_array.size()==0){
			std::cerr << "Error, no atoms generated for requested system shape - increase system dimensions or reduce particle size!" << std::endl;
			err::vexit();
		}

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
	// find centre unit cell -- unsafe for large unit cells
	//particle_origin[0] = double(vmath::iround(cs::system_dimensions[0]/(2.0*unit_cell.dimensions[0])))*unit_cell.dimensions[0];
	//particle_origin[1] = double(vmath::iround(cs::system_dimensions[1]/(2.0*unit_cell.dimensions[1])))*unit_cell.dimensions[1];
	//particle_origin[2] = double(vmath::iround(cs::system_dimensions[2]/(2.0*unit_cell.dimensions[2])))*unit_cell.dimensions[2];

	particle_origin[0] = cs::system_dimensions[0]*0.5; 
	particle_origin[1] = cs::system_dimensions[1]*0.5;
	particle_origin[2] = cs::system_dimensions[2]*0.5;

	// check for move in particle origin and that unit cell size < 0.5 system size
	if(cs::particle_creation_parity==1 && 
		(2.0*unit_cell.dimensions[0]<cs::system_dimensions[0]) && 
		(2.0*unit_cell.dimensions[1]<cs::system_dimensions[1]) &&
		(2.0*unit_cell.dimensions[2]<cs::system_dimensions[2])){
		particle_origin[0]+=unit_cell.dimensions[0]*0.5;
		particle_origin[1]+=unit_cell.dimensions[1]*0.5;
		particle_origin[2]+=unit_cell.dimensions[2]*0.5;
	}
	
	// Use particle type flags to determine which particle shape to cut
	switch(cs::system_creation_flags[1]){
		case 0: // Bulk
			bulk(catom_array);
			break;
		case 1: // Cube
			cube(particle_origin,catom_array,0);
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
		case 6: // Teardrop
			tear_drop(particle_origin,catom_array,0);
			break;
		default:
			std::cout << "Unknown particle type requested for single particle system" << std::endl;
			err::vexit();
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
	const double repeat_size = cs::particle_scale+cs::particle_spacing;
	int num_x_particle = vmath::iround(cs::system_dimensions[0]/repeat_size);
	int num_y_particle = vmath::iround(cs::system_dimensions[1]/repeat_size);
	

	// Loop to generate cubic lattice points
	int particle_number=0;
	
	for (int x_particle=0;x_particle < num_x_particle;x_particle++){
		for (int y_particle=0;y_particle < num_y_particle;y_particle++){

			double particle_origin[3];
			// find centre unit cell
			//particle_origin[0] = double(iround(cs::system_dimensions[0]/(2.0*cs::unit_cell_size[0])))*cs::unit_cell_size[0];
			//particle_origin[1] = double(iround(cs::system_dimensions[1]/(2.0*cs::unit_cell_size[1])))*cs::unit_cell_size[1];
			//particle_origin[2] = double(iround(cs::system_dimensions[2]/(2.0*cs::unit_cell_size[2])))*cs::unit_cell_size[2];
			// Determine particle origin
			particle_origin[0] = double(x_particle)*repeat_size + repeat_size;
			particle_origin[1] = double(y_particle)*repeat_size + repeat_size;
			particle_origin[2] = double(vmath::iround(cs::system_dimensions[2]/(2.0*cs::unit_cell_size[2])))*cs::unit_cell_size[2];

			if(cs::particle_creation_parity==1){
				particle_origin[0]+=unit_cell.dimensions[0]*0.5;
				particle_origin[1]+=unit_cell.dimensions[1]*0.5;
				particle_origin[2]+=unit_cell.dimensions[2]*0.5;
			}
			// Check to see if a complete particle fits within the system bounds
			if((particle_origin[0]<(cs::system_dimensions[0]-cs::particle_scale)) &&
				(particle_origin[1]<(cs::system_dimensions[1]-cs::particle_scale))){

				// Use particle type flags to determine which particle shape to cut
				switch(cs::system_creation_flags[1]){
					case 0: // Bulk
						bulk(catom_array);
						break;
					case 1: // Cube
						cube(particle_origin,catom_array,particle_number);
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
					case 6: // Teardrop
						tear_drop(particle_origin,catom_array,particle_number);
						break;
					default:
						std::cout << "Unknown particle type requested for single particle system" << std::endl;
						err::vexit();
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
	
	num_x_particle = 2+iround(material_parameters::system_dimensions[0]/material_parameters::particle_scale);
	num_y_particle = 2+iround(material_parameters::system_dimensions[1]/material_parameters::particle_scale);
	
	//----------------------------------------------------------
	// Set integer equivalents of system/particle dimensions
	//----------------------------------------------------------
	
	//calculate int particle radius
	int_particle_scale = iround(material_parameters::particle_scale/(material_parameters::lattice_constant[1]));
	//int_particle_spacing = iround(material_parameters::particle_spacing/(2.0*material_parameters::lattice_constant[1]));
	
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
				particle_coords[0] = iround((1.0+particle_parity)*delta_particle_x_parity + delta_particle_x*x_particle);
				particle_coords[1] = iround((1.0+particle_parity)*delta_particle_y_parity + delta_particle_y*y_particle);
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
					err::vexit();
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

	zlog<< zTs() << "Determining alloy concentrations" << std::endl; 

	// loop over all atoms
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		// if atom material is alloy master then reassign according to % chance
		int local_material=catom_array[atom].material;
		if(mp::material[local_material].alloy_master==true){
		  // now check for unordered alloy
			if(mp::material[local_material].alloy_class==-1){
				//loop over all potential alloy materials
				for(int mat=0;mat<mp::num_materials;mat++){
					double probability = mp::material[local_material].alloy[mat];
					if(mtrandom::grnd() < probability){
						catom_array[atom].material=mat;
					}
				}
			}
			// if not ordered, then assume ordered
			else{
				// loop over all alloy materials
				for(int mat=0;mat<mp::num_materials;mat++){
					// get class of alloy material
					int alloy_class = mp::material[mat].alloy_class;
					// check for matching class and uc
					// ----- DOES NOT WORK for > 1 alloy!! 
					// need to check for correct alloy master material ------------
					if(catom_array[atom].uc_category==alloy_class){
						// set material
						catom_array[atom].material=mat;
					}
				}
			}
		}
	}

	// Determine number of atoms of each class and output to log
	std::vector<unsigned int> MaterialNumbers(mp::num_materials,0);
	for(unsigned int atom=0;atom<catom_array.size();atom++) MaterialNumbers.at(catom_array[atom].material)++;
	for(int mat=0;mat<mp::num_materials;mat++) zlog << zTs() << "Material " << mat << " " << mp::material[mat].name << " makes up " << double(MaterialNumbers[mat])*100.0/double(catom_array.size()) << "% of all atoms." << std::endl;

	return EXIT_SUCCESS;	
}

int intermixing(std::vector<cs::catom_t> & catom_array){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::intermixing has been called" << std::endl;}      

	// loop over all atoms
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		// get current material
		int current_material=catom_array[atom].material;
		int final_material=current_material;

		//loop over all potential intermixing materials
		for(int mat=0;mat<mp::num_materials;mat++){
			if(mp::material[current_material].intermixing[mat]>0.0){
				// find which region atom is in and test for probability of different material
				double z=catom_array[atom].z;
				double min = mp::material[mat].min*cs::system_dimensions[2];
				double max = mp::material[mat].max*cs::system_dimensions[2];
				double mean = (min+max)/2.0;
				if(z<=min){
					double probability=0.5+0.5*tanh((z-min)/(mp::material[current_material].intermixing[mat]*cs::system_dimensions[2]));
					if(mtrandom::grnd() < probability) final_material=mat;
				}
				else if(z>min && z<=mean){
					double probability=0.5+0.5*tanh((z-min)/(mp::material[current_material].intermixing[mat]*cs::system_dimensions[2]));
					if(mtrandom::grnd() < probability) final_material=mat;
				}
				else if(z>mean && z<=max){
					double probability=0.5-0.5*tanh((z-max)/(mp::material[current_material].intermixing[mat]*cs::system_dimensions[2]));
					if(mtrandom::grnd() < probability) final_material=mat;
				}
				else if(z>max){
					double probability=0.5-0.5*tanh((z-max)/(mp::material[current_material].intermixing[mat]*cs::system_dimensions[2]));
					//std::cout << current_material << "\t" << mat << "\t" << atom << "\t" << z << "\t" << max << "\t" << probability << std::endl;
					if(mtrandom::grnd() < probability) final_material=mat;
				}
			}
		}
		
		// set final material
		catom_array[atom].material=final_material;
	}
	
	return EXIT_SUCCESS;    
}

void dilute (std::vector<cs::catom_t> & catom_array){
    // check calling of routine if error checking is activated
    if(err::check==true){std::cout << "cs::dilute has been called" << std::endl;}
    
    // loop over all atoms
    for(unsigned int atom=0;atom<catom_array.size();atom++){
      // if atom material is alloy master
      int local_material=catom_array[atom].material;
      double probability = mp::material[local_material].density;
      if(mtrandom::grnd() > probability) catom_array[atom].include=false;
    }

    return;
  }
  
void geometry (std::vector<cs::catom_t> & catom_array){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::geometry has been called" << std::endl;}
    
	// Check for any geometry
	bool cut=false;
	
	for(int mat=0; mat<mp::num_materials; mat++){
		if(mp::material[mat].geometry>0) cut=true;
	}
	
	// Return from function if no geometry is defined. 
	if(cut==false) return;
	
	// Otherwise proceed
	zlog << zTs() << "Cutting materials within defined geometry." << std::endl; 
	
	// Check for force material type by geometry 
	if(cs::SelectMaterialByGeometry==false){
		
		// loop over all atoms
		for(unsigned int atom=0;atom<catom_array.size();atom++){
			
			// check for geometry information
			const int geo=mp::material[catom_array[atom].material].geometry;

			// if exists, then remove atoms outside polygon
			if(geo!=0){
				double x = catom_array[atom].x;
				double y = catom_array[atom].y;
				std::vector<double> px(geo);
				std::vector<double> py(geo);
				// Initialise polygon points
				for(int p=0;p<geo;p++){
					px[p]=mp::material[catom_array[atom].material].geometry_coords[p][0]*cs::system_dimensions[0];
					py[p]=mp::material[catom_array[atom].material].geometry_coords[p][1]*cs::system_dimensions[1];
				}
				// check if point is outside of polygon, if so delete it 
				if(vmath::point_in_polygon2(x,y,px,py,geo)==false){
					catom_array[atom].include=false;
				}
			}
			
		}

	}
	else{
	
		// Re-identify all atoms as material 0 
		for(unsigned int atom=0;atom<catom_array.size();atom++){
			catom_array[atom].material=0;
		}
		
		// loop over all materials and include accoring to geometry

		// determine z-bounds for materials
		std::vector<double> mat_min(mp::num_materials);
		std::vector<double> mat_max(mp::num_materials);

		// initialise array to hold z-bound information in real space
		for(int mat=0;mat<mp::num_materials;mat++){
			mat_min[mat]=mp::material[mat].min*cs::system_dimensions[2];
			mat_max[mat]=mp::material[mat].max*cs::system_dimensions[2];
			// alloys generally are not defined by height, and so have max = 0.0
			if(mat_max[mat]<0.0000001) mat_max[mat]=-0.1;
		}
		
		// Loop over all materials	
		for(int mat=0;mat<mp::num_materials;mat++){

			// check for geometry information
			const int geo=mp::material[mat].geometry;

			// if geometry information exists, then include atoms inside polygon and within material height
			if(geo!=0){

				// create array to store geoemtric points
				std::vector<double> px(geo);
				std::vector<double> py(geo);
				// Initialise polygon points
				for(int p=0;p<geo;p++){
					px[p]=mp::material[mat].geometry_coords[p][0]*cs::system_dimensions[0];
					py[p]=mp::material[mat].geometry_coords[p][1]*cs::system_dimensions[1];
				}
			
				for(unsigned int atom=0;atom<catom_array.size();atom++){
					double x = catom_array[atom].x;
					double y = catom_array[atom].y;
					double z = catom_array[atom].z;

					const double cz=catom_array[atom].z;
					if((cz>=mat_min[mat]) && (cz<mat_max[mat]) && (vmath::point_in_polygon2(x,y,px,py,geo)==true)){
						catom_array[atom].material=mat;
						catom_array[atom].include=true;
					}
				}
			}
		}
	}
	
	return;
}
  
} // end of namespace
