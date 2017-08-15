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
//======================================================================
//                         create_system_type
//   Subroutine to set system size and create desired crystal structure
//
//======================================================================

// C++ standard library headers
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <list>

// Vampire headers
#include "errors.hpp"
#include "create.hpp"
#include "grains.hpp"
#include "material.hpp"
#include "random.hpp"
#include "unitcell.hpp"
#include "vio.hpp"
#include "vmath.hpp"

// Internal create header
#include "internal.hpp"

namespace cs{

   //----------------------------------------
   // function prototypes
   //----------------------------------------
   int particle(std::vector<cs::catom_t> &);
   int particle_array(std::vector<cs::catom_t> &);

   int intermixing(std::vector<cs::catom_t> &);
   void dilute(std::vector<cs::catom_t> &);
   void geometry(std::vector<cs::catom_t> &);
   void fill(std::vector<cs::catom_t> &);
   void roughness(std::vector<cs::catom_t> &);
   void calculate_atomic_composition(std::vector<cs::catom_t> &);
   void centre_particle_on_atom(std::vector<double>& particle_origin, std::vector<cs::catom_t>& catom_array);

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
	// Local variables
	//----------------------------------------

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

      // Check for voronoi construction and apply before csg operations
      if(create::internal::generate_voronoi_substructure) create::internal::voronoi_substructure(catom_array);

		// call fill function to fill in void
		fill(catom_array);

		// call geometry function
		geometry(catom_array);

		// call intermixing function - must be before alloy function
		intermixing(catom_array);

		// call surface roughness function
		// Formally called here but now moved to crystal structure generation
		//roughness(catom_array);

		// call alloy function
		create::internal::alloy(catom_array);

		// call dilution function
		dilute(catom_array);

		// Delete unneeded atoms
		clear_atoms(catom_array);

		// Calculate final atomic composition
		calculate_atomic_composition(catom_array);

		// Check for zero atoms generated
		if(catom_array.size()==0){
			terminaltextcolor(RED);
			std::cerr << "Error, no atoms generated for requested system shape - increase system dimensions or reduce particle size!" << std::endl;
			terminaltextcolor(WHITE);
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

	std::vector<double> particle_origin(3,0.0);

	particle_origin[0] = cs::system_dimensions[0]*0.5;
	particle_origin[1] = cs::system_dimensions[1]*0.5;
	particle_origin[2] = cs::system_dimensions[2]*0.5;

   centre_particle_on_atom(particle_origin, catom_array);

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
			create::internal::bulk(catom_array);
			break;
		case 1: // Cube
			create::internal::cube(particle_origin,catom_array,0);
			break;
		case 2: // Cylinder
			create::internal::cylinder(particle_origin,catom_array,0);
			break;
      case 3: // Ellipsoid
         create::internal::ellipsoid(particle_origin,catom_array,0);
         break;
		case 4: // Sphere
			create::internal::sphere(particle_origin,catom_array,0);
			break;
		case 5: // Truncated Octahedron
			create::internal::truncated_octahedron(particle_origin,catom_array,0);
			break;
		case 6: // Teardrop
			create::internal::teardrop(particle_origin,catom_array,0);
			break;
      case 7: // Faceted particle
   		create::internal::faceted(particle_origin,catom_array,0);
   		break;
		case 8: // Cone
			create::internal::cone(particle_origin,catom_array,0);
			break;
      case 9: // Bubble
         create::internal::bubble(particle_origin,catom_array,0);
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
	//
	//====================================================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::particle_array has been called" << std::endl;}

	// Set number of particles in x and y directions
	const double repeat_size = cs::particle_scale+cs::particle_spacing;
	int num_x_particle = vmath::iceil(cs::system_dimensions[0]/repeat_size);
	int num_y_particle = vmath::iceil(cs::system_dimensions[1]/repeat_size);

	// Loop to generate cubic lattice points
	int particle_number=0;

	std::vector<double> particle_origin(3,0.0);

	for (int x_particle=0;x_particle < num_x_particle;x_particle++){
		for (int y_particle=0;y_particle < num_y_particle;y_particle++){

			// Determine particle origin
			particle_origin[0] = double(x_particle)*repeat_size + cs::particle_scale*0.5 + cs::particle_array_offset_x;
			particle_origin[1] = double(y_particle)*repeat_size + cs::particle_scale*0.5 + cs::particle_array_offset_y;
			particle_origin[2] = double(vmath::iround(cs::system_dimensions[2]/(2.0*unit_cell.dimensions[2])))*unit_cell.dimensions[2];

         centre_particle_on_atom(particle_origin, catom_array);

			if(cs::particle_creation_parity==1){
				particle_origin[0]+=unit_cell.dimensions[0]*0.5;
				particle_origin[1]+=unit_cell.dimensions[1]*0.5;
				particle_origin[2]+=unit_cell.dimensions[2]*0.5;
			}
			// Check to see if a complete particle fits within the system bounds
			if((particle_origin[0]<=(cs::system_dimensions[0]-cs::particle_scale*0.5)) &&
				(particle_origin[1]<=(cs::system_dimensions[1]-cs::particle_scale*0.5))){

				// Use particle type flags to determine which particle shape to cut
				switch(cs::system_creation_flags[1]){
					case 0: // Bulk
						create::internal::bulk(catom_array);
						break;
					case 1: // Cube
						create::internal::cube(particle_origin,catom_array,particle_number);
						break;
					case 2: // Cylinder
						create::internal::cylinder(particle_origin,catom_array,particle_number);
						break;
               case 3: // Ellipsoid
                  create::internal::ellipsoid(particle_origin,catom_array,particle_number);
                  break;
					case 4: // Sphere
						create::internal::sphere(particle_origin,catom_array,particle_number);
						break;
					case 5: // Truncated Octahedron
						create::internal::truncated_octahedron(particle_origin,catom_array,particle_number);
						break;
					case 6: // Teardrop
						create::internal::teardrop(particle_origin,catom_array,particle_number);
						break;
               case 7: // Faceted particle
                  create::internal::faceted(particle_origin,catom_array,particle_number);
                  break;
		         case 8: // Cone
			         create::internal::cone(particle_origin,catom_array,particle_number);
			         break;
               case 9: // Bubble
                  create::internal::bubble(particle_origin,catom_array,particle_number);
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

	// Check for no generated particles and print error message
	if(particle_number==0){
		zlog << zTs() << "Error: no particles generated in particle array." << std::endl;
		zlog << zTs() << "Info: Particle arrays require that at least 1 complete particle fits within the system dimensions." << std::endl;
		zlog << zTs() << "Info: Increase x and y system dimensions to at least one particle-scale." << std::endl;
	}

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

void clear_atoms(std::vector<cs::catom_t> & catom_array){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "cs::clear_atoms has been called" << std::endl;}

   // Get original and new number of atoms
   const int num_atoms=catom_array.size();
   int num_included=0;
   for(int a=0;a<num_atoms;a++){
      if(catom_array[a].include == true && mp::material[catom_array[a].material].non_magnetic != 1){
         num_included++;
      }
   }

   // check if there are unneeded atoms
   if(num_atoms!=num_included){
      // create temporary copy for atoms
      std::vector<cs::catom_t> tmp_catom_array(num_atoms);
      tmp_catom_array=catom_array;
      // resize original array to new number of atoms
      catom_array.resize(num_included);
      int atom=0;
      // loop over all existing atoms
      for(int a=0;a<num_atoms;a++){
         // if atom is to be included and is non-magnetic copy to new array
         if(catom_array[a].include==true && mp::material[catom_array[a].material].non_magnetic != 1 ){
            catom_array[atom]=tmp_catom_array[a];
            atom++;
         }
         // if atom is part of a non-magnetic material to be removed then save to nm array
         else if(catom_array[a].include == true && mp::material[catom_array[a].material].non_magnetic == 1){
            cs::nm_atom_t tmp;
         	tmp.x = catom_array[a].x;
         	tmp.y = catom_array[a].y;
         	tmp.z = catom_array[a].z;
         	tmp.mat = catom_array[a].material;
            tmp.cat = catom_array[atom].lh_category;
         	// save atom to non-magnet array
         	cs::non_magnetic_atoms_array.push_back(tmp);
         }
      }
      tmp_catom_array.resize(0);

      zlog << zTs() << "Removed " << cs::non_magnetic_atoms_array.size() << " non-magnetic atoms from system" << std::endl;

   }

   return;

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

void calculate_atomic_composition(std::vector<cs::catom_t> & catom_array){

	zlog<< zTs() << "Determining atomic composition" << std::endl;

	// Determine number of atoms of each class and output to log
	std::vector<unsigned int> MaterialNumbers(mp::num_materials,0);
	for(unsigned int atom=0;atom<catom_array.size();atom++) MaterialNumbers.at(catom_array[atom].material)++;

	// Output composition to log file
	for(int mat=0;mat<mp::num_materials;mat++) zlog << zTs() << "Material " << mat+1 << " " << mp::material[mat].name << " makes up " << double(MaterialNumbers[mat])*100.0/double(catom_array.size()) << "% of all atoms." << std::endl;

	return;

}

int intermixing(std::vector<cs::catom_t> & catom_array){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::intermixing has been called" << std::endl;}

   // re-seed random number generator on each CPU with a different number
   uint64_t long_seed = create::internal::mixing_seed + vmpi::my_rank * vmpi::num_processors;
   uint32_t short_seed = long_seed;   // truncate seed to 32 bit integer with wrap around
	create::internal::grnd.seed(short_seed);

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
				double min = create::internal::mp[mat].min*cs::system_dimensions[2];
				double max = create::internal::mp[mat].max*cs::system_dimensions[2];
				double mean = (min+max)/2.0;
				if(z<=min){
					double probability=0.5+0.5*tanh((z-min)/(mp::material[current_material].intermixing[mat]*cs::system_dimensions[2]));
					if(create::internal::grnd() < probability) final_material=mat;
				}
				else if(z>min && z<=mean){
					double probability=0.5+0.5*tanh((z-min)/(mp::material[current_material].intermixing[mat]*cs::system_dimensions[2]));
					if(create::internal::grnd() < probability) final_material=mat;
				}
				else if(z>mean && z<=max){
					double probability=0.5-0.5*tanh((z-max)/(mp::material[current_material].intermixing[mat]*cs::system_dimensions[2]));
					if(create::internal::grnd() < probability) final_material=mat;
				}
				else if(z>max){
					double probability=0.5-0.5*tanh((z-max)/(mp::material[current_material].intermixing[mat]*cs::system_dimensions[2]));
					//std::cout << current_material << "\t" << mat << "\t" << atom << "\t" << z << "\t" << max << "\t" << probability << std::endl;
					if(create::internal::grnd() < probability) final_material=mat;
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

   // re-seed random number generator on each CPU with a different number
   uint64_t long_seed = create::internal::dilute_seed + vmpi::my_rank * vmpi::num_processors;
   uint32_t short_seed = long_seed;   // truncate seed to 32 bit integer with wrap around
   create::internal::grnd.seed(short_seed);

   // loop over all atoms
   for(unsigned int atom=0;atom<catom_array.size();atom++){
      // if atom material is alloy master
      int local_material=catom_array[atom].material;
      double probability = mp::material[local_material].density;
      if(create::internal::grnd() > probability) catom_array[atom].include=false;
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

		// Re-identify all atoms as material 0 and exclude by default
		for(unsigned int atom=0;atom<catom_array.size();atom++){
			catom_array[atom].material=0;
			catom_array[atom].include=false;
		}

		// loop over all materials and include according to geometry

		// determine z-bounds for materials
		std::vector<double> mat_min(mp::num_materials);
		std::vector<double> mat_max(mp::num_materials);

		// initialise array to hold z-bound information in real space
		for(int mat=0;mat<mp::num_materials;mat++){
			mat_min[mat]=create::internal::mp[mat].min*cs::system_dimensions[2];
			mat_max[mat]=create::internal::mp[mat].max*cs::system_dimensions[2];
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
					const double z = catom_array[atom].z;
					if((z>=mat_min[mat]) && (z<mat_max[mat]) && (vmath::point_in_polygon2(x,y,px,py,geo)==true)){
						catom_array[atom].material=mat;
						catom_array[atom].include=true;
					}
				}
			}
		}
	}

	return;
}

//-----------------------------------------------------------
//
///  Function to replace deleted atoms with in-fill material
//
///  Can be used to create embedded nanoparticle arrays or
///  granular recording media.
//
///  v1 18/09/2013
///  (c) R F L Evans
//
//-----------------------------------------------------------
void fill(std::vector<cs::catom_t> & catom_array){
   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "cs::fill has been called" << std::endl;}

   //loop over all potential intermixing materials
   for(int mat=0;mat<mp::num_materials;mat++){
      if(mp::material[mat].fill){
         double min = create::internal::mp[mat].min*cs::system_dimensions[2];
         double max = create::internal::mp[mat].max*cs::system_dimensions[2];

         // loop over all atoms selecting only deselected atoms within min/max
         for(unsigned int atom=0;atom<catom_array.size();atom++){
            if( (catom_array[atom].z < max) && (catom_array[atom].z >= min) && (catom_array[atom].include==false)){
               // set atom to fill material
               catom_array[atom].material=mat;
               // re-include atom
               catom_array[atom].include=true;
            }
         }
      }
   }

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to alter particle origin to be centred on an atom
//------------------------------------------------------------------------------------------------------
void centre_particle_on_atom(std::vector<double>& particle_origin, std::vector<cs::catom_t>& catom_array){

   vmpi::barrier();

   // set initial max range
   double max_range_sq = 1e123;
   unsigned int nearest; // nearest atom to initial particle origin

   // copy to temporary for speed
   const double prx = particle_origin[0];
   const double pry = particle_origin[1];
   const double prz = particle_origin[2];

   // loop over all atoms to find closest atom
   for(int atom=0;atom<catom_array.size();atom++){
      double dx = catom_array[atom].x-particle_origin[0];
      double dy = catom_array[atom].y-particle_origin[1];
      double dz = catom_array[atom].z-particle_origin[2];
      double r = dx*dx + dy*dy + dz*dz;
      if(r < max_range_sq){
         max_range_sq = r;
         nearest = atom;
      }
   }

   // set particle origin to nearest atom
   particle_origin[0] = catom_array[nearest].x;
   particle_origin[1] = catom_array[nearest].y;
   particle_origin[2] = catom_array[nearest].z;

   //-----------------------------------------------------
   // For parallel reduce on all CPUs
   //-----------------------------------------------------
   #ifdef MPICF

      // set up array to get ranges from all CPUs on rank 0
      std::vector<double> ranges;
      if(vmpi::my_rank == 0) ranges.resize(vmpi::num_processors, 1.e123);
      else ranges.resize(1,0.0); // one value sufficient on all other CPUs

      // gather max ranges from all cpus on root (1 data point from each process)
      MPI_Gather(&max_range_sq, 1, MPI_DOUBLE, &ranges[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      // variable to store rank of minimum range
      unsigned int rank_of_min_range=0;

      // work out minimum range on root
      if(vmpi::my_rank==0){

         double min_range = 1.e123;

         // loop over all ranges and determine minimum and cpu location
         for(int i=0; i<ranges.size(); i++){
            if(ranges[i] < min_range){
               min_range = ranges[i];
               rank_of_min_range = i;
            }
         }
      }

      // broadcast id of nearest to all cpus from root
      MPI_Bcast(&rank_of_min_range, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

      // broadcast position to all cpus
      MPI_Bcast(&particle_origin[0], 3, MPI_DOUBLE, rank_of_min_range, MPI_COMM_WORLD);

      vmpi::barrier();

   #endif

   return;

}

} // end of namespace
