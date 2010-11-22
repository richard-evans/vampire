///
/// @file
/// @brief Contains functions for cutting shapes from crystals. 
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    05/03/2010
/// @internal
///	Created:		05/03/2010
///	Revision:	  ---
///=====================================================================================
///

#include "errors.hpp"
#include "create.hpp"
#include "material.hpp"
#include "vmath.hpp"
#include <cmath>
#include <string>
#include <iostream>

namespace cs{

int bulk(std::vector<cs::catom_t> & catom_array, const int grain){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::bulk has been called" << std::endl;}

	const int max_vertices=50;

	//----------------------------------------------------
	// Loop over all atoms and mark atoms within geometry
	//----------------------------------------------------
	const int num_atoms = catom_array.size();
	
 	for(int atom=0;atom<num_atoms;atom++){

		const int geo=mp::material[catom_array[atom].material].geometry;

		if(geo==0){
			catom_array[atom].include=true;
		}
		else{
			double x = catom_array[atom].x;
			double y = catom_array[atom].y;
			double px[max_vertices];
			double py[max_vertices];
			// Initialise polygon points
			for(int p=0;p<geo;p++){
				px[p]=mp::material[catom_array[atom].material].geometry_coords[p][0]*mp::system_dimensions[0];
				py[p]=mp::material[catom_array[atom].material].geometry_coords[p][1]*mp::system_dimensions[1];
			}
			if(vmath::point_in_polygon(x,y,px,py,geo)==true){
				catom_array[atom].include=true;
				catom_array[atom].grain=grain;
			}
		}
	}
	return EXIT_SUCCESS;	
}

int cylinder(double particle_origin[],std::vector<cs::catom_t> & catom_array, const int grain){
	
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::cylinder has been called" << std::endl;}

	//-----------------------------------------
	// Set particle radius
	//-----------------------------------------
	double particle_radius_squared = (mp::particle_scale*0.5)*(mp::particle_scale*0.5);
	
	//-----------------------------------------------
	// Loop over all atoms and mark atoms in sphere
	//-----------------------------------------------
	const int num_atoms = catom_array.size();
	
 	for(int atom=0;atom<num_atoms;atom++){
		double range_squared = 	(catom_array[atom].x-particle_origin[0])*(catom_array[atom].x-particle_origin[0]) + 
										(catom_array[atom].y-particle_origin[1])*(catom_array[atom].y-particle_origin[1]);
		if(range_squared<=particle_radius_squared){
			catom_array[atom].include=true;
			catom_array[atom].grain=grain;
		}
	}
	return EXIT_SUCCESS;	
}

int sphere(double particle_origin[],std::vector<cs::catom_t> & catom_array, const int grain){
	//====================================================================================
	//
	//									cs_sphere
	//
	//					Subroutine to cut a spherical particle shape
	//
	//							Version 1.0 R Evans 22/09/2008
	//
	//====================================================================================
	
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::sphere has been called" << std::endl;}

	// Set particle radius
	double particle_radius_squared = (mp::particle_scale*0.5)*(mp::particle_scale*0.5);
	
	// Loop over all atoms and mark atoms in sphere
	const int num_atoms = catom_array.size();
	
 	for(int atom=0;atom<num_atoms;atom++){
		double range_squared = (catom_array[atom].x-particle_origin[0])*(catom_array[atom].x-particle_origin[0]) + 
							 (catom_array[atom].y-particle_origin[1])*(catom_array[atom].y-particle_origin[1]) +
							 (catom_array[atom].z-particle_origin[2])*(catom_array[atom].z-particle_origin[2]);
		if(mp::material[catom_array[atom].material].core_shell_size>0.0){
			// Reverse Loop over materials
			for(int mat=mp::num_materials-1;mat>-1;mat--){
				double my_radius = mp::material[mat].core_shell_size;
				double max_range = my_radius*my_radius*particle_radius_squared;
				if(range_squared<=max_range){
					catom_array[atom].include=true;
					catom_array[atom].material=mat;
					catom_array[atom].grain=grain;
				}
			}
		}
		else if(range_squared<=particle_radius_squared) {
			catom_array[atom].include=true;
			catom_array[atom].grain=grain;
		}
	}
	return EXIT_SUCCESS;	
}

int truncated_octahedron(double particle_origin[],std::vector<cs::catom_t> & catom_array, const int grain){
	//====================================================================================
	//
	//								cs_truncated_octahedron
	//
	//				Subroutine to cut a truncated octahedron particle shape
	//
	//							Version 1.0 R Evans 22/09/2008
	//
	//====================================================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::truncated_octahedron has been called" << std::endl;}

	// Set truncated octahedron parameters
	const double to_length = mp::particle_scale*0.5*3.0/2.0;
	const double to_height = mp::particle_scale*0.5;
	double x_vector[3];
	
	// Loop over all atoms and mark atoms in truncate octahedron
	const int num_atoms = catom_array.size();
				  	
	for(int atom=0;atom<num_atoms;atom++){
		x_vector[0] = fabs(catom_array[atom].x-particle_origin[0]);
		x_vector[1] = fabs(catom_array[atom].y-particle_origin[1]);
		x_vector[2] = fabs(catom_array[atom].z-particle_origin[2]);
		double range = x_vector[0] + x_vector[1] + x_vector[2];

		if(mp::material[catom_array[atom].material].core_shell_size>0.0){
			// Reverse Loop over materials
			for(int mat=mp::num_materials-1;mat>-1;mat--){
				double my_radius = mp::material[mat].core_shell_size;
				double my_to_height = my_radius*to_height;
				double my_to_length = my_radius*to_length;
				
				if((range<=my_to_length) && (x_vector[0] <= my_to_height) && (x_vector[1] <= my_to_height) && (x_vector[2] <= my_to_height)){
					catom_array[atom].include=true;
					catom_array[atom].material=mat;
					catom_array[atom].grain=grain;
				}
			}
		}
		else if((range<=to_length) && (x_vector[0] <= to_height) && (x_vector[1] <= to_height) && (x_vector[2] <= to_height)){
			catom_array[atom].include=true;
			catom_array[atom].grain=grain;
		}
	}
	
	return EXIT_SUCCESS;	
}

}
