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
#include <cstdlib>
#include <list>
#include <string>
#include <iostream>

struct core_radius_t{
   int mat;
   double radius;
};

/// comparison function for reverse order sorting
bool compare_radius(core_radius_t first,core_radius_t second){
   if(first.radius<second.radius) return false;
   else return true;
}

namespace cs{

int bulk(std::vector<cs::catom_t> & catom_array){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::bulk has been called" << std::endl;}

	// Loop over all atoms and mark as selected
	const int num_atoms = catom_array.size();
	
 	for(int atom=0;atom<num_atoms;atom++){
		catom_array[atom].include=true;
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
	double particle_radius_squared = (cs::particle_scale*0.5)*(cs::particle_scale*0.5);
	
	//-----------------------------------------------
	// Loop over all atoms and mark atoms in sphere
	//-----------------------------------------------
	const int num_atoms = catom_array.size();

   // determine order for core-shell particles
   std::list<core_radius_t> material_order(0);
   for(int mat=0;mat<mp::num_materials;mat++){
      core_radius_t tmp;
      tmp.mat=mat;
      tmp.radius=mp::material[mat].core_shell_size;
      material_order.push_back(tmp);
   }
   // sort by increasing radius
   material_order.sort(compare_radius);

 	for(int atom=0;atom<num_atoms;atom++){
		double range_squared = 	(catom_array[atom].x-particle_origin[0])*(catom_array[atom].x-particle_origin[0]) + 
										(catom_array[atom].y-particle_origin[1])*(catom_array[atom].y-particle_origin[1]);
		if(mp::material[catom_array[atom].material].core_shell_size>0.0){
         // Iterate over materials
         for(std::list<core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
            int mat = (it)->mat;
				double my_radius = mp::material[mat].core_shell_size;
				double max_range = my_radius*my_radius*particle_radius_squared;
				double maxz=mp::material[mat].max*cs::system_dimensions[2];
				double minz=mp::material[mat].min*cs::system_dimensions[2];
				double cz=catom_array[atom].z;
				// check for within core shell range
				if(range_squared<=max_range){
					if((cz>=minz) && (cz<maxz)){
					catom_array[atom].include=true;
					catom_array[atom].material=mat;
					catom_array[atom].grain=grain;
					}
					// if set to clear atoms then remove atoms within radius
					else if(cs::fill_core_shell==false){
						catom_array[atom].include=false;
					}
				}
			}
		}
		else if(range_squared<=particle_radius_squared){
			catom_array[atom].include=true;
			catom_array[atom].grain=grain;
		}
	}
	return EXIT_SUCCESS;	
}

void ellipsoid(double particle_origin[],std::vector<cs::catom_t> & catom_array, const int grain){
   //--------------------------------------------------------------------------------------------
   //
   ///  Function to cut an ellipsoid particle shape
   ///
   ///  (C) R F L Evans 23/04/2013
   ///
   ///  Equation for an ellipsoid is:
   ///
   ///  (x-h)^2   (y-k)^2    (z-l)^2
   ///  -------   -------    -------  <= 1
   ///   r_x^2     r_y^2      r_z^2
   ///
   ///  Evaluation of the above for particle coordinates h,k,l
   ///  and atomic positions x,y,z will include all atoms within
   ///  ellipse with vertices r_x, r_y and r_z.
   //
   //--------------------------------------------------------------------------------------------

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "cs::ellipsoid has been called" << std::endl;}

   // Set particle radius
   const double particle_radius_squared = (cs::particle_scale*0.5)*(cs::particle_scale*0.5);

   // Use shape modifiers to generate ellipsoid
   const double inv_rx_sq = 1.0/(particle_radius_squared*cs::particle_shape_factor_x*cs::particle_shape_factor_x);
   const double inv_ry_sq = 1.0/(particle_radius_squared*cs::particle_shape_factor_y*cs::particle_shape_factor_y);
   const double inv_rz_sq = 1.0/(particle_radius_squared*cs::particle_shape_factor_z*cs::particle_shape_factor_z);

   // Loop over all atoms and mark atoms in sphere
   const int num_atoms = catom_array.size();

   // determine order for core-shell particles
   std::list<core_radius_t> material_order(0);
   for(int mat=0;mat<mp::num_materials;mat++){
      core_radius_t tmp;
      tmp.mat=mat;
      tmp.radius=mp::material[mat].core_shell_size;
      material_order.push_back(tmp);
   }
   // sort by increasing radius
   material_order.sort(compare_radius);

   for(int atom=0;atom<num_atoms;atom++){
      const double range_x_sq = (catom_array[atom].x-particle_origin[0])*(catom_array[atom].x-particle_origin[0]);
      const double range_y_sq = (catom_array[atom].y-particle_origin[1])*(catom_array[atom].y-particle_origin[1]);
      const double range_z_sq = (catom_array[atom].z-particle_origin[2])*(catom_array[atom].z-particle_origin[2]);

      if(mp::material[catom_array[atom].material].core_shell_size>0.0){
         // Iterate over materials
         for(std::list<core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
            int mat = (it)->mat;
            const double my_radius_sq = mp::material[mat].core_shell_size*mp::material[mat].core_shell_size;
            double maxz=mp::material[mat].max*cs::system_dimensions[2];
            double minz=mp::material[mat].min*cs::system_dimensions[2];
            double cz=catom_array[atom].z;
            // check for within core shell range
            if(range_x_sq*inv_rx_sq + range_y_sq*inv_ry_sq +range_z_sq*inv_rz_sq<=my_radius_sq){
               if((cz>=minz) && (cz<maxz)){
                  catom_array[atom].include=true;
                  catom_array[atom].material=mat;
                  catom_array[atom].grain=grain;
               }
               // if set to clear atoms then remove atoms within radius
               else if(cs::fill_core_shell==false){
                  catom_array[atom].include=false;
               }
            }
         }
      }
      else if(range_x_sq*inv_rx_sq + range_y_sq*inv_ry_sq +range_z_sq*inv_rz_sq<=1.0) {
         catom_array[atom].include=true;
         catom_array[atom].grain=grain;
      }
   }
   return;
}

int sphere(double particle_origin[],std::vector<cs::catom_t> & catom_array, const int grain){
	//====================================================================================
	///
	///									cs_sphere
	///
	///					Subroutine to cut a spherical particle shape
	///
	///							Version 1.0 R Evans 22/09/2008
	//
	//====================================================================================
	
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::sphere has been called" << std::endl;}

	// Set particle radius
	double particle_radius_squared = (cs::particle_scale*0.5)*(cs::particle_scale*0.5);
	
	// Loop over all atoms and mark atoms in sphere
	const int num_atoms = catom_array.size();

   // determine order for core-shell particles
   std::list<core_radius_t> material_order(0);
   for(int mat=0;mat<mp::num_materials;mat++){
      core_radius_t tmp;
      tmp.mat=mat;
      tmp.radius=mp::material[mat].core_shell_size;
      material_order.push_back(tmp);
   }
   // sort by increasing radius
   material_order.sort(compare_radius);

 	for(int atom=0;atom<num_atoms;atom++){
		double range_squared = (catom_array[atom].x-particle_origin[0])*(catom_array[atom].x-particle_origin[0]) + 
							 (catom_array[atom].y-particle_origin[1])*(catom_array[atom].y-particle_origin[1]) +
							 (catom_array[atom].z-particle_origin[2])*(catom_array[atom].z-particle_origin[2]);
		if(mp::material[catom_array[atom].material].core_shell_size>0.0){
         // Iterate over materials
         for(std::list<core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
            int mat = (it)->mat;
				double my_radius = mp::material[mat].core_shell_size;
				double max_range = my_radius*my_radius*particle_radius_squared;
				double maxz=mp::material[mat].max*cs::system_dimensions[2];
				double minz=mp::material[mat].min*cs::system_dimensions[2];
				double cz=catom_array[atom].z;
				// check for within core shell range
				if(range_squared<=max_range){
					if((cz>=minz) && (cz<maxz)){
						catom_array[atom].include=true;
						catom_array[atom].material=mat;
						catom_array[atom].grain=grain;
					}
					// if set to clear atoms then remove atoms within radius
					else if(cs::fill_core_shell==false){
						catom_array[atom].include=false;
					}
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
	///								cs_truncated_octahedron
	///
	///				Subroutine to cut a truncated octahedron particle shape
	///
	///							Version 1.0 R Evans 22/09/2008
	///
	//====================================================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::truncated_octahedron has been called" << std::endl;}

	// Set truncated octahedron parameters
	const double to_length = cs::particle_scale*0.5*3.0/2.0;
	const double to_height = cs::particle_scale*0.5;
        //const double to_length = cs::particle_scale*0.5;
	//const double to_height = to_length*2.0/3.0;
	double x_vector[3];
	
	// Loop over all atoms and mark atoms in truncate octahedron
	const int num_atoms = catom_array.size();

   // determine order for core-shell particles
   std::list<core_radius_t> material_order(0);
   for(int mat=0;mat<mp::num_materials;mat++){
      core_radius_t tmp;
      tmp.mat=mat;
      tmp.radius=mp::material[mat].core_shell_size;
      material_order.push_back(tmp);
   }
   // sort by increasing radius
   material_order.sort(compare_radius);

	for(int atom=0;atom<num_atoms;atom++){
		x_vector[0] = fabs(catom_array[atom].x-particle_origin[0]);
		x_vector[1] = fabs(catom_array[atom].y-particle_origin[1]);
		x_vector[2] = fabs(catom_array[atom].z-particle_origin[2]);
		double range = x_vector[0] + x_vector[1] + x_vector[2];

		if(mp::material[catom_array[atom].material].core_shell_size>0.0){
         // Iterate over materials
         for(std::list<core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
            int mat = (it)->mat;
				double my_radius = mp::material[mat].core_shell_size;
				double my_to_height = my_radius*to_height;
				double my_to_length = my_radius*to_length;
				double maxz=mp::material[mat].max*cs::system_dimensions[2];
				double minz=mp::material[mat].min*cs::system_dimensions[2];
				double cz=catom_array[atom].z;
				// check for within core shell range
				if((range<=my_to_length) && (x_vector[0] <= my_to_height) && (x_vector[1] <= my_to_height) && (x_vector[2] <= my_to_height)){
					if((cz>=minz) && (cz<maxz)){
						catom_array[atom].include=true;
						catom_array[atom].material=mat;
						catom_array[atom].grain=grain;
					}
					// if set to clear atoms then remove atoms within radius
					else if(cs::fill_core_shell==false){
						catom_array[atom].include=false;
					}
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

int cube(double particle_origin[],std::vector<cs::catom_t> & catom_array, const int grain){
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::cube has been called" << std::endl;}

	// Set particle size
	double side_length=cs::particle_scale*0.5;

	// Loop over all atoms and mark atoms in cube
	const int num_atoms = catom_array.size();

   // determine order for core-shell particles
   std::list<core_radius_t> material_order(0);
   for(int mat=0;mat<mp::num_materials;mat++){
      core_radius_t tmp;
      tmp.mat=mat;
      tmp.radius=mp::material[mat].core_shell_size;
      material_order.push_back(tmp);
   }
   // sort by increasing radius
   material_order.sort(compare_radius);

 	for(int atom=0;atom<num_atoms;atom++){
		double dx=fabs(catom_array[atom].x-particle_origin[0]);
		double dy=fabs(catom_array[atom].y-particle_origin[1]);
		if(mp::material[catom_array[atom].material].core_shell_size>0.0){
         // Iterate over materials
         for(std::list<core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
            int mat = (it)->mat;
				double my_length = mp::material[mat].core_shell_size;
				double max_length = my_length*side_length;
				double maxz=mp::material[mat].max*cs::system_dimensions[2];
				double minz=mp::material[mat].min*cs::system_dimensions[2];
				double cz=catom_array[atom].z;
				// check for within core shell range
				if((dx<=max_length) && (dy<=max_length)){
					if((cz>=minz) && (cz<maxz)){
						catom_array[atom].include=true;
						catom_array[atom].material=mat;
						catom_array[atom].grain=grain;
					}
					// if set to clear atoms then remove atoms within radius
					else if(cs::fill_core_shell==false){
						catom_array[atom].include=false;
					}
				}
			}
		}
		else if((dx<=side_length) && (dy<=side_length)){
			catom_array[atom].include=true;
			catom_array[atom].grain=grain;
		}
	}
	return EXIT_SUCCESS;	

}

// Teardrop
int tear_drop(double particle_origin[],std::vector<cs::catom_t> & catom_array, const int grain){
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::tear_drop has been called" << std::endl;}

	// teapdrop dimensions
	// 0.01242725414 = 6nm/(6nm + 6nm +500nm)
	double TeardropMinZ=0.01242725414; // Frac system height
	double TeardropMaxZ=0.01242725414+0.01242725414;
	double TeardropRadius=TeardropMaxZ-TeardropMinZ;
	double TeardropMinRadius=1.5; // Angstroms
	
	// Set particle size
	double side_length=cs::particle_scale*0.5;

	// Loop over all atoms and mark atoms in cube
	const int num_atoms = catom_array.size();
	
 	for(int atom=0;atom<num_atoms;atom++){
		double dx=fabs(catom_array[atom].x-particle_origin[0]);
		double dy=fabs(catom_array[atom].y-particle_origin[1]);
		
		// check for atoms constrained by box
		if((dx<=side_length) && (dy<=side_length)){
			
			// // check for lower box
			if(catom_array[atom].z <= cs::system_dimensions[2]*TeardropMinZ){
			catom_array[atom].include=true;
			catom_array[atom].grain=grain;
			}
			else if(catom_array[atom].z >= cs::system_dimensions[2]*(1.0-TeardropMinZ)){
			catom_array[atom].include=true;
			catom_array[atom].grain=grain;
			}
			// check for teardrop part
			else{
				double Height;
				// z < 0.5
				if(catom_array[atom].z <= cs::system_dimensions[2]*0.5){
					Height=catom_array[atom].z-cs::system_dimensions[2]*TeardropMinZ;
				}
				else{
					Height=cs::system_dimensions[2]*(1.0-TeardropMinZ)-catom_array[atom].z;
				}
				double RadiusAtHeight=cs::particle_scale*0.5*exp(-Height/(TeardropRadius*cs::system_dimensions[2]))+TeardropMinRadius;
				double RadiusSquared=dx*dx+dy*dy;
				if(RadiusSquared<=RadiusAtHeight*RadiusAtHeight){
					catom_array[atom].include=true;
					catom_array[atom].grain=grain;
				}

			}
			
		}
	}
	return EXIT_SUCCESS;	

}

}
