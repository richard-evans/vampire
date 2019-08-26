//-----------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "create.hpp"

// Internal create header
#include "internal.hpp"

namespace create{

namespace internal{

   //------------------------------------------------------------
   // Function to cretate tear drop shape
   //------------------------------------------------------------
   void teardrop(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain){

   	// teapdrop dimensions
   	// 0.01242725414 = 6nm/(6nm + 6nm +500nm)
   	double teardrop_min_z = 0.01242725414; // Frac system height
   	double teardrop_max_z = 0.01242725414 + 0.01242725414;
   	double teardrop_radius = teardrop_max_z - teardrop_min_z;
   	double teardrop_min_radius = 1.5; // Angstroms

   	// Set particle size
   	double side_length=cs::particle_scale*0.5;

   	// Loop over all atoms and mark atoms in cube
   	const int num_atoms = catom_array.size();

    	for(int atom=0; atom < num_atoms; atom++){
   		double dx=fabs(catom_array[atom].x-particle_origin[0]);
   		double dy=fabs(catom_array[atom].y-particle_origin[1]);

   		// check for atoms constrained by box
   		if((dx<=side_length) && (dy<=side_length)){

   			// // check for lower box
   			if(catom_array[atom].z <= cs::system_dimensions[2]*teardrop_min_z){
   			catom_array[atom].include=true;
   			catom_array[atom].grain=grain;
   			}
   			else if(catom_array[atom].z >= cs::system_dimensions[2]*(1.0-teardrop_min_z)){
   			catom_array[atom].include=true;
   			catom_array[atom].grain=grain;
   			}
   			// check for teardrop part
   			else{
   				double height;
   				// z < 0.5
   				if(catom_array[atom].z <= cs::system_dimensions[2]*0.5){
   					height=catom_array[atom].z-cs::system_dimensions[2]*teardrop_min_z;
   				}
   				else{
   					height=cs::system_dimensions[2]*(1.0-teardrop_min_z)-catom_array[atom].z;
   				}
   				double radius_at_height=cs::particle_scale*0.5*exp(-height/(teardrop_radius*cs::system_dimensions[2]))+teardrop_min_radius;
   				double radius_squared = dx*dx + dy*dy;
   				if(radius_squared <= radius_at_height*radius_at_height){
   					catom_array[atom].include=true;
   					catom_array[atom].grain=grain;
   				}

   			}

   		}
   	}

   	return;

   }

} // end of internal namespace

} // end of create namespace
