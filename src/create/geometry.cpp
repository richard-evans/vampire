//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "create.hpp"
#include "vio.hpp"
#include "vmath.hpp"


// create module headers
#include "internal.hpp"

//------------------------------------------------------------------------------
//
// Function to implement a custom 2D geometry specified by the user as a set of
// points. For example, a square as follows:
//
//                  *------------------------*
//                  |                        |
//                  |      *          *      |
//                  |                        |
//                  |                        |
//                  |                        |
//                  |                        |
//                  |      *          *      |
//                  |                        |
//                  *------------------------*
//
//------------------------------------------------------------------------------
void create::internal::geometry(std::vector<cs::catom_t>& catom_array){

	//-----------------------------------------------
   // Check for any geometry in material parameters
   //-----------------------------------------------
	bool cut=false;
	for(int mat=0; mat<mp::num_materials; mat++){
		if(mp::material[mat].geometry>0) cut=true;
	}

	// Return from function if no geometry is defined.
	if( cut == false ) return;

   //-----------------------------------------------
   // Otherwise proceed
   //-----------------------------------------------
	zlog << zTs() << "Cutting materials within defined geometry." << std::endl;

   //-----------------------------------------------
   // Check for force material type by geometry
   //-----------------------------------------------
   if( !create::internal::select_material_by_geometry ){

		// loop over all atoms
      // c++11 style range loop with auto reference to class variable
      // for(unsigned int atom=0; atom < catom_array.size(); atom++){
      for(auto& atom : catom_array){

			// check for geometry information
			const int geo = mp::material[atom.material].geometry;

			// if exists, then remove atoms outside polygon
			if( geo > 0 ){
				double x = atom.x; // store polygon points in material class as tuples... TBD
				double y = atom.y;
				std::vector<double> px(geo);
				std::vector<double> py(geo);
				// Initialise polygon points
				for(int p=0;p<geo;p++){
					px[p]=mp::material[atom.material].geometry_coords[p][0]*cs::system_dimensions[0];
					py[p]=mp::material[atom.material].geometry_coords[p][1]*cs::system_dimensions[1];
				}
				// check if point is outside of polygon, if so delete it
				if( !vmath::point_in_polygon2(x,y,px,py,geo) ) atom.include = false;

			}

		}

	}
   //--------------------------------------------------
   // Otherwise use special algorithm to include atoms
   //-------------------------------------------------
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
