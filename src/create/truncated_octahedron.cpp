//-----------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <list>

// Vampire headers
#include "create.hpp"

// Internal create header
#include "internal.hpp"

namespace create{

namespace internal{

   void truncated_octahedron(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain){

   	// Set truncated octahedron parameters
      const double sfx = cs::particle_shape_factor_x;
      const double sfy = cs::particle_shape_factor_y;
      const double sfz = cs::particle_shape_factor_z;
   	const double to_length = cs::particle_scale*0.5*3.0/2.0;
   	const double to_height = cs::particle_scale*0.5;

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

      // Unroll min, max and fill for performance
      std::vector<double> mat_min(mp::num_materials);
      std::vector<double> mat_max(mp::num_materials);
      std::vector<double> mat_cssize(mp::num_materials);
      std::vector<int> uc_cat(mp::num_materials); // array of material -> unit cell material associations

      for(int mat=0;mat<mp::num_materials;mat++){
         mat_min[mat]=create::internal::mp[mat].min*cs::system_dimensions[2];
         mat_max[mat]=create::internal::mp[mat].max*cs::system_dimensions[2];
         mat_cssize[mat] = mp::material[mat].core_shell_size;
         // alloys generally are not defined by height, and so have max = 0.0
         if(mat_max[mat]<1.e-99) mat_max[mat]=-0.1;
         uc_cat[mat] = create::internal::mp[mat].unit_cell_category; // unit cell category of material
      }

   	for(int atom=0;atom<num_atoms;atom++){

   		x_vector[0] = fabs(catom_array[atom].x - particle_origin[0]);
   		x_vector[1] = fabs(catom_array[atom].y - particle_origin[1]);
   		x_vector[2] = fabs(catom_array[atom].z - particle_origin[2]);

         const double cz = catom_array[atom].z;
         const int atom_uc_cat = catom_array[atom].uc_category;

   		double range = x_vector[0] + x_vector[1] + x_vector[2];

   		if(mp::material[catom_array[atom].material].core_shell_size>0.0){
            // Iterate over materials
            for(std::list<core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
               int mat = (it)->mat;
   				double my_radius = mat_cssize[mat];
   				double my_to_height = my_radius*to_height;
   				double my_to_length = my_radius*to_length;
   				double maxz = mat_max[mat];
   				double minz = mat_min[mat];
   				// check for within core shell range
   				if((range<=my_to_length) && (x_vector[0] <= my_to_height*sfx) && (x_vector[1] <= my_to_height*sfy) && (x_vector[2] <= my_to_height*sfz)){
   					if( (cz >= minz) && (cz < maxz) && (atom_uc_cat == uc_cat[mat])){
   						catom_array[atom].include = true;
   						catom_array[atom].material = mat;
   						catom_array[atom].grain = grain;
   					}
   					// if set to clear atoms then remove atoms within radius
   					else if(cs::fill_core_shell == false){
   						catom_array[atom].include = false;
   					}
   				}
   			}
   		}
   		else if((range<=to_length) && (x_vector[0] <= to_height*sfx) && (x_vector[1] <= to_height*sfy) && (x_vector[2] <= to_height*sfz)){
   			catom_array[atom].include = true;
   			catom_array[atom].grain = grain;
   		}
   	}

   	return;
   }


} // end of internal namespace

} // end of create namespace
