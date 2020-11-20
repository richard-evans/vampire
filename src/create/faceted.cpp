//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <list>

// Vampire headers
#include "create.hpp"

// Internal sim header
#include "internal.hpp"


namespace create{
namespace internal{

void faceted(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain){
	//====================================================================================
	//
	//	Function to cut a wulff particle shape
	//
	// (c) R F L Evans (2016). All rights reserved.
	//
	//====================================================================================

	// Set faceted particle parameters
   const double sfx = cs::particle_shape_factor_x;
   const double sfy = cs::particle_shape_factor_y;
   const double sfz = cs::particle_shape_factor_z;
	const double rsize = cs::particle_scale*0.5; // reduced particle size
   const double fr100 = create::internal::faceted_particle_100_radius; // local constant for 100 facet radius
   const double fr110 = create::internal::faceted_particle_110_radius; // local constant for 110 facet radius
   const double fr111 = create::internal::faceted_particle_111_radius; // local constant for 111 facet radius
   const double root2 = sqrt(2.0);

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

      // calculate reduced atom position
		double rx = fabs(catom_array[atom].x-particle_origin[0]);
		double ry = fabs(catom_array[atom].y-particle_origin[1]);
		double rz = fabs(catom_array[atom].z-particle_origin[2]);

		if(cs::core_shell_particles){
         // Iterate over materials
         for(std::list<core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
            int mat = (it)->mat;
				const double my_rsize = rsize*mp::material[mat].core_shell_size;;

            // Calculate facet radii
            bool in100 = (rx <= fr100*my_rsize*sfx) && (ry <= fr100*my_rsize*sfy) && (rz <= fr100*my_rsize*sfz);
            bool in110 = (rx <= fr110*my_rsize*root2 + fr110*my_rsize*(sfx-1.0) + fr110*my_rsize*(sfy-1.0) - ry) &&
                         (rz <= fr110*my_rsize*root2 + fr110*my_rsize*(sfz-1.0) + fr110*my_rsize*(sfy-1.0) - ry) &&
                         (rz <= fr110*my_rsize*root2 + fr110*my_rsize*(sfz-1.0) + fr110*my_rsize*(sfx-1.0) - rx);
            bool in111 = rx + ry + rz < 1.5*my_rsize*fr111 + fr111*my_rsize*((sfx-1.0) + (sfy-1.0) + (sfz-1.0));

            if(in100 && in110 && in111){
               double maxz=create::internal::mp[mat].max*cs::system_dimensions[2];
               double minz=create::internal::mp[mat].min*cs::system_dimensions[2];
               double cz=catom_array[atom].z;
               const int atom_uc_cat = catom_array[atom].uc_category;
               const int mat_uc_cat = create::internal::mp[mat].unit_cell_category;
					if((cz>=minz) && (cz<maxz) && (atom_uc_cat == mat_uc_cat) ){
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
		else
      {

         // Calculate facet radii
         bool in100 = (rx <= fr100*rsize*sfx) && (ry <= fr100*rsize*sfy) && (rz <= fr100*rsize*sfz);
         bool in110 = (rx <= fr110*rsize*root2 + fr110*rsize*(sfx-1.0) + fr110*rsize*(sfy-1.0) - ry) &&
                      (rz <= fr110*rsize*root2 + fr110*rsize*(sfz-1.0) + fr110*rsize*(sfy-1.0) - ry) &&
                      (rz <= fr110*rsize*root2 + fr110*rsize*(sfz-1.0) + fr110*rsize*(sfx-1.0) - rx);
         bool in111 = rx + ry + rz < 1.5*rsize*fr111 + fr111*rsize*((sfx-1.0) + (sfy-1.0) + (sfz-1.0));

         if(in100 && in110 && in111){
            catom_array[atom].include=true;
            catom_array[atom].grain=grain;
		   }
      }
	}

	return;
}

} // end of namespace internal
} // end of create namespace
