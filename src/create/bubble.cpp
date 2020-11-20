//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <list>

// Vampire headers
#include "create.hpp"

// Internal sim header
#include "internal.hpp"

namespace create{
   namespace internal{

      void bubble(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain){

         // Set bubble nucleation height and radius
         const double nucleation_height = create::internal::bubble_nucleation_height;
         const double bubble_radius = create::internal::bubble_radius;

         // calculate reduced ranges for materials with small offset to prevent dangling atoms
         const double rminz = -0.01;
         const double rmaxz = 1.01;
         const double ssz = cs::system_dimensions[2];

         // Determine core shell radii
         std::list<core_radius_t> material_order(0);
         for(int mat=0;mat<mp::num_materials;mat++){
            core_radius_t tmp;
            tmp.mat=mat;
            tmp.radius=mp::material[mat].core_shell_size;
            material_order.push_back(tmp);
         }
         // sort by increasing radius
         material_order.sort(compare_radius);

			//-----------------------------------------------
			// Loop over all atoms and mark atoms in sphere
			//-----------------------------------------------
			const int num_atoms = catom_array.size();

         const double particle_radius_sq = 0.25 * cs::particle_scale * cs::particle_scale;

 	 	 	for(int atom=0;atom<num_atoms;atom++){

            const double range_sq = (catom_array[atom].x-particle_origin[0])*(catom_array[atom].x-particle_origin[0]) +
                                    (catom_array[atom].y-particle_origin[1])*(catom_array[atom].y-particle_origin[1]);

            const double z = catom_array[atom].z;
            const double frh = z/ssz;
            double factor_radius = 0.0;

            if(frh > nucleation_height){
               // multiply by small factor to ensure grains touch at boundary for zero spacing
               factor_radius = 1.04*pow(1.0+((nucleation_height-frh)/(rmaxz-nucleation_height)),bubble_radius);
            }
            else{
               factor_radius = 1.04*pow((1.0-(frh-nucleation_height)/(rminz-nucleation_height)),bubble_radius);
            }

            /*if(cs::core_shell_particles){
               // Iterate over materials
               for(std::list<core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
                  int mat = (it)->mat;
						double my_radius = mp::material[mat].core_shell_size;
						const double my_radius_sq = my_radius*my_radius;
						double maxz=mp::material[mat].max*cs::system_dimensions[2];
						double minz=mp::material[mat].min*cs::system_dimensions[2];
						//double max_range = my_radius*my_radius*particle_radius_squared;
						// check for within core shell range
				      //if(range_x_sq*inv_r_sq + range_y_sq*inv_r_sq <= (2.0*cs::system_dimensions[2]-(cz))*(2.0*cs::system_dimensions[2]-(cz))*inv_c_sq*my_radius_sq){
				      if(range_x_sq*inv_r_sq + range_y_sq*inv_r_sq <= (L_cone-(cz))*(L_cone-(cz))*inv_c_sq*my_radius_sq){
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
				}*/
            // determine if atom is within bubble
				if(range_sq < factor_radius*factor_radius*particle_radius_sq){
					catom_array[atom].include=true;
					catom_array[atom].grain=grain;
				}
			} // end for loop

	      return;

      }  // end of bubble function

   }  // end internal namespace
}  // end create namespace
