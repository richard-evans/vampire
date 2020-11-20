//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
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

      void cone(std::vector<double>& particle_origin, std::vector<cs::catom_t> & catom_array, const int grain){
         //----------------------------------------------------------------------
         // Function to cut a truncated cone
         //
         // (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
         //
         // (x-h)^2  +  (y-k)^2    (2*L-z-l)^2
         // --------  --------- = -----------  (z-l-aL)^2
         //  (r)^2      (r)^2       (2*L)^2
         //
         // Evaluation of the above for particle coordinates h,k,l
         // and atomic positions x,y,z will include all atoms within
         // cone with vertex L, radius r
         //
         //----------------------------------------------------------------------
         double PI=3.14159265358979323846264338327;

			//-----------------------------------------
			// Set particle radius
			//-----------------------------------------
			double particle_radius_squared = (cs::particle_scale*0.5)*(cs::particle_scale*0.5);
         double alpha  = 90.0 - create::internal::cone_angle;
         double L_cone = (cs::particle_scale*0.5)*tan(alpha*PI/180.);
         // Use shape modifiers to generate ellipsoid
         const double inv_r_sq = 1.0/(particle_radius_squared);
         //const double inv_c_sq = 1.0/(cs::system_dimensions[2]*cs::system_dimensions[2]*2.0*2.0);
         const double inv_c_sq = 1.0/(L_cone*L_cone);

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
            const double range_x_sq = (catom_array[atom].x-particle_origin[0])*(catom_array[atom].x-particle_origin[0]);
            const double range_y_sq = (catom_array[atom].y-particle_origin[1])*(catom_array[atom].y-particle_origin[1]);
            //const double range_z_sq = (catom_array[atom].z-particle_origin[2])*(catom_array[atom].z-particle_origin[2])*4.0;
				double cz=catom_array[atom].z;

				if(mp::material[catom_array[atom].material].core_shell_size>0.0){
               // Iterate over materials
               for(std::list<core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
                  int mat = (it)->mat;
						double my_radius = mp::material[mat].core_shell_size;
						const double my_radius_sq = my_radius*my_radius;
						double maxz=create::internal::mp[mat].max*cs::system_dimensions[2];
						double minz=create::internal::mp[mat].min*cs::system_dimensions[2];
                  const int atom_uc_cat = catom_array[atom].uc_category;
                  const int mat_uc_cat = create::internal::mp[mat].unit_cell_category;
						//double max_range = my_radius*my_radius*particle_radius_squared;
						// check for within core shell range
				      //if(range_x_sq*inv_r_sq + range_y_sq*inv_r_sq <= (2.0*cs::system_dimensions[2]-(cz))*(2.0*cs::system_dimensions[2]-(cz))*inv_c_sq*my_radius_sq){
				      if(range_x_sq*inv_r_sq + range_y_sq*inv_r_sq <= (L_cone-(cz))*(L_cone-(cz))*inv_c_sq*my_radius_sq){
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
				//else if(range_x_sq*inv_r_sq + range_y_sq*inv_r_sq <= (2.0*cs::system_dimensions[2]-(cz))*(2.0*cs::system_dimensions[2]-(cz))*inv_c_sq){
				else if(range_x_sq*inv_r_sq + range_y_sq*inv_r_sq <= (L_cone-(cz))*(L_cone-(cz))*inv_c_sq){
					catom_array[atom].include=true;
					catom_array[atom].grain=grain;
				}
			} // end for loop
	      return;
      }  // end cone function
   }  // end internal namespace
}  // end create namespace
