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

   void ellipsoid(std::vector<double>& particle_origin,std::vector<cs::catom_t> & catom_array, const int grain){
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

      // Set particle radius
      const double particle_radius_squared = (cs::particle_scale * 0.5) * (cs::particle_scale * 0.5);

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
         const double range_x_sq = (catom_array[atom].x-particle_origin[0])*(catom_array[atom].x-particle_origin[0]);
         const double range_y_sq = (catom_array[atom].y-particle_origin[1])*(catom_array[atom].y-particle_origin[1]);
         const double range_z_sq = (catom_array[atom].z-particle_origin[2])*(catom_array[atom].z-particle_origin[2]);

         if(mp::material[catom_array[atom].material].core_shell_size>0.0){
            // Iterate over materials
            for(std::list<core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
               int mat = (it)->mat;
               const double my_radius_sq = mp::material[mat].core_shell_size*mp::material[mat].core_shell_size;
               double maxz = mat_max[mat];
               double minz = mat_min[mat];
               double cz = catom_array[atom].z;
               const int atom_uc_cat = catom_array[atom].uc_category;
               // check for within core shell range
               if(range_x_sq * inv_rx_sq + range_y_sq * inv_ry_sq + range_z_sq * inv_rz_sq <= my_radius_sq){
                  if( (cz >= minz) && (cz < maxz) && (atom_uc_cat == uc_cat[mat]) ){
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
         else if(range_x_sq * inv_rx_sq + range_y_sq * inv_ry_sq + range_z_sq * inv_rz_sq <= 1.0) {
            catom_array[atom].include = true;
            catom_array[atom].grain = grain;
         }
      }

      return;

   }


} // end of internal namespace

} // end of create namespace
