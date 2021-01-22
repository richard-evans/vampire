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

//------------------------------------------------------------------------------
//   Function to determine atom material allocation based on z-height
//   For unit cells with multiple materials this allocation respects the
//   initial material allocation
//------------------------------------------------------------------------------
void layers(std::vector<cs::catom_t> & catom_array){

   // If z-height material selection is enabled then do so
	if(create::internal::select_material_by_z_height){

		// Check for interfacial roughness and call custom material assignment routine
		if(cs::interfacial_roughness) create::internal::roughness(catom_array);

      // Check for multilayer system and if required generate multilayers
      else if(cs::multilayers) cs::generate_multilayers(catom_array);

		// Otherwise perform normal assignment of materials
		else{

			// determine z-bounds for materials
			std::vector<double> mat_min(mp::num_materials);
			std::vector<double> mat_max(mp::num_materials);
         std::vector<bool> mat_fill(mp::num_materials);
         std::vector<int> uc_cat(mp::num_materials); // array of material -> unit cell material associations

         // Unroll min, max and fill for performance
			for(int mat=0;mat<mp::num_materials;mat++){
				mat_min[mat]=create::internal::mp[mat].min*cs::system_dimensions[2];
				mat_max[mat]=create::internal::mp[mat].max*cs::system_dimensions[2];
				// alloys generally are not defined by height, and so have max = 0.0
				if(mat_max[mat]<1.e-99) mat_max[mat]=-0.1;
            mat_fill[mat]=mp::material[mat].fill;
            uc_cat[mat] = create::internal::mp[mat].unit_cell_category; // unit cell category of material
			}

			// Assign materials to generated atoms
			for(unsigned int atom=0;atom<catom_array.size();atom++){
            const double cz=catom_array[atom].z;
            const int atom_uc_cat = catom_array[atom].uc_category;
				for(int mat=0;mat<mp::num_materials;mat++){
               // check atom within bounds, is not a fill material and has correct unit cell category
					if( (cz>=mat_min[mat]) && (cz<mat_max[mat]) && (mat_fill[mat]==false) && (atom_uc_cat == uc_cat[mat]) ){
						catom_array[atom].material=mat;
						catom_array[atom].include=true;
					}
				}
			}
		}

		// Delete unneeded atoms
      // Edit - (not sure why it does that here - causes unfiller atoms to be removed from system)
      // removed this call for now (June 2019)
		//clear_atoms(catom_array);

	}

   return;

}


} // end of internal namespace

} // end of create namespace
