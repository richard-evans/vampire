//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "create.hpp"
#include "material.hpp"

// Internal create header
#include "internal.hpp"

namespace create{
namespace internal{

void alloy(std::vector<cs::catom_t> & catom_array){

	// loop over all atoms
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		// if atom material is alloy master then reassign according to % chance
		int local_material=catom_array[atom].material;
		if(create::internal::mp[local_material].alloy_master==true){
			//loop over all potential alloy materials
			for(int mat=0;mat<mp::num_materials;mat++){
				double probability = create::internal::mp[local_material].slave_material[mat].fraction;
				if(create::internal::grnd() < probability){
					catom_array[atom].material=mat;
				}
			}
      }
   }

	return;

}

} // end of internal namespace
} // end of create namespace
