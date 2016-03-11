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
		if(mp::material[local_material].alloy_master==true){
		  // now check for unordered alloy
			if(mp::material[local_material].alloy_class==-1){
				//loop over all potential alloy materials
				for(int mat=0;mat<mp::num_materials;mat++){
					double probability = mp::material[local_material].alloy[mat];
					if(create::internal::grnd() < probability){
						catom_array[atom].material=mat;
					}
				}
			}
			// if not ordered, then assume ordered
			else{
				// loop over all alloy materials
				for(int mat=0;mat<mp::num_materials;mat++){
					// get class of alloy material
					int alloy_class = mp::material[mat].alloy_class;
					// check for matching class and uc
					// ----- DOES NOT WORK for > 1 alloy!!
					// need to check for correct alloy master material ------------
					if(catom_array[atom].uc_category==alloy_class){
						// set material
						catom_array[atom].material=mat;
					}
				}
			}
		}
	}

	return;
}

} // end of internal namespace
} // end of create namespace
