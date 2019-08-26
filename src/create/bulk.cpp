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
   // Function to include all atoms in simulation
   //------------------------------------------------------------
   void bulk(std::vector<cs::catom_t> & catom_array){

   	// Loop over all atoms and mark as selected
   	const int num_atoms = catom_array.size();

    	for(int atom=0; atom < num_atoms; atom++){
   		catom_array[atom].include=true;
   	}

   	return;

   }

} // end of internal namespace

} // end of create namespace
