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
#include "vio.hpp"

// Internal create header
#include "internal.hpp"

namespace create{
namespace internal{

//--------------------------------------------------------------------
// Function to remove atoms that are outside of the desired gemoetry
//--------------------------------------------------------------------
void clear_atoms(std::vector<cs::catom_t> & catom_array){

   // Get original and new number of atoms
   const int num_atoms=catom_array.size();
   int num_included=0;
   for(int a=0;a<num_atoms;a++){
      if(catom_array[a].include == true && mp::material[catom_array[a].material].non_magnetic != 1){
         num_included++;
      }
   }

   // check if there are unneeded atoms
   if(num_atoms!=num_included){
      // create temporary copy for atoms
      std::vector<cs::catom_t> tmp_catom_array(num_atoms);
      tmp_catom_array=catom_array;
      // resize original array to new number of atoms
      catom_array.resize(num_included);
      int atom=0;
      // loop over all existing atoms
      for(int a=0;a<num_atoms;a++){
         // if atom is to be included and is non-magnetic copy to new array
         if(catom_array[a].include==true && mp::material[catom_array[a].material].non_magnetic != 1 ){
            catom_array[atom]=tmp_catom_array[a];
            atom++;
         }
         // if atom is part of a non-magnetic material then save to nm array
         else if(catom_array[a].include == true && mp::material[catom_array[a].material].non_magnetic == 1){
            cs::nm_atom_t tmp;
         	tmp.x = catom_array[a].x;
         	tmp.y = catom_array[a].y;
         	tmp.z = catom_array[a].z;
         	tmp.mat = catom_array[a].material;
         	tmp.element = mp::material[catom_array[a].material].element;
         	// save atom to non-magnet array
         	cs::non_magnetic_atoms_array.push_back(tmp);
         }
      }
      tmp_catom_array.resize(0);

      zlog << zTs() << "Removed " << cs::non_magnetic_atoms_array.size() << " non-magnetic atoms from system" << std::endl;

   }

   return;

}

}
}
