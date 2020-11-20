//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <list>

// Vampire headers
#include "create.hpp"

// Internal create header
#include "internal.hpp"

// comparison function
bool compare(cs::catom_t first,cs::catom_t second){
	if(first.grain<second.grain) return true;
	else return false;
}

namespace create{
namespace internal{

//------------------------------------------------------------------------------
// Function to sort atoms by grain number (for improved performance)
//------------------------------------------------------------------------------
void sort_atoms_by_grain(std::vector<cs::catom_t> & catom_array){

   // Get number of atoms
   const int num_atoms=catom_array.size();

   // Create list object
   std::list <cs::catom_t> catom_list(num_atoms);

   // copy data to list
   copy(catom_array.begin(), catom_array.end(), catom_list.begin());

   // sort date in list
   catom_list.sort(compare);

   // copy list to data
   copy(catom_list.begin(), catom_list.end(), catom_array.begin());

   return;

}

} // end of namespace internal
} // end of namespace create
