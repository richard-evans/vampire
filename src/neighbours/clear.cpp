//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "neighbours.hpp"

namespace neighbours{

   //----------------------------------------------------------------------------------
   // Function to empty data from neighbourlist
   //----------------------------------------------------------------------------------
void list_t::clear(){

   // force deallocation by making main object data go out of scope
   // Everybody who loves C++ scoping rules say woo!

   // a simple unallocated array of neighbours
   std::vector< std::vector <neighbours::neighbour_t> > tmp;

   // swap the pointers
   tmp.swap(list);

   // leaving unloved memory behind
   return;

}

} // end of namespace neighbours
