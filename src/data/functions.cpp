//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo 2022. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"

// // exchange module headers
// #include "internal.hpp"

namespace atoms{

   //--------------------------------------------------------------------------------
   // Function to initialise previous spin configuration arrays 
   //--------------------------------------------------------------------------------
   void initialise_old_spins(){
      const uint16_t size = atoms::x_spin_array.size();
		atoms::x_old_spin_array.resize(size, 0.0);
		atoms::y_old_spin_array.resize(size, 0.0);
		atoms::z_old_spin_array.resize(size, 0.0);
      return;
   }

   //--------------------------------------------------------------------------------
   // Function to save previous spin configuration  
   //--------------------------------------------------------------------------------
   void store_old_spins(){
		for(int atom=0;atom<num_atoms;atom++){
			atoms::x_old_spin_array[atom] = atoms::x_spin_array[atom];
			atoms::y_old_spin_array[atom] = atoms::y_spin_array[atom];
			atoms::z_old_spin_array[atom] = atoms::z_spin_array[atom];
		}
      return;
   }

} // end of atom namespace