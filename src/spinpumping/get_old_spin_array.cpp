//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard Evans 2022. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "spinpumping.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace spin_pumping{
	namespace internal{
   	//--------------------------------------------------------------------------------
   	// Function to recover old spin configuration x component
   	//--------------------------------------------------------------------------------
   	std::vector<double> get_old_spins_x(const unsigned num_atoms){
      	std::vector<double> out;
			for(int atom=0;atom<num_atoms;atom++){ out.push_back(atoms::x_old_spin_array[atom]);}
      	return out;
   	}

   	//--------------------------------------------------------------------------------
   	// Function to recover old spin configuration y component
   	//--------------------------------------------------------------------------------
   	std::vector<double> get_old_spins_y(const unsigned num_atoms){
      	std::vector<double> out;
			for(int atom=0;atom<num_atoms;atom++){ out.push_back(atoms::y_old_spin_array[atom]);}
      	return out;
   	}

   	//--------------------------------------------------------------------------------
   	// Function to recover old spin configuration z component
   	//--------------------------------------------------------------------------------
   	std::vector<double> get_old_spins_z(const unsigned num_atoms){
      	std::vector<double> out;
			for(int atom=0;atom<num_atoms;atom++){ out.push_back(atoms::z_old_spin_array[atom]);}
      	return out;
   	}

   } // end of namespace internal
} // end of namespace
