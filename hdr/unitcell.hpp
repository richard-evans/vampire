//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef UNITCELL_H_
#define UNITCELL_H_

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers
#include "unitcell.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for unitcell module
//--------------------------------------------------------------------------------
namespace unitcell{

   // Unit cell atom class definition
   class atom_t {
	public:
		double x; /// atom x-coordinate
		double y; /// atom y-coordinate
		double z; /// atom z-coordinate
		unsigned int mat; /// material
		unsigned int lc; /// lattice category
		unsigned int hc; /// height category
		unsigned int ni; /// number of interactions
      unsigned int nbqi; /// number of interactions
	};

   // Unit cell interaction class definition
	class interaction_t {
	public:
		unsigned int i; /// atom unit cell id
		unsigned int j; /// neighbour atom unit cell id
		int dx; /// delta x in unit cells
		int dy; /// delta y in unit cells
		int dz; /// delta z in unit cells
      double rij; // interaction range (unit cells)
		double Jij[3][3]; /// Exchange tensor
	};

   // Unit cell class definition
	class unit_cell_t {
	public:

		double dimensions[3];
		double shape[3][3];
      double cutoff_radius; // nearest neighbours

		unsigned int lcsize; /// number of local categories
		unsigned int hcsize; /// number of height categories
		unsigned int interaction_range; /// maximum range in unit cells
		unsigned int surface_threshold; /// threshold for surface atoms

		// list of atoms in each unit cell
		std::vector <unitcell::atom_t> atom;

		// list of interactions in each unit cell
		std::vector <unitcell::interaction_t> interaction;

      // list of biquadratic interactions in each unit cell
		std::vector <unitcell::interaction_t> biquadratic_interaction;

	};

   //-----------------------------------------------------------------------------
   // Function to initialise unitcell module
   //-----------------------------------------------------------------------------
   void initialise(unit_cell_t& unit_cell);

   //---------------------------------------------------------------------------
   // Function to process input file parameters for unitcell module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

   void set_crystal_structure_to_simple_cubic();

} // end of unitcell namespace

// alias namespace
namespace uc = unitcell;

#endif //UNITCELL_H_
