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

#ifndef NEIGHBOURS_H_
#define NEIGHBOURS_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "create_atoms_class.hpp"
#include "exchange_types.hpp" // needed for exchange interaction type definition
#include "neighbours.hpp"
#include "unitcell.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for neighbours module
//--------------------------------------------------------------------------------
namespace neighbours{

   //-----------------------------------------------------------------------------
   // Simple class of neighbour type to store information about the neighbour
   //-----------------------------------------------------------------------------
   class neighbour_t {
	public:

		int nn; // atom id of neighbour
		int i; // interaction type of neighbour

      double vx; // real coordinate vector between atoms i->j
      double vy;
      double vz;

	};

   //-----------------------------------------------------------------------------
   // Simple class of neighbour list definining a set of interactions
   //-----------------------------------------------------------------------------
   class list_t{
   public:

      // list of neighbours in terms of atom IDs
      std::vector<std::vector <neighbours::neighbour_t> > list;

      // generate neighbour list from interaction template and list of atoms
      void generate(std::vector<cs::catom_t>& atoms,
                    unitcell::exchange_template_t& exchange,
                    const unsigned int num_atoms_in_unit_cell,
                    double ucdx, double ucdy, double ucdz);

      // release neighbour list
      void clear();

   };

} // end of neighbours namespace

#endif //NEIGHBOURS_H_
