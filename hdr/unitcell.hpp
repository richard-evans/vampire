//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016, Jack B. Collings 2021. All rights reserved.
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
#include "exchange_types.hpp" // needed for exchange interaction type definition
#include "unitcell.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for unitcell module
//--------------------------------------------------------------------------------
namespace unitcell{

   //---------------------------------------------------------------------------
   // Unit cell atom class definition
   //---------------------------------------------------------------------------
   class atom_t {
	public:
      double x; /// atom x-coordinate
      double y; /// atom y-coordinate
      double z; /// atom z-coordinate
      unsigned int mat; /// material
      unsigned int lc; /// lattice category
      unsigned int hc; /// height category
      unsigned int ni; /// number of interactions
      bool nm; // non-magnetic atom (no interactions are calculated)

      // constructor
      atom_t():
         x(0.0),
         y(0.0),
         z(0.0),
         mat(0),
         lc(0),
         hc(0),
         ni(0),
         nm(false)
      {
      };

	};

   //---------------------------------------------------------------------------
   // Unit cell interaction class definition
   //---------------------------------------------------------------------------
	class interaction_t {
	public:
      unsigned int i; /// atom unit cell id
      unsigned int j; /// neighbour atom unit cell id
      unsigned int mat_i; /// atom material category
      unsigned int mat_j; /// neighbour material category
      unsigned int shell; // shell number of interaction
      int dx; /// delta x in unit cells
      int dy; /// delta y in unit cells
      int dz; /// delta z in unit cells
      double rij; // interaction range (unit cells)
      double Jij[3][3]; /// Exchange tensor
	};

   //------------------------------------------------------------------------
   // Class to hold exchange function parameters
   //------------------------------------------------------------------------
   class exchange_parameters_t{
      public:
         double decay_length;
         double decay_multiplier;
         double decay_shift;

      exchange_parameters_t():
         decay_length(0.4),
         decay_multiplier(0.0),
         decay_shift(0.0)
      {
      };
   };

   //---------------------------------------------------------------------------
   // Unit cell exchange template class definition
   //---------------------------------------------------------------------------
   class exchange_template_t {
   public:

      exchange::exchange_t exchange_type; // exchange type to use in simulation
      bool use_material_exchange_constants; // flag to enable material exchange parameters
      int num_unit_cell_atoms; // number of atoms in unit cell

      // list of interactions in each unit cell
      std::vector <unitcell::interaction_t> interaction;

      // list of number of interactions from template for each atom in unit cell
      std::vector <int> ni;

      // Class constructor with rational defaults
      exchange_template_t():
         exchange_type(exchange::isotropic),
         use_material_exchange_constants(true)
      {
         // Do nothing
         return;
      };

      void read_interactions(
         const int num_atoms, // num atoms in unit cell
         std::stringstream& ucf_file,
         std::istringstream& ucf_ss,
         std::string& filename,
         unsigned int& line_counter,
         unsigned int& interaction_range);

      // function to set exchange type
      unsigned int set_exchange_type(std::string exchange_type_string);

      // function to verify exchange interactions are reciprocal
      void verify(std::string filename);

      // normalisation function to achieve same exchange sum as nn approximation
      void normalise_exchange(std::vector < std::vector <double> > &nn_cutoff_range);
      void normalise_functional_exchange(std::vector < std::vector <double> > &nn_cutoff_range);

      // function to find crystal shells
      void find_shells();

   };

   //---------------------------------------------------------------------------
   // Unit cell class definition
   //---------------------------------------------------------------------------
	class unit_cell_t {
	public:

		double dimensions[3];
		double shape[3][3];
      double cutoff_radius; // nearest neighbours
      unsigned int interaction_range; /// maximum range in unit cells

		unsigned int lcsize; /// number of local categories
		unsigned int hcsize; /// number of height categories
		unsigned int surface_threshold; /// threshold for surface atoms

		// list of atoms in each unit cell
		std::vector <unitcell::atom_t> atom;

      unitcell::exchange_template_t bilinear;
      unitcell::exchange_template_t biquadratic;
      //exchange_template_t fourspin_interaction; // tbc

	};

   //-----------------------------------------------------------------------------
   // Function to initialise unitcell module
   //-----------------------------------------------------------------------------
   void initialise(unit_cell_t& unit_cell);

   //---------------------------------------------------------------------------
   // Function to process input file parameters for unitcell module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line, int const superIndex, int const subIndex);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

   void set_crystal_structure_to_simple_cubic();

} // end of unitcell namespace

// alias namespace
namespace uc = unitcell;

#endif //UNITCELL_H_
