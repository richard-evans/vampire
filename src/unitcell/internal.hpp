//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016, Jack B Collings 2021. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef UNITCELL_INTERNAL_H_
#define UNITCELL_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the unitcell module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <sstream>

// Vampire headers
#include "unitcell.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------
      enum exchange_function_t { nearest_neighbour, shell, exponential, material_exponential, RKKY};

      //-------------------------------------------------------------------------
      // Internal shared variables set by input file
      //-------------------------------------------------------------------------
	   extern std::string crystal_structure;
	   extern std::string unit_cell_filename;

      extern double unit_cell_size_x;
      extern double unit_cell_size_y;
      extern double unit_cell_size_z;


      extern exchange_function_t exchange_function; // exchange function type
      extern double exchange_interaction_range; // multiples of interaction radius used
      extern double exchange_decay; // exponential function parameter
      extern double exchange_multiplier; // exponential function parameter
      extern double exchange_shift; // exponential function parameter
      extern double RKKYkf; // RKKY function parameter
      extern std::vector <std::vector <exchange_parameters_t> > material_exchange_parameters; // holds exponential exchange parameters for material-wise exchange
      extern std::vector <std::vector <double> > nn_cutoff_range; // holds nearest neighbour cut-off ranges between different unit cell categories.
      extern std::vector <std::vector <double> > interaction_cutoff_range; // holds interaction cut-off ranges.

      extern bool sublattice_materials; // flag to enable identification of atoms in simple crystals by material

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      void build_simple_cubic(unitcell::unit_cell_t& unit_cell);
      void build_body_centred_cubic(unitcell::unit_cell_t& unit_cell);
      void build_body_centred_cubic_110(unitcell::unit_cell_t& unit_cell);
      void build_face_centred_cubic(unitcell::unit_cell_t& unit_cell);
      void build_face_centred_cubic_111(unitcell::unit_cell_t& unit_cell);
      void build_hexagonal_close_packed(unitcell::unit_cell_t& unit_cell);
      void build_honeycomb(unitcell::unit_cell_t& unit_cell);
      void build_honeycomb_alpha(unitcell::unit_cell_t& unit_cell);
      void build_honeycomb_beta(unitcell::unit_cell_t& unit_cell);
      void build_kagome(unitcell::unit_cell_t& unit_cell);
      void build_heusler(unitcell::unit_cell_t& unit_cell);
      void build_mn2au(unitcell::unit_cell_t& unit_cell);
      void build_NdFeB(unitcell::unit_cell_t& unit_cell);
      void build_rock_salt(unitcell::unit_cell_t& unit_cell);
      void build_SmFeN(unitcell::unit_cell_t& unit_cell);
      void build_spinel(unitcell::unit_cell_t& unit_cell);
      void build_spinel_layered(unitcell::unit_cell_t& unit_cell);

      void calculate_interactions(unit_cell_t& unit_cell);
      void read_unit_cell(unit_cell_t & unit_cell, std::string filename);
      void read_biquadratic_interactions(unit_cell_t & unit_cell,
                                         std::stringstream& ucf,
                                         std::istringstream& ucf_ss,
                                         std::string& filename,
                                         unsigned int& line_counter,
                                         int& interaction_range);
      void verify_exchange_interactions(unit_cell_t & unit_cell, std::string filename);
      double exchange(double range, double cutoff, int mat_i, int mat_j);

      void write_unit_cell_file(unit_cell_t & uc);

   } // end of internal namespace

} // end of unitcell namespace

#endif //UNITCELL_INTERNAL_H_
