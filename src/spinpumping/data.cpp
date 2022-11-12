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
#include "spinpumping.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_pumping{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside spintransport module
      //------------------------------------------------------------------------
      bool enabled = false; // bool to enable spin transport calculation
      bool output_atomistic_spin_pumping_flag = false; // flag to toggle output of atomic resolution spin current
      bool output_cells_spin_pumping_flag = false; // flag to toggle output of cell-resolved spin current

      unsigned int update_rate  = 100;  // number of timesteps between updates
      unsigned int time_counter = 100;  // number of timesteps since last update (initially set to same as update rate to ensure calculation at start)
      uint64_t config_counter = 0;  // counter for update of spin pumping configs to file

      std::vector<internal::mp_t> mp; // array of material properties

      double cell_size_x = 10.0; // cell size along x-direction (1 nm default size)
      double cell_size_y = 10.0; // cell size along y-direction (1 nm default size)
      double cell_size_z = 10.0; // cell size along z-direction (1 nm default size)

      unsigned int total_num_cells  = 0; // number of cells

      // arrays to store cell properties
      std::vector <bool> cell_magnetic;               // boolean array to determine if cell is magnetic or not
      std::vector <double> cell_alpha;                // cell magnetization (average of constituent atoms)
      std::vector <double> cell_magnetization;        // 3N normalised magnetization in each cell
      std::vector <double> cell_isaturation;          // inverse magnetic saturation at T=0 in each cell
      std::vector <double> cell_position;             // 3N array of cell positions (origin)
      std::vector <double> cell_spin_mix_conductance; // array to store spin mixing conductance in each cell

      // array to store which cell each atom is in
      std::vector <unsigned int> atom_in_cell;

      std::vector <bool> material_magnetic;                      // boolean array to determine if material is magnetic or not
      std::vector <int> atoms_type_array;                    // array to store material type of atoms

      std::vector <double> x_atom_spin_pumping_array;        // arrays to store atom spin pumping x-component
      std::vector <double> y_atom_spin_pumping_array;        // arrays to store atom spin pumping y-component
      std::vector <double> z_atom_spin_pumping_array;        // arrays to store atom spin pumping z-component

      std::vector <double> x_cell_spin_pumping_array;        // arrays to store cell spin pumping x-component
      std::vector <double> y_cell_spin_pumping_array;        // arrays to store cell spin pumping y-component
      std::vector <double> z_cell_spin_pumping_array;        // arrays to store cell spin pumping z-component

   } // end of internal namespace

} // end of spin_transport namespace
