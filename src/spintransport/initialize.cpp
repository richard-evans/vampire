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

// Vampire headers
#include "create.hpp"
#include "spintransport.hpp"
#include "vio.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{

   // simple struct to store 3D cell info
   struct cell3D_t{
      uint64_t id; // id of cell
      std::vector<uint64_t> atom; // list of atoms in each cell
      std::vector<uint64_t> nm_atom; // list of non-magnetic atoms in each cell
   };

   //----------------------------------------------------------------------------
   // Function to initialize spin transport module
   //----------------------------------------------------------------------------
   void initialize(const double system_size_x, // maximum dimensions of system along x-direction (angstroms)
                   const double system_size_y, // maximum dimensions of system along y-direction (angstroms)
                   const double system_size_z, // maximum dimensions of system along z-direction (angstroms)
                   const int num_materials,    // number of materials
                   const uint64_t num_atoms,   // number of atoms
                   const std::vector<int>& atoms_type_array, // material types of atoms
                   const std::vector<double>& atoms_x_coord_array, // x-coordinates of atoms
                   const std::vector<double>& atoms_y_coord_array, // y-coordinates of atoms
                   const std::vector<double>& atoms_z_coord_array, // z-coordinates of atoms
                   const std::vector<cs::nm_atom_t> non_magnetic_atoms_array // list of non-magnetic atoms
   ){

      //-------------------------------------------------------------------------
      // create 3D array of cells
      //
      //                 |       |       |       |       |
      //                 |       |       |       |       |
      //                 |       |       |       |       |
      //                 |       |       |       |       |
      //                 |       |       |       |       |
      //      stacks_z   |       |       |       |       |
      //                 |       |       |       |       |
      //    (current-    |       |   ^   |       |       |
      //     direction)  |       |   |   |       |       |
      //                 |       |   |   |       |       |
      //                 |       |   |   |       |       |
      //
      //                              stacks_x
      //
      //
      //
      //
      //
      //
      //
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // check that module is needed - if not do nothing
      //-------------------------------------------------------------------------
      if( st::internal::enabled == false ) return;

      // Inform user that module is initializing
      zlog << zTs() << "Initializing spin transport module" << std::endl;

      // If st material parameters uninitialised then initialise with default parameters
      if(spin_transport::internal::mp.size() == 0){
         spin_transport::internal::mp.resize(num_materials);
      }

      // copy cell sizes to local variable to enable arbitrary current direction
      const double cell_size[3] = { st::internal::cell_size_x,
                                    st::internal::cell_size_y,
                                    st::internal::cell_size_z };

      // calculate number of cells in each direction x,y,z (rounding up)
      const int num_cells[3] = { ceil(system_size_x/cell_size[0]),
                                 ceil(system_size_y/cell_size[1]),
                                 ceil(system_size_z/cell_size[2]) };


      int stack_x = 0; // spatial direction of stack arrays  x,y,z
      int stack_y = 1; // (each mapping into a physical direction)
      int stack_z = 2;

      // calculate stack order based on current direction (always calculated along the stack-z direction)
      if(st::internal::current_direction == st::internal::px){
         stack_x = 1; // y-direction
         stack_y = 2; // z-direction
         stack_z = 0; // x-direction
      }

      //---------------------------------------------------------------------------------
      // Calculate cell IDs in 3D to calculate atom-cell associations
      //---------------------------------------------------------------------------------
      std::vector< std::vector < std::vector <cell3D_t> > > cells3D; // 3D list of cell IDs

      int cell_id =0;

      // resize cell array to hold stack x_cells
      cells3D.resize(num_cells[stack_x]);

      // loop over all x-cells in stack
      for(unsigned int i = 0; i < cells3D.size(); i++){

         // resize cell array[i] to hold stack y_cells
         cells3D[i].resize(num_cells[stack_y]);

         // loop over all y-cells in stack
         for(unsigned int j = 0; j < cells3D[i].size(); j++){

            // resize cell array[i][j] to hold stack z_cells (each stack linear in memory along z (current) direction)
            cells3D[i][j].resize(num_cells[stack_z]); //

            // loop over all z-cells to set linear cell ID for cell arrays
            for(unsigned int k = 0; k < cells3D[i][j].size(); k++){
               cells3D[i][j][k].id = cell_id;
               cell_id++; // increment cell number
            }

         }

      }

      //---------------------------------------------------------------------------------
      // Calculate atom-cell associations and cell positions
      //---------------------------------------------------------------------------------

      // resize array to store which cell each atom is in
      st::internal::atom_in_cell.resize(num_atoms);

      for(uint64_t atom = 0; atom < num_atoms; atom++){

         // get atomic coordinates
         const double cx = atoms_x_coord_array[atom];
         const double cy = atoms_y_coord_array[atom];
         const double cz = atoms_z_coord_array[atom];

         // calculate 3D cell ID in atoms coordinate system
         unsigned int cellID[3] = { cx/cell_size[0], cy/cell_size[1], cz/cell_size[2]};

         // now calculate in stack coordinate system
         const uint64_t i = cellID[stack_x];
         const uint64_t j = cellID[stack_y];
         const uint64_t k = cellID[stack_z];

         // associate atom with stack cell ID
         st::internal::atom_in_cell[atom] = cells3D[i][j][k].id;

         // add atom to atoms in cell list
         cells3D[i][j][k].atom.push_back(atom);

      }

      // include non-magnetic atoms
      const uint64_t num_nm_atoms = non_magnetic_atoms_array.size();
      for(uint64_t atom = 0; atom < num_nm_atoms; atom++){

         // get atomic coordinates
         const double cx = non_magnetic_atoms_array[atom].x;
         const double cy = non_magnetic_atoms_array[atom].y;
         const double cz = non_magnetic_atoms_array[atom].z;

         // calculate 3D cell ID in atoms coordinate system
         unsigned int cellID[3] = { cx/cell_size[0], cy/cell_size[1], cz/cell_size[2]};

         // now calculate in stack coordinate system
         const uint64_t i = cellID[stack_x];
         const uint64_t j = cellID[stack_y];
         const uint64_t k = cellID[stack_z];

         // add atom to atoms in non-magnetic cell list
         cells3D[i][j][k].nm_atom.push_back(atom);

      }

      //---------------------------------------------------------------------------------
      // Calculate stack data
      //---------------------------------------------------------------------------------

      // calculate number of stacks (x*y)
      st::internal::num_stacks = num_cells[stack_x]*num_cells[stack_y];
      const uint64_t num_cells_in_stack = num_cells[stack_z]; // save number of cells in stack along current direction

      // resize limits for each stack
      st::internal::stack_start_index.resize(st::internal::num_stacks);
      st::internal::stack_final_index.resize(st::internal::num_stacks);

      // determine initial start and end cell of each stack
      for(int s = 0; s < st::internal::num_stacks; s++){
         st::internal::stack_start_index[s] = s*num_cells_in_stack;
         st::internal::stack_final_index[s] = s*num_cells_in_stack + num_cells_in_stack; // loop to less than this number
      }

      // calculate total number of cells
      st::internal::total_num_cells = num_cells[0]*num_cells[1]*num_cells[2];

      // resize array to store sequential list of magnetic cells to calculate the resistance over
      st::internal::next_cell_in_stack.resize(st::internal::total_num_cells); // list of next cell in stack to account for tunnel barrier

      // resize arrays to store average resistance in each cell
      st::internal::cell_resistance.resize(st::internal::total_num_cells,0.0);
      st::internal::cell_spin_resistance.resize(st::internal::total_num_cells,0.0);

      // cell size parameters for resitivity to resistance calculation
      const double iA = 1.0 / (cell_size[stack_x] * cell_size[stack_y]); // Angstroms^-2
      const double l = cell_size[stack_z];

      // loop over all xy-cells (stacks)
      for(unsigned int i = 0; i < cells3D.size(); i++){
         for(unsigned int j = 0; j < cells3D[i].size(); j++){

            // loop over all cells in stack
            for(unsigned int k = 0; k < cells3D[i][j].size(); k++){

               // determine cell ID
               const uint64_t cell = cells3D[i][j][k].id;

               // determine total number of atoms
               uint64_t num_atoms_in_cell = cells3D[i][j][k].atom.size() + cells3D[i][j][k].nm_atom.size();

               // if cell is empty assume uniform padding resistance
               if(num_atoms_in_cell == 0){
                  st::internal::cell_resistance[cell]      = st::internal::environment_resistivity * l * iA;
                  st::internal::cell_spin_resistance[cell] = st::internal::environment_resistivity * l * iA;
               }
               // otherwise add contributions from atoms to calculate average resistivity
               else{

                  // variables to accumulate resistances for cell
                  double resistivity = 0.0;
                  double spin_resistivity = 0.0;

                  // magnetic atoms
                  for(unsigned int atom = 0; atom < cells3D[i][j][k].atom.size(); atom++ ){
                     int mat = atoms_type_array[atom]; // get material type
                     resistivity += st::internal::mp[mat].resistivity; // add resistivity to total
                     spin_resistivity += st::internal::mp[mat].spin_resistivity; // add spin resistivity to total
                  }

                  // non-magnetic atoms
                  for(unsigned int atom = 0; atom < cells3D[i][j][k].nm_atom.size(); atom++ ){
                     int mat = non_magnetic_atoms_array[atom].mat; // get material type
                     resistivity += st::internal::mp[mat].resistivity; // add resistivity to total
                     spin_resistivity += st::internal::mp[mat].spin_resistivity; // add spin resistivity to total
                  }

                  // calculate average resistivity
                  double mean_resistivity      = resistivity      / double( num_atoms_in_cell );
                  double mean_spin_resistivity = spin_resistivity / double( num_atoms_in_cell );

                  // set cell resistance
                  st::internal::cell_resistance[cell]      = mean_resistivity * l * iA;
                  st::internal::cell_spin_resistance[cell] = mean_spin_resistivity * l * iA;

               }


            }
         }
      }

      st::internal::cell_magnetization.resize(st::internal::total_num_cells);
      st::internal::cell_position.resize(3*st::internal::total_num_cells);

      return;

   }

} // end of spin_transport namespace
