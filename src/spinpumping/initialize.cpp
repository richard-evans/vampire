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
#include "create.hpp"
#include "spintransport.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_pumping{

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
                   const uint64_t num_atoms,   // number of local atoms
                   const std::vector<int>& atoms_type_array, // material types of atoms
                   const std::vector<double>& atoms_x_coord_array, // x-coordinates of atoms
                   const std::vector<double>& atoms_y_coord_array, // y-coordinates of atoms
                   const std::vector<double>& atoms_z_coord_array, // z-coordinates of atoms
                   const std::vector<double>& atoms_m_spin_array,  // moments of atoms (muB)
                   const std::vector<double>& material_damping_array, // array of material level damping constants
                   const std::vector<bool>& is_magnetic_material, // array of size num_mat to state whether material is magnetic (true) or not (false)
                   const std::vector<cs::nm_atom_t> non_magnetic_atoms_array // list of non-magnetic atoms
   ){

      //-------------------------------------------------------------------------
      // check that module is needed - if not do nothing
      //-------------------------------------------------------------------------
      if( spin_pumping::internal::enabled == false ) return;

      // Inform user that module is initializing
      zlog << zTs() << "Initializing spin pumping module" << std::endl;

      // If st material parameters uninitialised then initialise with default parameters
      if(spin_pumping::internal::mp.size() != static_cast<size_t>(num_materials)){
         spin_pumping::internal::mp.resize(num_materials);
      }

      // Inform user that module is initializing structures for atomistic calculation
      zlog << zTs() << "Initializing structures for atomistic calculation of spin pumping" << std::endl;
      spin_pumping::internal::material_magnetic = is_magnetic_material;
      spin_pumping::internal::atoms_type_array = atoms_type_array;
      // Resize arrays to store cross product spin with time derivative of spin
      spin_pumping::internal::x_atom_spin_pumping_array.resize(num_atoms, 0.0);
      spin_pumping::internal::y_atom_spin_pumping_array.resize(num_atoms, 0.0);
      spin_pumping::internal::z_atom_spin_pumping_array.resize(num_atoms, 0.0);

      // Inform user that module is initializing structures for cells calculation
      zlog << zTs() << "Initializing structures for cell-resolved calculation of spin pumping" << std::endl;
      // Resize arrays to store cross product spin with time derivative of spin
      spin_pumping::internal::x_cell_spin_pumping_array.resize(spin_pumping::internal::total_num_cells, 0.0);
      spin_pumping::internal::y_cell_spin_pumping_array.resize(spin_pumping::internal::total_num_cells, 0.0);
      spin_pumping::internal::z_cell_spin_pumping_array.resize(spin_pumping::internal::total_num_cells, 0.0);

      // If enabled, output calculated atomistic coordinates and moments (passing local values)
      if(spin_pumping::internal::output_atomistic_spin_pumping_flag) {
         // zlog << zTs() << "Outputting atomistic coordinates for atomistic calculation of spin pumping" << std::endl;
         spin_pumping::internal::output_atomistic_coordinates(num_atoms, atoms_x_coord_array, atoms_y_coord_array, atoms_z_coord_array, atoms_m_spin_array);
      }
      // If enabled, output calculated atomistic coordinates and moments (passing local values)
      if(spin_pumping::internal::output_cells_spin_pumping_flag) {
         // zlog << zTs() << "Outputting atomistic coordinates for atomistic calculation of spin pumping" << std::endl;
         // spin_pumping::internal::output_atomistic_coordinates(num_atoms, atoms_x_coord_array, atoms_y_coord_array, atoms_z_coord_array, atoms_m_spin_array);
      }

      return;

   }

} // end of spin_pumping namespace
