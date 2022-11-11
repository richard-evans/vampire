//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard Evans 2022. All rights reserved.
//
//------------------------------------------------------------------------------
//

#ifndef SPINPUMPING_INTERNAL_H_
#define SPINPUMPING_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the spintransport module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "spinpumping.hpp"

// spinpumping module headers

namespace spin_pumping{

   namespace internal{
      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------
      // simple initialised class for set variables
      class set_double_t{

      private:
         double value; // value
         bool setf; // flag specifiying variable has been set

      public:
         // class functions
         // constructor
         set_double_t() : value(0.0), setf(false) { }

         // setting function
         void set(double in_value){
            value = in_value;
            setf = true;
         };

         // get value function
         double get(){ return value; };
         // check if variable is set
         bool is_set(){ return setf; };

      };
      //-----------------------------------------------------------------------------
      // materials class for storing exchange material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

         private:

         public:

            //----------------
            // variables
            //----------------
            set_double_t spin_mix_conductance; // spin mixing conductance (1 / Ohm m^2 )

            // constructor
            mp_t (const unsigned int max_materials = 100) {
               spin_mix_conductance.set(0.0); // default value is for generic non magnetic material
            }; // end of constructor

      }; // end of exchange::internal::mp class

      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
      extern bool enabled; // bool to enable spin transport calculation
      extern bool output_atomistic_spin_pumping_flag; // flag to toggle output of atomic resolution spin current
      extern bool output_cells_spin_pumping_flag; // flag to toggle output of cells-reseolved spin current

      extern unsigned int update_rate;  // number of timesteps between updates
      extern unsigned int time_counter; // number of timesteps since last update
      extern uint64_t config_counter; // counter for update of spin pumping configs to file

      extern std::vector<internal::mp_t> mp; // array of material properties

      extern double cell_size_x; // cell size along x-direction
      extern double cell_size_y; // cell size along y-direction
      extern double cell_size_z; // cell size along z-direction

      extern unsigned int total_num_cells; // number of cells

      // arrays to store cell properties
      extern std::vector <bool> cell_magnetic;               // boolean array to determine if cell is magnetic or not
      extern std::vector <double> cell_alpha;                // cell magnetization (average of constituent atoms)
      extern std::vector <double> cell_magnetization;        // 3N normalised magnetization in each cell
      extern std::vector <double> cell_isaturation;          // inverse magnetic saturation at T=0 in each cell
      extern std::vector <double> cell_position;             // 3N array of cell positions (origin)
      extern std::vector <double> cell_spin_mix_conductance; // array to store spin mixing conductance in each cell

      // array to store which cell each atom is in
      extern std::vector <unsigned int> atom_in_cell;

      extern std::vector <bool> material_magnetic;                      // boolean array to determine if material is magnetic or not
      extern std::vector <int> atoms_type_array;                    // array to store material type of atoms

      extern std::vector <double> x_atom_spin_pumping_array;        // arrays to store atom spin pumping x-component
      extern std::vector <double> y_atom_spin_pumping_array;        // arrays to store atom spin pumping y-component
      extern std::vector <double> z_atom_spin_pumping_array;        // arrays to store atom spin pumping z-component

      extern std::vector <double> x_cell_spin_pumping_array;        // arrays to store cell spin pumping x-component
      extern std::vector <double> y_cell_spin_pumping_array;        // arrays to store cell spin pumping y-component
      extern std::vector <double> z_cell_spin_pumping_array;        // arrays to store cell spin pumping z-component

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Functions to exctract spin array belonging to previous time step
      //-------------------------------------------------------------------------
      std::vector<double> get_old_spins_x(const unsigned num_atoms);
      std::vector<double> get_old_spins_y(const unsigned num_atoms);
      std::vector<double> get_old_spins_z(const unsigned num_atoms);

      //---------------------------------------------------------------------------
      // Function to compute atomistic spin pumping as : s_i x ds_i/dt
      //---------------------------------------------------------------------------
      void calculate_spin_pumping(const unsigned int num_local_atoms,            // number of local atoms
                           const std::vector<double>& atoms_x_spin_array, // x-spin vector of atoms
                           const std::vector<double>& atoms_y_spin_array, // y-spin vector of atoms
                           const std::vector<double>& atoms_z_spin_array, // z-spin-vector of atoms
                           const std::vector<double>& atoms_x_old_spin_array, // old x-spin vector of atoms
                           const std::vector<double>& atoms_y_old_spin_array, // old y-spin vector of atoms
                           const std::vector<double>& atoms_z_old_spin_array, // old z-spin-vector of atoms
                           const std::vector<double>& atoms_m_spin_array  // moment of atoms
      );
 
      //---------------------------------------------------------------------------
      // Function to calculate cells magnetisation
      //---------------------------------------------------------------------------
      void calculate_cell_magnetization(const unsigned int num_local_atoms,            // number of local atoms
                           const std::vector<double>& atoms_x_spin_array, // x-spin vector of atoms
                           const std::vector<double>& atoms_y_spin_array, // y-spin vector of atoms
                           const std::vector<double>& atoms_z_spin_array, // z-spin-vector of atoms
                           const std::vector<double>& atoms_m_spin_array  // moment of atoms
      );

      //---------------------------------------------------------------------------
      // Function to compute cells spin pumping 
      //---------------------------------------------------------------------------
      void calculate_cells_spin_pumping();

      //-----------------------------------------------------------------
      // Function to initialise atomic resolution output of coordinates
      //-----------------------------------------------------------------
      void output_atomistic_coordinates(const int num_atoms,                      // number of atoms (only correct in serial)
                                        const std::vector<double>& x_coord_array, // atomic coordinates (angstroms)
                                        const std::vector<double>& y_coord_array,
                                        const std::vector<double>& z_coord_array,
                                        const std::vector<double>& moments_array);

      //-----------------------------------------------------------------
      // Function to output atomic resolution dipole field
      //-----------------------------------------------------------------
      void output_atomistic_spin_pumping(const uint64_t config_file_counter);

      //-----------------------------------------------------------------
      // Function to output atomic resolution dipole field
      //-----------------------------------------------------------------
      void output_cells_spin_pumping(const uint64_t config_file_counter);

   } // end of internal namespace

} // end of spin_transport namespace

#endif //SPINPUMPING_INTERNAL_H_
