
//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------
   bool enabled = false; // Flag to determine if micromagnetic simulation is enabled

   //boolean to determine whether the simulation is micromagnetic
   discretisation_t discretisation_type = atomistic;

   double MR_resistance = 0.0;

   bool enable_resistance = false;
   //lists to store atomistic/microamgnetic cells/atoms
   std::vector <int> list_of_atomistic_atoms(0);
   std::vector <int> list_of_none_atomistic_atoms(0);
   std::vector <int> list_of_micromagnetic_cells(0);
   std::vector <int> list_of_empty_micromagnetic_cells(0);

   //sets initial values to 0
   int number_of_atomistic_atoms = 0;
   int number_of_none_atomistic_atoms = 0;
   int number_of_micromagnetic_cells = 0;

   //is the discretisation of each cell microamgnetic or atomistic
   std::vector < bool > cell_discretisation_micromagnetic;

   //number of atomic steps per micromagnetic step
   int num_atomic_steps_mm = 1;

   //sets the integrator to be LLB if none is specified
   int integrator = 1;


   //variables for boltzman
   std::vector < std::vector < double > > P;

   std::vector < double > P1D(1001,0.0);

   bool boltzman = false;

   double mean_M=0.0;
   double counter=0.0;



   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside micromagnetic module
      //------------------------------------------------------------------------

      int my_num_micromagnetic_cells;
      int my_start_index; // first cell to intergrate on local (my) cpu
      int my_end_index;  // last cell +1 to intergrate on local (my) cpu


      int resistance_layer_1 = 0;
      int resistance_layer_2 = 0;

      std::vector <double> fmr_H(3,0.0);

      //stores the micromagnetic properties of the macrocells
      std::vector<double> A;
      std::vector<double> alpha;
      std::vector<double> one_o_chi_perp;
      std::vector<double> one_o_chi_para;
      std::vector<double> gamma;
      std::vector<double> ku;
      std::vector<double> ms;
      std::vector<double> Tc;
      std::vector<double> m_e;
      std::vector<double> alpha_para;
      std::vector<double> alpha_perp;

      std::vector <double> cell_material_array;

      // Thermal field array
      std::vector<double> thermal_field_array_x(0);
      std::vector<double> thermal_field_array_y(0);
      std::vector<double> thermal_field_array_z(0);

      bool mm_correction;
      double pinning_field_height;
      std::vector <double> pinning_field_x;
      std::vector <double> pinning_field_y;
      std::vector <double> pinning_field_z;

      //stores the external fields (x,y,z)
      std::vector<double> ext_field;

      //start and end index arrays for the neighbouring atoms for field calcualtions.
      std::vector<double> fields_neighbouring_atoms_begin;
      std::vector<double> fields_neighbouring_atoms_end;

      //macrocell neighbourlists
      std::vector<double> macro_neighbour_list_start_index;
      std::vector<double> macro_neighbour_list_end_index;
      std::vector<double> macro_neighbour_list_array;


   } // end of internal namespace

} // end of micromagnetic namespace
