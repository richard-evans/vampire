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

   //boolean to determine whether the simulation is micromagnetic
   int discretisation_type = 0;

   double MR_resistance = 0.0;

   bool enable_resistance = false;
   //lsits to store atomistic/microamgnetic cells/atoms
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

      double bias_magnets_max_height = 1.0;
      double bias_magnets_min_height= 0.0;
      double bias_magnet_ms_input = 1.0;

      double bias_magnets_max_width = 1.0;
      double bias_magnets_min_width = 0.0;

      int bias_magnets_gap = 10;
      double overlap_area = 0.0;

      double res_GMR = 1.0;
      double res_RA = 0.3;
      double shield_Ms = -1;
      int resistance_layer_1 = 0;
      int resistance_layer_2 = 0;

      std::vector <double> fmr_H(3,0.0);

      std::vector < double > bias_field_x;
      std::vector < double > bias_field_y;
      std::vector < double > bias_field_z;

      bool bias_magnets = false;
      //stores the micromagnetic properties of the macrocells
      std::vector<double> A;
      std::vector<double> alpha;
      std::vector<double> one_o_chi_perp;
      std::vector<double> one_o_chi_para;
      std::vector<double> gamma;
      std::vector<double> ku;
      std::vector<double> ku_x;
      std::vector<double> ku_y;
      std::vector<double> ku_z;
      std::vector<double> ms;
      std::vector<double> T;
      std::vector<double> Tc;
      std::vector<double> m_e;
      std::vector<double> alpha_para;
      std::vector<double> alpha_perp;

      std::vector <double> mat_vol;
      std::vector <double> mat_ms;
      std::vector <double> prefactor;

      std::vector <double> cell_material_array;

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
