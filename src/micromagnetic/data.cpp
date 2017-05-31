
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

   std::vector < double > list_of_atomistic_atoms;
   std::vector < double > list_of_none_atomistic_atoms;
   std::vector < double > list_of_micromagnetic_cells;

   int number_of_atomistic_atoms = 0;
   int number_of_none_atomistic_atoms = 0;
   int number_of_micromagnetic_cells = 0;


   std::vector < bool > cell_discretisation_micromagnetic;

   int num_atomic_steps_mm = 1;

   int integrator = 1;

   namespace internal{


      //------------------------------------------------------------------------
      // Shared variables inside micromagnetic module
      //------------------------------------------------------------------------


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

      //stores the external fields (x,y,z)
      std::vector<double> ext_field;


      std::vector<double> fields_neighbouring_atoms_begin;
      std::vector<double> fields_neighbouring_atoms_end;

      //macrocell neighbourlists
      std::vector<double> macro_neighbour_list_start_index;
      std::vector<double> macro_neighbour_list_end_index;
      std::vector<double> macro_neighbour_list_array;


   } // end of internal namespace

} // end of micromagnetic namespace
