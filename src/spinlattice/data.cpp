//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "sld.hpp"
#include <vector>
#include "material.hpp"

// sld module headers
#include "internal.hpp"

namespace sld{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------
   double var_test=0;
   double spin_temperature;
   double lattice_temperature;
   double J_eff;
   double C_eff;
   bool suzuki_trotter_parallel_initialized = false;
   
    

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside sld module
      //------------------------------------------------------------------------

      bool enabled; // bool to enable module

      std::vector<internal::mp_t> mp; // array of material properties

      double r_cut_pot; // mechanical potential cutoff
      double r_cut_fields;
      
      double dr_init; // initial conditions
      double th_velo;

      bool harmonic; //flag for harmonic potential
      bool pseudodipolar;

      //initial sld neighbor list
      //std::vector<int> sld_neighbour_list_start_index;
      //std::vector<int> sld_neighbour_list_end_index;
      //std::vector<int> sld_neighbour_list_array;

      std::vector<double> x0_coord_array;
      std::vector<double> y0_coord_array;
      std::vector<double> z0_coord_array;
     
     
      std::vector <double> x_coord_storage_array;
      std::vector <double> y_coord_storage_array;
      std::vector <double> z_coord_storage_array;

      std::vector<double> forces_array_x;
      std::vector<double> forces_array_y;
      std::vector<double> forces_array_z;

      std::vector<double> fields_array_x;
      std::vector<double> fields_array_y;
      std::vector<double> fields_array_z;

      std::vector<double> velo_array_x;
      std::vector<double> velo_array_y;
      std::vector<double> velo_array_z;

      std::vector<double> potential_eng;
      std::vector<double> sumJ;
      std::vector<double> sumC;
      std::vector<double> exch_eng;
      std::vector<double> coupl_eng;
      
      std::vector<int> test_atom_list; //Core atoms of each octant

      
      //MPI variables
      std::vector<std::vector<int> > c_octants; //Core atoms of each octant
      std::vector<std::vector<int> > b_octants; //Boundary atoms of each octant
      std::vector <int> all_atoms_octant_start_index; 
      std::vector <int> all_atoms_octant_end_index;
      std::vector <int> all_atoms_octant;


   } // end of internal namespace

} // end of sld namespace
