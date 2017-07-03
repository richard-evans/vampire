//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "dipole.hpp"

// dipole module headers
#include "internal.hpp"

namespace dipole{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------
   int update_rate=100; /// timesteps between updates
   bool activated=false;

   // define arrays for B-field
   std::vector < double > cells_field_array_x;
   std::vector < double > cells_field_array_y;
   std::vector < double > cells_field_array_z;
   std::vector < double > atom_dipolar_field_array_x;
   std::vector < double > atom_dipolar_field_array_y;
   std::vector < double > atom_dipolar_field_array_z;
   // define arrays for mu_0*Hdemag - field
   std::vector < double > cells_mu0Hd_field_array_x;
   std::vector < double > cells_mu0Hd_field_array_y;
   std::vector < double > cells_mu0Hd_field_array_z;
   std::vector < double > atom_mu0demag_field_array_x;
   std::vector < double > atom_mu0demag_field_array_y;
   std::vector < double > atom_mu0demag_field_array_z;

   double cutoff = 2.0;  /// cutoff distance between cells over which bare macro cell model can be applied
                         /// N.B.: after 12 cells inter-intra method is equivalent to bare macrocell method.
                         /// Although, 2 cells is enough because there are other error sources limiting the accuracy.

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside dipole module
      //------------------------------------------------------------------------
      bool initialised=false;

      int update_time=-1; /// last update time

      // solver to be used for dipole method
      dipole::internal::solver_t solver = dipole::internal::tensor; // default is tensor method

      const double prefactor=1.0e+23; // 1e-7/1e30

      std::vector <std::vector < double > > rij_tensor_xx;
      std::vector <std::vector < double > > rij_tensor_xy;
      std::vector <std::vector < double > > rij_tensor_xz;

      std::vector <std::vector < double > > rij_tensor_yy;
      std::vector <std::vector < double > > rij_tensor_yz;
      std::vector <std::vector < double > > rij_tensor_zz;

      int num_atoms;
      std::vector < int > atom_type_array;
      std::vector < int > atom_cell_id_array;

      int cells_num_cells;
      int cells_num_local_cells;
      std::vector <int>  cells_local_cell_array;
      std::vector <int>  cells_num_atoms_in_cell;
      std::vector < double > cells_volume_array;

      std::vector<double> cells_pos_and_mom_array;
      std::vector < int > proc_cell_index_array1D;
      //------------------------------------------------------------------------
      // Shared functions inside dipole module
      //------------------------------------------------------------------------
      void update_field();

   } // end of internal namespace

} // end of dipole namespace
