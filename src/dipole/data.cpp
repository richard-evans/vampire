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
   bool fft = false;
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

      fftw_plan MxP,MyP,MzP;
      fftw_plan HxP,HyP,HzP;

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

      //only initialise if FFT has been enabled in makefile
      #ifdef FFT
      fftw_complex *N2xx0; //3D Array for dipolar field
      fftw_complex *N2xy0;
      fftw_complex *N2xz0;

      fftw_complex *N2yx0; //3D Array for dipolar field
      fftw_complex *N2yy0;
      fftw_complex *N2yz0;

      fftw_complex *N2zx0; //3D Array for dipolar field
      fftw_complex *N2zy0;
      fftw_complex *N2zz0;

      fftw_complex *N2xx; //3D Array for dipolar field
      fftw_complex *N2xy;
      fftw_complex *N2xz;

      fftw_complex *N2yx; //3D Array for dipolar field
      fftw_complex *N2yy;
      fftw_complex *N2yz;

      fftw_complex *N2zx; //3D Array for dipolar field
      fftw_complex *N2zy;
      fftw_complex *N2zz;


      fftw_complex *Mx_in; //3D Array for dipolar field
      fftw_complex *My_in;
      fftw_complex *Mz_in;

      fftw_complex *Hx_in; //3D Array for dipolar field
      fftw_complex *Hy_in;
      fftw_complex *Hz_in;

      fftw_complex *Mx_out; //3D Array for dipolar field
      fftw_complex *My_out;
      fftw_complex *Mz_out;

      fftw_complex *Hx_out; //3D Array for dipolar field
      fftw_complex *Hy_out;
      fftw_complex *Hz_out;

      //stores number of macrocells in x,y,z
      unsigned int num_macro_cells_x;
      unsigned int num_macro_cells_y;
      unsigned int num_macro_cells_z;
      unsigned int eight_num_cells;
      #endif

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
      void update_field_fft();

   } // end of internal namespace

} // end of dipole namespace
