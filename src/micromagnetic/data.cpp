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

//   bool micro;
      bool discretisation_micromagnetic = false;
      bool initialised = false;
   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   std::vector < std::vector <int > > P;
   std::vector < int > P1D;
   double mean_M;
   int counter;

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside micromagnetic module
      //------------------------------------------------------------------------


      std::vector<double> Ax;
      std::vector<double> Ay;
      std::vector<double> Az;
      std::vector<double> alpha;
      std::vector<double> chi_perp;
      std::vector<double> chi_para;
      std::vector<double> gamma;
      std::vector<double> ku;
      std::vector<double> ms;
      std::vector<double> Tc;

      std::vector<double> x_array;
      std::vector<double> y_array;
      std::vector<double> z_array;

      std::vector<double> ext_field;

      std::vector<double> x_euler_array;
      std::vector<double> y_euler_array;
      std::vector<double> z_euler_array;

      std::vector<double> x_heun_array;
      std::vector<double> y_heun_array;
      std::vector<double> z_heun_array;

      std::vector<double> mx_store;
      std::vector<double> my_store;
      std::vector<double> mz_store;

      std::vector<double> mx_init;
      std::vector<double> my_init;
      std::vector<double> mz_init;

      std::vector<double> macro_neighbour_list_start_index;
      std::vector<double> macro_neighbour_list_end_index;
      std::vector<double> macro_neighbour_list_array;
      /*
   	const double prefactor=1.0e+23; // 1e-7/1e30

      Array3D<fftw_complex> Nxx; // creates the stencil complex array Nxx
      Array3D<fftw_complex> Nyx; // creates the stencil complex array Nyx
      Array3D<fftw_complex> Nzx; // creates the stencil complex array Nzx

      Array3D<fftw_complex> Nxy; // creates the stencil complex array Nxy
      Array3D<fftw_complex> Nyy; // creates the stencil complex array Nyy
      Array3D<fftw_complex> Nzy; // creates the stencil complex array Nzy

      Array3D<fftw_complex> Nxz; // creates the stencil complex array Nxz
      Array3D<fftw_complex> Nyz; // creates the stencil complex array Nyz
      Array3D<fftw_complex> Nzz; // creates the stencil complex array Nzz
*/
   //   int num_macro_cells_x;
   //   int num_macro_cells_y;
   //   int num_macro_cells_z;


   } // end of internal namespace

} // end of micromagnetic namespace
