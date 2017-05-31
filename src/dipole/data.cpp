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


   std::vector < double > cells_field_array_x;
   std::vector < double > cells_field_array_y;
   std::vector < double > cells_field_array_z;
   std::vector < double > atom_dipolar_field_array_x;
   std::vector < double > atom_dipolar_field_array_y;
   std::vector < double > atom_dipolar_field_array_z;

   double cutoff = 12.0; //12.0; /// cutoff distance between cells over which bare macro cell model can be applied

   //uint64_t sim_time;

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside dipole module
      //------------------------------------------------------------------------
      bool initialised=false;


      int update_time=-1; /// last update time

      const double prefactor=1.0e+23; // 1e-7/1e30

      std::vector <std::vector < double > > rij_inter_xx;
      std::vector <std::vector < double > > rij_inter_xy;
      std::vector <std::vector < double > > rij_inter_xz;

      std::vector <std::vector < double > > rij_inter_yy;
      std::vector <std::vector < double > > rij_inter_yz;
      std::vector <std::vector < double > > rij_inter_zz;

      std::vector <std::vector < double > > rij_intra_xx;
      std::vector <std::vector < double > > rij_intra_xy;
      std::vector <std::vector < double > > rij_intra_xz;

      std::vector <std::vector < double > > rij_intra_yy;
      std::vector <std::vector < double > > rij_intra_yz;
      std::vector <std::vector < double > > rij_intra_zz;

      Array3D<fftw_complex> Nxx0; //3D Array for dipolar field
	   Array3D<fftw_complex> Nxy0;
	   Array3D<fftw_complex> Nxz0;

		Array3D<fftw_complex> Nyx0; //3D Array for dipolar field
		Array3D<fftw_complex> Nyy0;
		Array3D<fftw_complex> Nyz0;

		Array3D<fftw_complex> Nzx0; //3D Array for dipolar field
		Array3D<fftw_complex> Nzy0;
		Array3D<fftw_complex> Nzz0;

      Array3D<fftw_complex> Nxx; //3D Array for dipolar field
	   Array3D<fftw_complex> Nxy;
	   Array3D<fftw_complex> Nxz;

		Array3D<fftw_complex> Nyx; //3D Array for dipolar field
		Array3D<fftw_complex> Nyy;
		Array3D<fftw_complex> Nyz;

		Array3D<fftw_complex> Nzx; //3D Array for dipolar field
		Array3D<fftw_complex> Nzy;
		Array3D<fftw_complex> Nzz;

      int num_macro_cells_x;
      int num_macro_cells_y;
      int num_macro_cells_z;
      int eight_num_cells;


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
