//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

#ifndef DIPOLE_INTERNAL_H_
#define DIPOLE_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the dipole module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "dipole.hpp"
#ifdef FFT
#include <fftw3.h>
#endif
// dipole module headers
#include "internal.hpp"

namespace dipole{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
      extern bool initialised;
      extern bool output_atomistic_dipole_field; // flag to toggle output of atomic resolution dipole field

      // enumerated list of different dipole solvers
      enum solver_t{
         macrocell      = 0, // original bare macrocell method (cheap but inaccurate)
         tensor         = 1, // new macrocell with tensor including local corrections
         //multipole    = 2, // bare macrocell but with multipole expansion
         hierarchical   = 3, // new macrocell with tensor including local corrections and nearfield multipole
         atomistic      = 4, // atomistic dipole dipole (too slow for anything over 1000 atoms)
         fft            = 5, // fft method wit tranlational invariance
         atomisticfft   = 6   // atomistic dipole dipole with fft
      };
      extern std::vector < int > cell_dx;
      extern std::vector < int > cell_dy;
      extern std::vector < int > cell_dz;

      extern std::vector < std::vector < std::vector<int> > > idarray;

      extern solver_t solver;

      extern int update_time; /// last update time

      extern const double prefactor; // 1e-7/1e30

      extern std::vector <std::vector < double > > rij_tensor_xx;
      extern std::vector <std::vector < double > > rij_tensor_xy;
      extern std::vector <std::vector < double > > rij_tensor_xz;

      extern std::vector <std::vector < double > > rij_tensor_yy;
      extern std::vector <std::vector < double > > rij_tensor_yz;
      extern std::vector <std::vector < double > > rij_tensor_zz;

      extern int num_atoms;
      extern std::vector < int > atom_type_array;
      extern std::vector < int > atom_cell_id_array;

      extern int cells_num_cells;
      extern int cells_num_local_cells;
      extern std::vector <int>  cells_local_cell_array;
      extern std::vector <int>  cells_num_atoms_in_cell;
      extern std::vector < double > cells_volume_array;

      extern std::vector<double> cells_pos_and_mom_array;
      extern std::vector < int > proc_cell_index_array1D;

      //------------------------------------------------------------------------
      // data structures for atomistic solver
      // (copy of all atom positions and spins on all processors)
      //------------------------------------------------------------------------

      extern int num_local_atoms; // number of local atoms (my processor)
      extern int total_num_atoms; // number of total atoms (all processors)

      // arrays to store atomic coordinates
      extern std::vector <double> cx;
      extern std::vector <double> cy;
      extern std::vector <double> cz;

      // arrays to store atomic spins
      extern std::vector <double> sx;
      extern std::vector <double> sy;
      extern std::vector <double> sz;
      extern std::vector <double> sm;

      // arrays for calculating displacements for parallelisation
      extern std::vector <int> receive_counts;
      extern std::vector <int> receive_displacements;

      //FFT arrays

      #ifdef FFT
         extern fftw_plan MxP,MyP,MzP;
         extern fftw_plan HxP,HyP,HzP;
         extern fftw_complex *N2xx0; //3D Array for dipolar field
         extern fftw_complex *N2xy0;
         extern fftw_complex *N2xz0;

         extern fftw_complex *N2yx0; //3D Array for dipolar field
         extern fftw_complex *N2yy0;
         extern fftw_complex *N2yz0;

         extern fftw_complex *N2zx0; //3D Array for dipolar field
         extern fftw_complex *N2zy0;
         extern fftw_complex *N2zz0;

         extern fftw_complex *N2xx; //3D Array for dipolar field
         extern fftw_complex *N2xy;
         extern fftw_complex *N2xz;

         extern fftw_complex *N2yx; //3D Array for dipolar field
         extern fftw_complex *N2yy;
         extern fftw_complex *N2yz;

         extern fftw_complex *N2zx; //3D Array for dipolar field
         extern fftw_complex *N2zy;
         extern fftw_complex *N2zz;

         extern fftw_complex *Mx_in; //3D Array for dipolar field
         extern fftw_complex *My_in;
         extern fftw_complex *Mz_in;

         extern fftw_complex *Hx_in; //3D Array for dipolar field
         extern fftw_complex *Hy_in;
         extern fftw_complex *Hz_in;

         extern fftw_complex *Mx_out; //3D Array for dipolar field
         extern fftw_complex *My_out;
         extern fftw_complex *Mz_out;

         extern fftw_complex *Hx_out; //3D Array for dipolar field
         extern fftw_complex *Hy_out;
         extern fftw_complex *Hz_out;

         extern unsigned int num_macro_cells_x;
         extern unsigned int num_macro_cells_y;
         extern unsigned int num_macro_cells_z;
         extern unsigned int eight_num_cells;
      #endif

      extern void update_field_fft();
      void initialize_fft_solver();

      namespace atomistic_fft{
          void initialize_atomistic_fft_solver();
          void update_field_atomistic_fft();
          void finalize_atomistic_fft_solver();
      }


      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      //void write_macrocell_data();
      int hierarchical_mag();
      extern void update_field();

      void allocate_memory(const int cells_num_local_cells, const int cells_num_cells);

      void initialize_tensor_solver(const int cells_num_atoms_in_unit_cell,
                                    int cells_num_cells, /// number of macrocells
                                    int cells_num_local_cells, /// number of local macrocells
                                    const double cells_macro_cell_size,
                                    std::vector <int>& cells_local_cell_array,
                                    std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                                    std::vector <int>& cells_num_atoms_in_cell_global, /// number of atoms in each cell
                                    std::vector < std::vector <int> >& cells_index_atoms_array,
                                    std::vector<double>& cells_volume_array,
                                    std::vector<double>& cells_pos_and_mom_array, // array to store positions and moment of cells
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,
                                    std::vector<int>& atom_type_array,
                                    std::vector<int>& atom_cell_id_array,
                                    std::vector<double>& atom_coords_x, //atomic coordinates
                                    std::vector<double>& atom_coords_y,
                                    std::vector<double>& atom_coords_z,
                                    int num_atoms);

      // new version of inter tensor method
      void compute_inter_tensor(const int celli,                                                // global ID of cell i
                                const int cellj,                                                // global ID of cell i
                                const int lc,                                                   // local cell ID
                                const double cutoff,                                            // cutoff range for dipole tensor construction (Angstroms)
                                const std::vector<int>& global_atoms_in_cell_count,             // number of atoms in each cell (all CPUs)
                                const std::vector<double>& cells_pos_and_mom_array,             // array of positions and cell moments
                                const std::vector<int>& list_of_cells_with_atoms,               // list of cells to access atoms
                                const std::vector< std::vector<double> >& atoms_in_cells_array  // output array of positions and moments of atoms in cells
                               );

      void compute_intra_tensor(const int celli,                                                // global ID of cell i
                                const int cellj,                                                // global ID of cell i
                                const int lc,                                                   // local cell ID
                                const std::vector<int>& global_atoms_in_cell_count,             // number of atoms in each cell (all CPUs)
                                const std::vector<int>& list_of_cells_with_atoms,               // list of cells to access atoms
                                const std::vector< std::vector<double> >& atoms_in_cells_array  // output array of positions and moments of atoms in cells
                               );

      void initialize_macrocell_solver(const int cells_num_atoms_in_unit_cell,
                                       int cells_num_cells, /// number of macrocells
                                       int cells_num_local_cells, /// number of local macrocells
                                       const double cells_macro_cell_size,
                                       std::vector <int>& cells_local_cell_array,
                                       std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                                       std::vector <int>& cells_num_atoms_in_cell_global, /// number of atoms in each cell
                                       std::vector < std::vector <int> >& cells_index_atoms_array,
                                       std::vector<double>& cells_volume_array,
                                       std::vector<double>& cells_pos_and_mom_array, // array to store positions and moment of cells
                                       std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                                       std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                                       std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,
                                       std::vector<int>& atom_type_array,
                                       std::vector<int>& atom_cell_id_array,
                                       std::vector<double>& atom_coords_x, //atomic coordinates
                                       std::vector<double>& atom_coords_y,
                                       std::vector<double>& atom_coords_z,
                                       int num_atoms);

      void initialize_atomistic_solver(int num_atoms,                      // number of atoms (only correct in serial)
                                       std::vector<double>& x_coord_array, // atomic corrdinates (angstroms)
                                       std::vector<double>& y_coord_array,
                                       std::vector<double>& z_coord_array,
                                       std::vector<double>& moments_array, // atomistic magnetic moments (bohr magnetons)
                                       std::vector<int>& mat_id_array);    // atom material ID

      void calculate_atomistic_dipole_field(std::vector<double>& x_spin_array, // atomic spin directions
                                            std::vector<double>& y_spin_array,
                                            std::vector<double>& z_spin_array);

      //-----------------------------------------------------------------------------
      // Function to send receive cells data to other cpus
      //-----------------------------------------------------------------------------
      int send_recv_cells_data(std::vector<int>& proc_cell_index_array1D,
                               std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_x,
                               std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_y,
                               std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_z,
                               std::vector< std::vector <int> >& cells_index_atoms_array,
                               std::vector<double>& cells_pos_and_mom_array,
                               std::vector<int>& cells_num_atoms_in_cell,
                               std::vector<int>& cells_cell_id_array,
                               std::vector<int>& cells_local_cell_array,
                               int cells_num_local_cells,
                               int cells_num_cells);

      //-----------------------------------------------------------------------------
      // Function to send receive atoms data to other cpus
      //-----------------------------------------------------------------------------
      int send_recv_atoms_data(std::vector<int>& proc_cell_index_array2D,
                               std::vector<int>& cell_id_array,
                               std::vector<int>& cells_local_cell_array,
                               std::vector<double>& atom_pos_x,
                               std::vector<double>& atom_pos_y,
                               std::vector<double>& atom_pos_z,
                               std::vector<int>& atom_type_array,
                               std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_x,
                               std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_y,
                               std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_z,
                               std::vector< std::vector <int> >& cells_index_atoms_array,
                               std::vector<double>& cells_pos_and_mom_array,
                               std::vector<int>& cells_num_atoms_in_cell,
                               int cells_num_local_cells,
                               int cells_num_cells,
                               double cells_macro_cell_size);

      //----------------------------------------------------------------
      //Function to sort cells/atoms data after sharing
      //----------------------------------------------------------------
      int sort_data(std::vector<int>& proc_cell_index_array1D,
                  std::vector<int>& cells_cell_id_array,
                  std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_x,
                  std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_y,
                  std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_z,
                  std::vector< std::vector <int> >& cells_index_atoms_array,
                  std::vector<double>& cells_pos_and_mom_array,
                  std::vector<int>& cells_num_atoms_in_cell,
                  int cells_num_local_cells,
                  int cells_num_cells);

      void initialise_atomistic_cell_data(const int num_cells,
                                          const int num_local_cells,
                                          const double cutoff,                                         // cutoff range for dipole tensor construction (Angstroms)
                                          const std::vector<int>& num_atoms_in_cell,                   // number of atoms in each cell (local CPU)
                                          const std::vector<int>& list_of_local_cells,                 // numerical list of cells containing atoms on local processor
                                          const std::vector<int>& global_atoms_in_cell_count,          // number of atoms in each cell (all CPUs)
                                          const std::vector<double>& pos_and_mom_array,                // array of positions and cell moments
                                          const std::vector < std::vector <int> >& index_atoms_array,  // 2D array of [cells][atomID]
                                          const std::vector<double>& atoms_coords_x,                   // input arrays of atom coordinates
                                          const std::vector<double>& atoms_coords_y,                   //
                                          const std::vector<double>& atoms_coords_z,                   //
                                          const std::vector<double>& atoms_moments,                    // input array of atom moments (Bohr magnetons)
                                          std::vector<int>& list_of_cells_with_atoms,                  // list of cells to access atoms
                                          std::vector< std::vector<double> >& atoms_in_cells_array     // output array of positions and moments of atoms in cells
                                         );

      /*--------------------------------------------------------*/
      /*Function to send cells field to be output in cfg file   */
      /*--------------------------------------------------------*/
      int send_cells_demag_factor(std::vector<int>& cells_cell_id_array,
                                 std::vector<double>& N_tensor_array,
                                 int cells_num_local_cells);

      //-----------------------------------------------------------------
      // Function to initialise atomic resolution output of dipole field
      //-----------------------------------------------------------------
      void output_atomistic_coordinates(int num_atoms,                      // number of atoms (only correct in serial)
                                        std::vector<double>& x_coord_array, // atomic coordinates (angstroms)
                                        std::vector<double>& y_coord_array,
                                        std::vector<double>& z_coord_array,
                                        std::vector<double>& moments_array);

       //-----------------------------------------------------------------
       // Function to output atomic resolution dipole field
       //-----------------------------------------------------------------
       void output_atomistic_dipole_fields();

       //-----------------------------------------------------------------
       // Function to output solver memory information
       //-----------------------------------------------------------------
       void output_dipole_solver_mem_info(int num_cells, int num_local_cells);

   } // end of internal namespace

} // end of dipole namespace

#endif //DIPOLE_INTERNAL_H_
