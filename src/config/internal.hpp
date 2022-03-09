//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Rory Pond and Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef CONFIG_INTERNAL_H_
#define CONFIG_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the config module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <cstdint>

// Vampire headers
#include "config.hpp"
#include "vmpi.hpp"

namespace config
{

namespace internal
{

   // enumerated integers for option selection
   enum format_t{ binary = 0, text = 1};
   enum mode_t{ legacy = 0, mpi_io = 1, fpprocess = 2, fpnode = 3};

   //-------------------------------------------------------------------------
   // Internal data type definitions
   //-------------------------------------------------------------------------

   extern format_t format; // format for data output (text, binary)
   extern mode_t mode; // output mode (legacy, mpi_io, file per process, file per io node)

   extern bool initialised; // flag to signify if config has been initialised

   extern bool output_atoms_config; // flag to enable atoms output
   extern bool output_atoms_config_continuous; // flag to enable continuous output of atomic configurations
   extern bool output_atoms_config_end; // flag to enable output of atomic configurations at end of simulation
   extern int output_atoms_config_rate; // rate to output atoms

   extern bool output_cells_config; // flag to enable cells output
   extern bool output_cells_config_continuous; // flag to enable continuous output of cells configurations
   extern bool output_cells_config_end; // flag to enable output of cells configurations at end of simulation
   extern int output_cells_config_rate; // rate to output cells

   extern int output_rate_counter_coords;

   extern bool identify_surface_atoms; // flag to identify surface atoms in config coordinate file

   // Field ranges for hysteresis ouput
   extern double field_output_min_1;
   extern double field_output_max_1;
   extern double field_output_min_2;
   extern double field_output_max_2;

   // fraction ranges for atom output
   extern double atoms_output_min[3];
   extern double atoms_output_max[3];

   // implementation variables
   extern std::vector<uint64_t> local_output_atom_list; // list of atom numbers to output to disk

   extern uint64_t total_output_atoms; // total number of atoms to be outputted (all processors)
   extern uint64_t total_output_cells; // total number of cells to be outputted

   // Data buffers for parallel i/o
   extern std::vector<double> local_buffer;
   extern std::vector<double> collated_buffer;

   // variables for collated data output
   extern int num_io_groups; // number of processors to output data
   extern int io_group_size; // number of processors in my io_comm group
   extern int io_group_id; // the group id to which I belong
   extern int io_group_rank; // rank of procesor in io_comm group
   extern int io_group_master_id; // rank of processor in io_comm group responsible for data output
   extern bool io_group_master; // flag to specify master process
   extern double io_data_size; // data size outputted to disk in GB (for binary mode)
   extern std::vector<int> io_group_recv_counts; // data to receive from each process in io group
   extern std::vector<int> io_group_displacements; // offsets in obuf to receive from each process in io group

   #ifdef MPICF
      extern MPI_Offset linear_offset; // offset for mpi-io collective routines for integer data (bytes)
      extern MPI_Offset buffer_offset; // offset for mpi-io collective routines for 3 vector double data (bytes)
      extern MPI_Comm io_comm; // MPI IO communicator specifying a group of processors who output as a group
   #endif

   //-------------------------------------------------------------------------
   // Internal shared variables
   //-------------------------------------------------------------------------

   //-------------------------------------------------------------------------
   // Internal function declarations
   //-------------------------------------------------------------------------
   void initialize();

   void atoms();
   void atoms_coords();
   void atoms_non_magnetic();

   double legacy_atoms();
   double legacy_atoms_coords();
   void legacy_cells();
   void legacy_cells_coords();

   double write_data(std::string, const std::vector<double> &buffer);
   double write_coord_data(std::string filename, const std::vector<double>& buffer, const std::vector<int>& type_buffer, const std::vector<int>& category_buffer);

   void copy_data_to_buffer(const std::vector<double> &x, // vector data
                            const std::vector<double> &y,
                            const std::vector<double> &z,
                            const std::vector<uint64_t> &mask,
                            std::vector<double> &buffer);

   void write_coordinate_meta();
   void write_non_magnetic_meta(const uint64_t num_data);
   void write_meta(const double simulation_time, // time (seconds)
                   const double temperature,     // system temperature (Kelvin)
                   const double applied_field_x, // applied field components (Tesla)
                   const double applied_field_y,
                   const double applied_field_z,
                   const double applied_field_mag,
                   const double magnetization_x, // magnetization components (normalized)
                   const double magnetization_y,
                   const double magnetization_z);

   std::string data_filename(bool coords);

} // end of internal namespace

} // end of config namespace

#endif //CONFIG_INTERNAL_H_
