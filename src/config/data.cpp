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

// C++ standard library headers

// Vampire headers
#include "config.hpp"

// config module headers
#include "internal.hpp"

namespace config{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside config module
      //------------------------------------------------------------------------

      // interface and selection variables
      format_t format = text; // format for data output (text, binary)
      mode_t mode = fpnode; // output mode (legacy, mpi_io, file per process, file per io node)

      bool initialised = false; // flag to signify if config has been initialised

      bool output_atoms_config = false; // flag to enable atoms output
      bool output_atoms_config_continuous = false; // flag to enable continuous output of atomic configurations
      bool output_atoms_config_end = false; // flag to enable atoms output at the end of simulation
      int output_atoms_config_rate = 1000; // rate to output atoms

      bool output_cells_config = false; // flag to enable cells output
      bool output_cells_config_continuous = false; // flag to enable cells output
      bool output_cells_config_end = false; // flag to enable cells output at the end of simulation
      int output_cells_config_rate = 1000; // rate to output cells

      int output_rate_counter_coords = 0;

      bool identify_surface_atoms = false; // flag to identify surface atoms in config coordinate file

      // Field ranges for hysteresis ouput
      double field_output_min_1 = -10000.0;
      double field_output_max_1 = 10000.0;
      double field_output_min_2 = -10000.0;
      double field_output_max_2 = 10000.0;

      // fraction ranges for atom output
      double atoms_output_min[3]={ 0.0, 0.0, 0.0};
      double atoms_output_max[3]={ 1.0, 1.0, 1.0};

      // implementation variables
      std::vector<uint64_t> local_output_atom_list(0); // list of atom numbers to output to disk

      uint64_t total_output_atoms = 0; // total number of atoms to be outputted (all processors)
      uint64_t total_output_cells = 0; // total number of cells to be outputted

      // Data buffers for parallel i/o
      std::vector<double> local_buffer(0);
      std::vector<double> collated_buffer(0);

      // variables for collated data output
      int num_io_groups = 1; // number of processors to output data
      int io_group_size = 1; // number of processors in my io_comm group
      int io_group_id = 0; // the group id to which I belong
      int io_group_rank = 0; // rank of procesor in io_comm group
      int io_group_master_id = 0; // rank of processor in io_comm group responsible for data output
      bool io_group_master = false; // flag to specify master process
      double io_data_size = 0.0; // data size outputted to disk in GB (for binary mode)
      std::vector<int> io_group_recv_counts(0); // data to receive from each process in io group
      std::vector<int> io_group_displacements(0); // offsets in obuf to receive from each process in io group

      #ifdef MPICF
         MPI_Offset linear_offset; // offset for mpi-io collective routines for integer data (bytes)
         MPI_Offset buffer_offset; // offset for mpi-io collective routines for 3 vector double data (bytes)
         MPI_Comm io_comm; // MPI IO communicator specifying a group of processors who output as a group
      #endif


   } // end of internal namespace

} // end of config namespace
