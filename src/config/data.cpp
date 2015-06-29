//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "config.hpp"

// config headers
#include "internal.hpp"

namespace config{

   //-----------------------------------------------------------------------------------------------
   // Externally visible variables
   //-----------------------------------------------------------------------------------------------

   namespace internal{
      //-----------------------------------------------------------------------------
      // Shared variables used for configuration output
      //-----------------------------------------------------------------------------
      bool output_atoms = false; // flag to control atomic output
      bool output_cells = false; // flag to control macrocell output
      bool output_meta = false; // flag to control meta data output
      bool output_coords = false; // flag to control atomic coordinate output
      int output_rate = 1000; // relative rate of file output
      int output_file_counter = 0; // configuration file number
      int output_rate_counter = 0; // configuration output counter

      data_format_t output_data_format = text;

      double atoms_output_min[3] = {0.0,0.0,0.0}; // Spatial range fr atomic output
      double atoms_output_max[3] = {1.0,1.0,1.0};
      std::vector<int> local_output_atom_list(0); // list of atoms to be output according to spatial range      
      int total_output_atoms = 0; // Total number of atoms to be output on local node

      // Buffer variables to store copies of data in float format for reduced file size
      std::vector<float> output_spin_buffer(0); // buffer to store float cast spin array for output to disk

      namespace mpi{         
         std::vector<float> storage_buffer(0); // temporary buffer to store output in parallel version
         int buffer_x_offset = 0; // offsets for mpi buffer locations
         int buffer_y_offset = 0;
         int buffer_z_offset = 0;

         int my_io_rank = 0; // my id in i/o communicator
         int my_io_master = 0; // my master id in i/o communicator who I send data to
         bool io_master = true; // flag set to true if I output data
         int num_io_nodes = 0; // total number of i/o processes
      } // end of mpi namespace

   } // end of internal namespace
} // end of config namespace

