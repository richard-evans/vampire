//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) rory.pond 2016. All rights reserved.
//
//   Email: rory.pond@york.ac.uk
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
        std::vector<int> local_output_atom_list(0);
      //------------------------------------------------------------------------
      // Shared variables inside config module
      //------------------------------------------------------------------------
      bool output_atoms_config=false;
      int output_atoms_config_rate=1000;

      bool output_cells_config=false;
      int output_cells_config_rate=1000;

      double field_output_min_1=-10000.0;
      double field_output_max_1=-0.0;
      double field_output_min_2=0.0;
      double field_output_max_2=10000.0;

      double atoms_output_min[3]={0.0,0.0,0.0};
      double atoms_output_max[3]={1.0,1.0,1.0};

      data_format output_data_format = text;
      bool output_gather = false;
      bool output_all = false;
      bool output_mpi_io = false;
      bool output_lagacy =true;
   } // end of internal namespace

} // end of config namespace

