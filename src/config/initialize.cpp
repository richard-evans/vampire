//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iomanip>
#include <sstream>

// Vampire headers
#include "config.hpp"
#include "material.hpp"
#include "vio.hpp"

// config headers
#include "internal.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions to output configuration files
//--------------------------------------------------------------------------------
namespace config{

   // Forward declaration of functions


   //-----------------------------------------------------------------------------
   //
   //
   //
   //
   //
   //
   //
   //
   //
   //-----------------------------------------------------------------------------
   void initialize(const int num_atoms, // number of local atoms
                   const double system_dimensions_x, // system size
                   const double system_dimensions_y,
                   const double system_dimensions_z,
                   const std::vector<int>& material, // material id
                   const std::vector<int>& category, // category id
                   const std::vector<double>& spins_cx, // spin coordinates (Angstroms)
                   const std::vector<double>& spins_cy,
                   const std::vector<double>& spins_cz,
                   const std::vector<double>& cells_cx, // cell coordinates (Angstroms)
                   const std::vector<double>& cells_cy,
                   const std::vector<double>& cells_cz)
   {

      // Check if config output is needed
      if(config::internal::output_meta || config::internal::output_coords){

      zlog << zTs() << "Initializing configuration file output" << std::endl;

      //--------------------------------------------------------------------
      // Determine atoms to output
      //--------------------------------------------------------------------

      // resize atom list to zero
      config::internal::local_output_atom_list.resize(0);

      // get output bounds
      const double minB[3]={config::internal::atoms_output_min[0]*system_dimensions_x,
                            config::internal::atoms_output_min[1]*system_dimensions_y,
                            config::internal::atoms_output_min[2]*system_dimensions_z};

      const double maxB[3]={config::internal::atoms_output_max[0]*system_dimensions_x,
                            config::internal::atoms_output_max[1]*system_dimensions_y,
                            config::internal::atoms_output_max[2]*system_dimensions_z};

      // loop over all local atoms and record output list
      for(int atom=0;atom<num_atoms;atom++){

         const double cc[3] = {spins_cx[atom],spins_cy[atom],spins_cz[atom]};

         // check atom within output bounds
         if((cc[0] >= minB[0]) && (cc[0]<=maxB[0])){
            if((cc[1] >= minB[1]) && (cc[1]<=maxB[1])){
               if((cc[2] >= minB[2]) && (cc[2]<=maxB[2])){
                  config::internal::local_output_atom_list.push_back(atom);
               }
            }
         }
      }

      // initialise mpi here
      config::internal::mpi::initialize();

      // resize output buffer to 3*num_output_atoms
      config::internal::output_spin_buffer.resize(3*config::internal::total_output_atoms);

      // cell buffer...

      //-------------------------------------------------------
      // Output spin coordinate meta data
      //-------------------------------------------------------
      if(config::internal::output_atoms || config::internal::output_coords){

         // write coordinate meta data
         config::internal::write_coordinate_meta();

         // write coordinate and id data
         config::internal::write_coordinate_data(spins_cx, spins_cy, spins_cz, material, category);

      }

      }

      return;

   }



} // end of namespace config
