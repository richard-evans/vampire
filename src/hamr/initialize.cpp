//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo 2022. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "errors.hpp"
#include "hamr.hpp"
#include "vio.hpp"

// hamr headers
#include "internal.hpp"

namespace hamr{

   //-----------------------------------------------------------------------------------------------
   // Externally visible variables
   //-----------------------------------------------------------------------------------------------

   //---------------------------------------------------------------------------------
   // Function for initialising hamr data structures and variables
   //---------------------------------------------------------------------------------
   void initialize(const double Hmin,
                   const double Hmax,
                   const double system_dimensions[3],
                   const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z,
                   const std::vector<int>& atom_type_array,
                   const int num_local_atoms
                  ){
   
      // output informative message
      zlog << zTs() << "Initialising data structures for hamr calculation." << std::endl;

      // check for prior initialisation
      if(hamr::internal::initialised){
         zlog << zTs() << "Warning: hamr module already initialised. Continuing." << std::endl;
         return;
      }

      //-------------------------------------------------------
      // Save value of local num atoms and resize field arrays
      //-------------------------------------------------------
      hamr::internal::num_local_atoms = num_local_atoms;
      // hamr::internal::x_field_array.resize(num_local_atoms); // arrays to store atomic spin torque field
      // hamr::internal::y_field_array.resize(num_local_atoms);
      // hamr::internal::z_field_array.resize(num_local_atoms);

      // hamr::internal::system_dimensions = system_dimensions;
      hamr::internal::system_dimensions[0] = system_dimensions[0];
      hamr::internal::system_dimensions[1] = system_dimensions[1];
      hamr::internal::system_dimensions[2] = system_dimensions[2];
      hamr::internal::atom_coords_x = atom_coords_x;
      hamr::internal::atom_coords_y = atom_coords_y;
      hamr::internal::atom_coords_z = atom_coords_z;

      hamr::internal::atom_type_array = atom_type_array;
      hamr::internal::Hmin = Hmin;
      hamr::internal::Hmax = Hmax;

      // Set initialised flag
      hamr::internal::initialised = true;

      return;
   }
   
} // end of hamr namespace
