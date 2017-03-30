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

#ifndef CONFIG_INTERNAL_H_
#define CONFIG_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the config module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "config.hpp"

namespace config
{
extern int output_rate_counter_coords;
extern int total_output_atoms;

namespace internal
{

   //-------------------------------------------------------------------------
   // Internal data type definitions
   //-------------------------------------------------------------------------
   extern bool output_atoms_config;
   extern int output_atoms_config_rate;

   extern bool output_cells_config;
   extern int output_cells_config_rate;

   extern double field_output_min_1;
   extern double field_output_max_1;
   extern double field_output_min_2;
   extern double field_output_max_2;

   extern double atoms_output_min[3];
   extern double atoms_output_max[3];

   extern std::vector<int> local_output_atom_list;

   enum data_format
   {
      binary = 1,
      text = 2
   };
   extern data_format output_data_format;
   extern bool output_new;
   extern bool mpi_io;
   //-------------------------------------------------------------------------
   // Internal shared variables
   //-------------------------------------------------------------------------

   //-------------------------------------------------------------------------
   // Internal function declarations
   //-------------------------------------------------------------------------
   void atoms();
   void atoms_coords();
   void cells();
   void cells_coords();

   void atoms_new();
   void atoms_coords_new();
   void cells_new();
   void cells_coords_new();

   void write_data(const std::vector<float> &buffer, bool coord);
   void write_data_text(std::string filename, const std::vector<float> &buffer);
   void write_data_binary(std::string filename, const std::vector<float> &buffer);
   void copy_data_to_buffer(const std::vector<double> &x, // vector data
                           const std::vector<double> &y,
                           const std::vector<double> &z,
                           const std::vector<int> &mask,
                           std::vector<float> &buffer);
   void write_coordinate_meta();
   void write_meta(const double simulation_time, // time (seconds)
                   const double temperature,     // system temperature (Kelvin)
                   const double applied_field_x, // applied field components (Tesla)
                   const double applied_field_y,
                   const double applied_field_z,
                   const double magnetization_x, // magnetization components (normalized)
                   const double magnetization_y,
                   const double magnetization_z,
                   const int num_files);
   std::string data_filename(bool coords);

} // end of internal namespace

} // end of config namespace

#endif //CONFIG_INTERNAL_H_
