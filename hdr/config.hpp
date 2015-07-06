//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------
//
//    Functions to output atomic coordinate and spin configurations to disk
//
//    Module enables bulk data output of atomic and macrocell level
//    coordinate and spin data to a series of configuration files. 
//    
//    Binary and text data formats are available to maximize peformance and 
//    portability respectively, with other specific data formats planned in
//    future. The data is designed to be read by the cfg2povray and cfg2rasmol
//    utilities included with the vampire distribution, as well as the vampire
//    spin viewer (vsv).
//
//    Data output is fully parallelized in MPI mode with a customizable number
//    of output processes per node, allowing for high scalability output for
//    parallel file systems.
//
//-----------------------------------------------------------------------------

// System headers
#include <string>
#include <vector>

// Program headers

#ifndef CONFIG_H_
#define CONFIG_H_

//--------------------------------------------------------------------------------
// Namespace for variables and functions to output configuration files
//--------------------------------------------------------------------------------
namespace config{

   //-----------------------------------------------------------------------------
   // Function to initialise configuration output
   //-----------------------------------------------------------------------------
   void initialize(const int num_atoms, // number of local atoms
                   const double system_dimensions_x, // system size
                   const double system_dimensions_y,
                   const double system_dimensions_z,
                   const std::vector<int>& material,
                   const std::vector<int>& category,
                   const std::vector<double>& spins_cx, // spin coordinates (Angstroms)
                   const std::vector<double>& spins_cy,
                   const std::vector<double>& spins_cz,
                   const std::vector<double>& cells_cx, // cell coordinates (Angstroms)
                   const std::vector<double>& cells_cy,
                   const std::vector<double>& cells_cz);

   //-----------------------------------------------------------------------------
   // Function to output spin and macrocell magnetic configurations
   //
   //
   //
   //
   //-----------------------------------------------------------------------------
   void output(const std::vector<double>& spins_x, // spin magnetization vectors (unit)
               const std::vector<double>& spins_y,
               const std::vector<double>& spins_z,
               const std::vector<double>& cells_x, // cell magnetization vectors (Bohr magnetons)
               const std::vector<double>& cells_y,
               const std::vector<double>& cells_z,
               const double simulation_time, // time (seconds)
               const double temperature, // system temperature (Kelvin)
               const double applied_field_x, // applied field components (Tesla)
               const double applied_field_y,
               const double applied_field_z,
               const double magnetization_x, // magnetization components (normalized)
               const double magnetization_y,
               const double magnetization_z);

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for config settings
   //-----------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

} // end of config namespace

#endif // CONFIG_H_
