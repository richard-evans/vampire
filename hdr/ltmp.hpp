//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------
//
//   Functions to calculate dynamic lateral and vertical temperature
//   profiles with the two temperature model. Laser has gaussian absorption
//   profile in x,y and an exponential decrease in absorption with increasing
//   depth (z).
//
//                               |
//                               |
//                               |
//                               |
//                               |
//                              _-_
//                             /   \
//                         __--     --__
//                     ----------------------
//                     |  |  |  |  |  |  |  |
//                     ----------------------
//                     |  |  |  |  |  |  |  |
//                     ----------------------
//
//   The deposited laser energy is modified laterally by a gaussian heat profile
//   given by
//
//   P(r) = exp (-4 ln(2.0) (r**2) / (fwhm**2) )
//
//   and/or vertically considering an exponential depth dependence of the laser
//   energy given by
//
//   P(z) = exp(-z/penetration-depth)
//
//   Both of these can be combined to give a full 3D solution of the two
//   temperature model, including dynamic heat distribution within the sample.

// System headers
#include <string>
#include <vector>

// Program headers

#ifndef LOCALTEMPERATURE_H_
#define LOCALTEMPERATURE_H_

//--------------------------------------------------------------------------------
// Namespace for variables and functions to calculate localised temperature pulse
//--------------------------------------------------------------------------------
namespace ltmp{

   //--------------------------------------------------------------------
   // Class to contain parameters for tabulated absorption profile
   //
   // Tabulated values are read from a file and added point-wise to
   // the class. During initialisation interpolating functions
   // are determined to calculate A(z)
   //
   class abs_t{

      private:

         int z_max; // maximum array value in tabulated function
         double A_max; // value of absorption at z_max (used for all z > z_max)
         bool profile_loaded; // flag indicating profile has been loaded from file

         std::vector<int> z; // input z-height from surface values
         std::vector<double> A; // input absorption values
         std::vector<double> m; // calculated m value
         std::vector<double> c; // calculated c value

      public:
         abs_t();
         bool is_set();
         void add_point(double height, double absorption);
         void set_interpolation_table();
         double get_absorption_constant(double height);
         void output_interpolated_function(int height);

   };

   //-----------------------------------------------------------------------------
   // Variables used for the localised temperature pulse calculation
   //-----------------------------------------------------------------------------
   extern abs_t absorption_profile; // class variable containing tabulated absorption profile

   //-----------------------------------------------------------------------------
   // Function to check local temperature pulse is enabled
   //-----------------------------------------------------------------------------
   bool is_enabled();

   //-----------------------------------------------------------------------------
   // Function to initialise localised temperature pulse calculation
   //-----------------------------------------------------------------------------
   void initialise(const double system_dimensions_x,
                  const double system_dimensions_y,
                  const double system_dimensions_z,
                  const std::vector<double>& atom_coords_x,
                  const std::vector<double>& atom_coords_y,
                  const std::vector<double>& atom_coords_z,
                  const std::vector<int>& atom_type_array,
                  const int num_local_atoms,
                  const double starting_temperature,
                  const double pump_power,
                  const double pump_time,
                  const double TTG,
                  const double TTCe,
                  const double TTCl,
                  const double dt,
                  const double Tmin,
                  const double Tmax);

   //-----------------------------------------------------------------------------
   // Function to copy localised thermal fields to external field array
   //-----------------------------------------------------------------------------
   void get_localised_thermal_fields(std::vector<double>& x_total_external_field_array,
                               std::vector<double>& y_total_external_field_array,
                               std::vector<double>& z_total_external_field_array,
                               const int start_index,
                               const int end_index);

   //-----------------------------------------------------------------------------
   // Function for updating localised temperature
   //-----------------------------------------------------------------------------
   void update_localised_temperature(const double start_from_start);

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for ltmp settings
   //-----------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

} // end of ltmp namespace

#endif // LOCALTEMPERATURE_H_
