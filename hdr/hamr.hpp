//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) Andrea Meo 2022. All rights reserved.
//
//------------------------------------------------------------------------------

// Program headers
#ifndef HAMR_H_
#define HAMR_H_

// System headers
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace hamr{

   //-----------------------------------------------------------------------------
   // Function to initialise HAMR calculation
   //-----------------------------------------------------------------------------
   void initialize(const double Hmin,
                   const double Hmax,
                   const double Tmin,
                   const double Tmax,
                   const double system_dimensions_x,
                   const double system_dimensions_y,
                   const double system_dimensions_z,
                   const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z,
                   const std::vector<int>& atom_type_array,
                   const int num_local_atoms);

   //-----------------------------------------------------------------------------
   // Function to check local temperature pulse is enabled
   //-----------------------------------------------------------------------------
   void fields(const int start_index,
               const int end_index,
               double H_applied,
               const double temperature,
               const double Hvecx,
               const double Hvecy,
               const double Hvecz,
               std::vector<double>& x_total_external_field_array,
               std::vector<double>& y_total_external_field_array,
               std::vector<double>& z_total_external_field_array);

   //-----------------------------------------------------------------------------
   // Function to run hamr continuous simulation
   //-----------------------------------------------------------------------------
   void hamr_continuous();

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for hamr settings
   //-----------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //-----------------------------------------------------------------------------
   // Function to HAMR parameters
   //-----------------------------------------------------------------------------
   double get_head_position_x();
   double get_head_position_y();
   double get_head_speed();
   double get_FWHM_x();
   double get_FWHM_y();
   double get_laser_sigma_x();
   double get_laser_sigma_y();
   double get_field_bounds_x();
   double get_field_bounds_y();
   double get_bit_size();
   double get_track_size();
   double get_field_rise_time();
   double get_field_fall_time();
   double get_NPS();
   double get_track_padding();
   bool get_initialisation_state();

   //-----------------------------------------------------------------------------
   // Functions to write and read the state of the write head in checkpoint files,
   // so that a HAMR simulation can be continued from a checkpoint
   //-----------------------------------------------------------------------------
   void save_checkpoint(std::ostream& chkfile);
   void load_checkpoint(std::istream& chkfile);

   //-----------------------------------------------------------------------------
   // Active window (hamr:active-window)
   //-----------------------------------------------------------------------------
   bool active_window_requested(); // hamr:active-window set for the HAMR program (for the decomposition)
   bool active_window_enabled(); // true once the write has started with the window
   const std::vector<int>& update_active_window(); // update for the current head position; returns the [first, last+1) pairs of atoms to integrate
   bool atom_in_window(const int atom); // the atom (local or halo) belongs to a grain of the window
   uint64_t window_changes();           // number of times the grains of the window have changed (same on every process)

   //-----------------------------------------------------------------------------
   // Spin snapshots with the active window: whether the next snapshot may hold
   // only the atoms that changed since the previous one, whether an atom has
   // changed, and notice that a snapshot has been written
   //-----------------------------------------------------------------------------
   bool window_snapshot_possible();
   bool spin_changed_since_snapshot(const int atom);
   void snapshot_written();

   //-----------------------------------------------------------------------------------------------
   // Externally visible variables used for HAMR calculation
   //-----------------------------------------------------------------------------------------------
   extern bool head_laser_on;


} // end of hamr namespace

#endif // HAMR_H_
