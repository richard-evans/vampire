//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
#ifndef SIM_H_
#define SIM_H_

//Headers
#include <cstdint>
#include <fstream>
#include <stdint.h>
#include <string>
#include <valarray>
#include <vector>

/// Enumerated lists for code readability
enum pump_functions_t {square=0, two_temperature, double_pump_two_temperature, double_pump_square};

namespace sim{

   extern std::ofstream mag_file;
   extern std::uint64_t time;
   extern std::uint64_t total_time;
   extern std::uint64_t loop_time;
   extern std::uint64_t partial_time;
   extern std::uint64_t equilibration_time;
   extern int runs;
   extern std::int64_t parity;
   extern std::uint64_t output_atoms_file_counter;
   extern std::uint64_t output_cells_file_counter;
   extern std::uint64_t output_rate_counter;

   extern bool ext_demag;

   extern double Tmax;
   extern double Tmin;
   extern double Teq;
   extern double temperature;
   extern double delta_temperature;
   extern double H_applied;
   extern double H_vec[3];
   extern double Hmin; // T
   extern double Hmax; // T
   extern double Hinc; // T
   extern double Heq; // T
   extern double applied_field_angle_phi;
   extern double applied_field_angle_theta;
   extern bool applied_field_set_by_angle;
   extern double fmr_field_strength; // Oscillating field strength (Tesla)
   extern double fmr_field_frequency; // Oscillating field frequency (GHz)
   extern std::vector<double> fmr_field_unit_vector; // Oscillating field direction
   extern double fmr_field; // Instantaneous value of the oscillating field strength H sin(wt)
   extern bool enable_fmr; // Flag to enable fmr field calculation
   extern std::int64_t iH;
   extern double H; // T

   extern double demag_factor[3];

   extern double constraint_phi; /// Constrained minimisation vector (azimuthal) [degrees]
   extern double constraint_theta; /// Constrained minimisation vector (rotational) [degrees]

   extern bool constraint_rotation; /// enables rotation of spins to new constraint direction
   extern bool constraint_phi_changed; /// flag to note change in phi
   extern bool constraint_theta_changed; /// flag to note change in theta

   extern double constraint_phi_min; /// loop angle min [degrees]
   extern double constraint_phi_max; /// loop angle max [degrees]
   extern double constraint_phi_delta; /// loop angle delta [degrees]

   extern double constraint_theta_min; /// loop angle min [degrees]
   extern double constraint_theta_max; /// loop angle max [degrees]
   extern double constraint_theta_delta; /// loop angle delta [degrees]

   // Monte Carlo variables
   extern int num_monte_carlo_preconditioning_steps;

   // extern double head_position[2];
   // extern double head_speed;
   // extern bool   head_laser_on;

   extern double cooling_time;
   extern int cooling_function_flag;
   extern pump_functions_t pump_function;
   extern double pump_time;
   extern double pump_power;
   extern double double_pump_time;
   extern double double_pump_power;
   extern double double_pump_Tmax;
   extern double double_pump_delay;
   extern double HeatSinkCouplingConstant;
   extern double TTCe; ///electron specific heat
   extern double TTCl; ///phonon specific heat
   extern double TTG;  ///electron coupling constant
   extern double TTTe; /// electron temperature
   extern double TTTp; /// phonon temperature

   extern bool applied_field_enabled; // include the applied field in the spin field (default true)
   extern bool thermal_field_enabled; // include the thermal field in the spin field (default true)

   // Local system variables
   extern bool local_temperature; /// flag to enable material specific temperature
   extern bool local_applied_field; /// flag to enable material specific applied field
   extern bool local_fmr_field; /// flag to enable material specific fmr field

   // Checkpoint flags and variables
   extern bool checkpoint_loaded_flag;  // Flag to determine if it is first step after loading checkpoint (true).
   extern bool load_checkpoint_flag; // Load spin configurations
   extern bool load_checkpoint_continue_flag; // Continue simulation from checkpoint time
   extern bool save_checkpoint_flag; // Save checkpoint
   extern bool save_checkpoint_continuous_flag; // save checkpoints during simulations
   extern int save_checkpoint_rate; // Default increment between checkpoints

   // Initialization functions
   extern void initialize(int num_materials);

   // User interface functions
   extern bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index);
   extern bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   // Wrapper Functions
   extern int run();
   extern int initialise();
   extern int integrate(std::uint64_t);

   // Field and energy functions
   extern double calculate_spin_energy(const int atom);
   extern double spin_applied_field_energy(const double, const double, const double);
   extern double spin_magnetostatic_energy(const int, const double, const double, const double);

   void calculate_spin_fields(const int start_index,const int end_index);
   void calculate_external_fields(const int start_index,const int end_index);
   void calculate_deterministic_external_fields(const int start_index,const int end_index); // external fields excluding thermal
   void calculate_thermal_fields_only(const int start_index,const int end_index); // thermal fields only (not HAMR)
   //spin temperature
    extern double compute_spin_temperature(const int start_index, // first atom for exchange interactions to be calculated
                 const int end_index,
                 const std::vector<int>& type_array, // type for atom
                 std::vector<double>& x_spin_array, // coord vectors for atoms
                 std::vector<double>& y_spin_array,
                 std::vector<double>& z_spin_array,
                 std::vector<double>& fields_array_x, //  vectors for fields
                 std::vector<double>& fields_array_y,
                 std::vector<double>& fields_array_z,
                 std::vector<double>& mu_s_array);

   extern double spin_temperature;

   //spin temperature
   extern double compute_spin_temperature(const int start_index, // first atom for exchange interactions to be calculated
                                          const int end_index,
                                          const std::vector<int>& type_array, // type for atom
                                          std::vector<double>& x_spin_array, // coord vectors for atoms
                                          std::vector<double>& y_spin_array,
                                          std::vector<double>& z_spin_array,
                                          std::vector<double>& fields_array_x, //  vectors for fields
                                          std::vector<double>& fields_array_y,
                                          std::vector<double>& fields_array_z,
                                          std::vector<double>& mu_s_array);

   extern double spin_temperature;

   // Monte Carlo statistics counters
   extern double mc_statistics_moves;
   extern double mc_statistics_reject;

   //------------------------------------------------------------------------
   // getter functions to give access to internal sim variables
   //------------------------------------------------------------------------
   std::vector<double> get_stt_polarization_unit_vector(); // unit vector spin polarization
   std::vector<double> get_stt_rj(); // array of stt relaxation constants
   std::vector<double> get_stt_pj(); // array of stt precession constants

}

/*namespace ckp{
        extern std::uint64_t parity;
        extern std::uint64_t output_atoms_file_counter;
        extern double H; // T
        extern std::uint64_t iH; // uT
} // end namespace ckp
*/

#endif /*SIM_H_*/
