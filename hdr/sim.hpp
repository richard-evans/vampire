//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
#ifndef SIM_H_
#define SIM_H_

//Headers
#include <fstream>
#include <stdint.h>
#include <string>
#include <valarray>
#include <vector>

/// Enumerated lists for code readability
enum pump_functions_t {square=0, two_temperature, double_pump_two_temperature, double_pump_square};

namespace sim{

	//track parameters
	extern double track_Ms;
	extern double track_bit_width;
	extern double track_bit_depth;
	extern double track_bit_size;

	extern double track_bit_gap;
	extern double track_track_gap;

	extern double cross_track_velocity;
	extern double down_track_velocity;

	extern int track_num_bits_per_track;
	extern int track_num_tracks;
	extern double LFA_scan_field_step;
	extern double Ms;
	extern double track_pos_x;
	extern double track_pos_z;
	extern bool LFA;

	// distance of tracks from read head
	extern double track_fly_height; // Angstroms

	extern double initial_down_track_position;
	extern double initial_cross_track_position;

	extern bool track_ms_file;

	extern std::vector < double > track_field_x;
	extern std::vector < double > track_field_y;
	extern std::vector < double > track_field_z;

	// enumerated list for integrators
	enum integrator_t{ llg_heun = 0, monte_carlo = 1, llg_midpoint = 2,
							 cmc = 3, hybrid_cmc = 4, llg_quantum = 5};

	extern std::ofstream mag_file;
	extern uint64_t time;
	extern uint64_t total_time;
	extern uint64_t loop_time;
	extern uint64_t partial_time;
	extern uint64_t equilibration_time;
	extern int runs;
	extern int64_t parity;
	extern uint64_t output_atoms_file_counter;
	extern uint64_t output_cells_file_counter;
	extern uint64_t output_rate_counter;

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
   extern int64_t iH;
   extern double H; // T

	extern double demag_factor[3];

	extern double constraint_phi; /// Constrained minimisation vector (azimuthal) [degrees]
	extern double constraint_theta; /// Constrained minimisation vector (rotational) [degrees]

	extern bool constraint_rotation; /// enables rotation of spins to new constraint direction
	extern bool constraint_phi_changed; /// flag to note change in phi
	extern bool constraint_theta_changed; /// flag to note change in theta

	extern double constraint_phi; /// Constrained minimisation vector (azimuthal) [degrees]
	extern double constraint_phi_min; /// loop angle min [degrees]
	extern double constraint_phi_max; /// loop angle max [degrees]
	extern double constraint_phi_delta; /// loop angle delta [degrees]

	extern double constraint_theta; /// Constrained minimisation vector (rotational) [degrees]
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

	extern int system_simulation_flags;
	extern int hamiltonian_simulation_flags[10];

	extern integrator_t integrator; // variable to specify integrator
	extern int program;

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
	extern int integrate(uint64_t);

	// Legacy integrators
	extern int LLB(int);
	extern int LLG(int);
	extern int LLG_relax(int);

	// New Integrators
	extern int LLG_Heun();
	extern int LLG_Heun_mpi();
	extern int LLG_Heun_cuda();
	extern int LLG_Midpoint();
	extern int LLG_Midpoint_mpi();
	extern int LLG_Midpoint_cuda();


	// Integrator initialisers
	extern int LLGinit();

	// Field and energy functions
	extern double calculate_spin_energy(const int atom);
   extern double spin_applied_field_energy(const double, const double, const double);
   extern double spin_magnetostatic_energy(const int, const double, const double, const double);

	void calculate_spin_fields(const int start_index,const int end_index);
	void calculate_external_fields(const int start_index,const int end_index);

   // LaGrange multiplier variables
   extern double lagrange_lambda_x;
   extern double lagrange_lambda_y;
   extern double lagrange_lambda_z;
   extern double lagrange_m;
   extern double lagrange_N;
   extern bool   lagrange_multiplier;
   extern void   update_lagrange_lambda();

	// Monte Carlo statistics counters
   extern double mc_statistics_moves;
   extern double mc_statistics_reject;

	extern int domain_wall_axis;
	extern double domain_wall_position;
	extern double domain_wall_discretisation;
	extern double domain_wall_centre;
	extern double domain_wall_width;
	extern std::vector < bool > anti_PBC;

	extern std::vector < double > domain_wall_second_vector_x;
	extern std::vector < double > domain_wall_second_vector_y;
	extern std::vector < double > domain_wall_second_vector_z;

	//------------------------------------------------------------------------
   // getter functions to give access to internal sim variables
   //------------------------------------------------------------------------
   std::vector<double> get_stt_polarization_unit_vector(); // unit vector spin polarization
   std::vector<double> get_stt_rj(); // array of stt relaxation constants
   std::vector<double> get_stt_pj(); // array of stt precession constants

}

/*namespace ckp{
        extern uint64_t parity;
        extern uint64_t output_atoms_file_counter;
        extern double H; // T
        extern uint64_t iH; // uT
} // end namespace ckp
*/

#endif /*SIM_H_*/
