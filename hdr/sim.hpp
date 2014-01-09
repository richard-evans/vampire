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
#include <valarray>
#include <vector>

/// Enumerated lists for code readability
enum pump_functions_t {square=0, two_temperature, double_pump_two_temperature, double_pump_square};

namespace sim{
	extern std::ofstream mag_file;
	extern int time;
	extern int total_time;
	extern int loop_time;
	extern int partial_time;
	extern int equilibration_time;
	extern int runs;
	
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
	extern double Heq; //T
	extern double applied_field_angle_phi;
	extern double applied_field_angle_theta;
	extern bool applied_field_set_by_angle;
	
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
	extern double mc_delta_angle; /// Tuned angle for Monte Carlo trial move
	enum mc_algorithms { spin_flip, uniform, angle, hinzke_nowak};
   extern mc_algorithms mc_algorithm; /// Selected algorith for Monte Carlo simulations

	extern double head_position[2];
	extern double head_speed;
	extern bool   head_laser_on;
	
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
	extern double TTG;///electron coupling constant    
	extern double TTTe; /// electron temperature
	extern double TTTp; /// phonon temperature
	
	extern int system_simulation_flags;
	extern int hamiltonian_simulation_flags[10];
	
	extern int integrator;
	extern int program;
	extern int AnisotropyType;
	
	extern bool surface_anisotropy;
	extern bool identify_surface_atoms;
	extern unsigned int surface_anisotropy_threshold;
	extern bool NativeSurfaceAnisotropyThreshold;
	
	// Anisotropy control booleans
	extern bool UniaxialScalarAnisotropy; /// Enables scalar uniaxial anisotropy
	extern bool TensorAnisotropy; /// Overrides vector uniaxial anisotropy (even slower)
	extern bool second_order_uniaxial_anisotropy; /// Enables second order uniaxial anisotropy
   extern bool sixth_order_uniaxial_anisotropy; // Enables sixth order uniaxial anisotropy
	extern bool CubicScalarAnisotropy; // Enables scalar cubic anisotropy
	extern bool EnableUniaxialAnisotropyUnitVector; /// enables anisotropy tensor if any material has non z-axis K
   extern bool lattice_anisotropy_flag; /// Enables lattice anisotropy

	// Local system variables
	extern bool local_temperature; /// flag to enable material specific temperature
	extern bool local_applied_field; /// flag to enable material specific applied field
	extern bool local_fmr_field; /// flag to enable material specific fmr field

	
	// Wrapper Functions
	extern int run();
	extern int initialise();
	extern int integrate(int);
	
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
	extern int MonteCarlo();
	extern int ConstrainedMonteCarlo();
	extern int ConstrainedMonteCarloMonteCarlo();
	extern std::valarray<double> mc_move(std::valarray<double>&);

	// Integrator initialisers
	extern void CMCinit();
	extern int LLGinit();
	extern void CMCMCinit();

	// Field and energy functions
	extern double calculate_spin_energy(const int, const int);
   extern double spin_exchange_energy_isotropic(const int, const double, const double , const double );
   extern double spin_exchange_energy_vector(const int, const double, const double, const double);
   extern double spin_exchange_energy_tensor(const int, const double, const double, const double);
   extern double spin_scalar_anisotropy_energy(const int, const double);
   extern double spin_second_order_uniaxial_anisotropy_energy(const int, const double, const double, const double);
   extern double spin_sixth_order_uniaxial_anisotropy_energy(const int, const double, const double, const double);   
   extern double spin_lattice_anisotropy_energy(const int, const double, const double, const double);
   extern double spin_cubic_anisotropy_energy(const int, const double, const double, const double);
   extern double spin_tensor_anisotropy_energy(const int, const double, const double, const double);
   extern double spin_surface_anisotropy_energy(const int, const int, const double, const double, const double);
   extern double spin_applied_field_energy(const double, const double, const double);
   extern double spin_magnetostatic_energy(const int, const double, const double, const double);
   extern double lattice_anisotropy_function(const double, const int);

   // LaGrange multiplier variables
   extern double lagrange_lambda_x;
   extern double lagrange_lambda_y;
   extern double lagrange_lambda_z;
   extern double lagrange_m;
   extern double lagrange_N;
   extern bool   lagrange_multiplier;
   extern void   update_lagrange_lambda();
}

namespace cmc{
	
	class cmc_material_t {
	public:

		double constraint_phi; /// Constrained minimisation vector (azimuthal) [degrees]
		double constraint_phi_min; /// loop angle min [degrees]
		double constraint_phi_max; /// loop angle max [degrees]
		double constraint_phi_delta; /// loop angle delta [degrees]

		double constraint_theta; /// Constrained minimisation vector (rotational) [degrees]
		double constraint_theta_min; /// loop angle min [degrees]
		double constraint_theta_max; /// loop angle max [degrees]
		double constraint_theta_delta; /// loop angle delta [degrees]
		
		// performance optimised rotational matrices
		double ppolar_vector[3];
		double ppolar_matrix[3][3];
		double ppolar_matrix_tp[3][3];
		
		// vector magnetisation
		double M_other[3];
		
	cmc_material_t():
		constraint_phi(0.0),
		constraint_phi_min(0.0),
		constraint_phi_max(0.0),
		constraint_phi_delta(5.0),
		constraint_theta(0.0),
		constraint_theta_min(0.0),
		constraint_theta_max(0.0),
		constraint_theta_delta(5.0)
	
	{

	//for(int i=0;i<100;i++){
	//	geometry_coords[i][0]=0.0;
	//	geometry_coords[i][1]=0.0;
	//}	
}
	};
	
	extern std::vector<cmc_material_t> cmc_mat;
	
	extern bool is_initialised;
	
	extern int active_material; /// material in current hybrid loop
	
	extern std::vector<std::vector< int > > atom_list;
	extern double mc_success;
	extern double mc_total;
	extern double sphere_reject;
	extern double energy_reject;
}

#endif /*SIM_H_*/
