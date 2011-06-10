#ifndef SIM_H_
#define SIM_H_

//Headers
#include <fstream>

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
	extern double demag_factor[3];
	
	extern double constraint_phi; // Constrained minimisation vector (azimuthal) [degrees]
	extern double constraint_theta; // Constrained minimisation vector (rotational) [degrees]

	extern double head_position[2];
	extern double head_speed;
	extern bool   head_laser_on;
	
	extern double cooling_time;
	extern int cooling_function_flag;
	extern double pump_time;
	extern double pump_power;
	
	extern int system_simulation_flags;
	extern int hamiltonian_simulation_flags[10];
	
	extern int integrator;
	extern int program;
	
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
	
	// Integrator initialisers
	extern void CMCinit();
	extern int LLGinit();
	
	// Field and energy functions
	extern double calculate_spin_energy(const int);

}

#endif /*SIM_H_*/
