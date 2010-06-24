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
	extern double Tmax;
	extern double temperature;
	extern double H_applied;
	extern double H_vec[3];
	extern double Hmin; // T
	extern double Hmax; // T
	extern double Hinc; // T
	
	extern int system_simulation_flags;
	extern int hamiltonian_simulation_flags[10];
	
	extern int run();
	extern int initialise();
	
	extern int LLB(int);
	extern int LLG(int);
	extern int LLG_relax(int);
}

#endif /*SIM_H_*/
