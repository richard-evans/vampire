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

	extern double constraint_phi; // Constrained minimisation vector (azimuthal) [degrees]
	extern double constraint_phi_min; // loop angle min [degrees]
	extern double constraint_phi_max; // loop angle max [degrees]
	extern double constraint_phi_delta; // loop angle delta [degrees]

	extern double constraint_theta; // Constrained minimisation vector (rotational) [degrees]
	extern double constraint_theta_min; // loop angle min [degrees]
	extern double constraint_theta_max; // loop angle max [degrees]
	extern double constraint_theta_delta; // loop angle delta [degrees]
	
	extern double head_position[2];
	extern double head_speed;
	extern bool   head_laser_on;
	
	extern double cooling_time;
	extern int cooling_function_flag;
	extern double pump_time;
	extern double pump_power;
        extern double HeatSinkCouplingConstant;
        extern double TTCe; //electron specific heat
        extern double TTCl; //phonon specific heat
        extern double TTG;//electron coupling constant    

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
	extern bool UniaxialScalarAnisotropy; // Enables scalar uniaxial anisotropy
	extern bool TensorAnisotropy; // Overrides vector uniaxial anisotropy (even slower)
	extern bool CubicScalarAnisotropy; // Enables scalar cubic anisotropy
	extern bool EnableUniaxialAnisotropyUnitVector; // enables anisotropy tensor if any material has non z-axis K

	
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
	
	// Integrator initialisers
	extern void CMCinit();
	extern int LLGinit();
	extern void CMCMCinit();
	
	// Field and energy functions
	extern double calculate_spin_energy(const int, const int);

}

namespace cmc{
	
	class cmc_material_t {
	public:

		double constraint_phi; // Constrained minimisation vector (azimuthal) [degrees]
		double constraint_phi_min; // loop angle min [degrees]
		double constraint_phi_max; // loop angle max [degrees]
		double constraint_phi_delta; // loop angle delta [degrees]

		double constraint_theta; // Constrained minimisation vector (rotational) [degrees]
		double constraint_theta_min; // loop angle min [degrees]
		double constraint_theta_max; // loop angle max [degrees]
		double constraint_theta_delta; // loop angle delta [degrees]
		
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
	
	extern std::vector<std::vector< int > > atom_list;
	extern double mc_success;
	extern double mc_total;
	extern double sphere_reject;
	extern double energy_reject;
}

#endif /*SIM_H_*/
