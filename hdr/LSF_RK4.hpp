#ifndef LSF_RK4_H_
#define LSF_RK4_H_

/// Header file for LSF-RK4 namespace
namespace LSF_RK4_arrays{

//==========================================================
// Namespace to store persistant LSF-RK4 integration arrays
//==========================================================

	extern std::vector<double> x_lsf_array;
	extern std::vector<double> y_lsf_array;
	extern std::vector<double> z_lsf_array;

	extern std::vector<double> x_initial_spin_array;
	extern std::vector<double> y_initial_spin_array;
	extern std::vector<double> z_initial_spin_array;

	extern std::vector<double> x_k1_array;
	extern std::vector<double> y_k1_array;
   extern std::vector<double> z_k1_array;

	extern std::vector<double> x_k2_array;
	extern std::vector<double> y_k2_array;
	extern std::vector<double> z_k2_array;

	extern std::vector<double> x_k3_array;
	extern std::vector<double> y_k3_array;
	extern std::vector<double> z_k3_array;

	extern std::vector<double> x_k4_array;
	extern std::vector<double> y_k4_array;
	extern std::vector<double> z_k4_array;

	extern bool LSF_RK4_set;

	extern std::vector<double> mod_S;

	extern std::vector<double> tx;
	extern std::vector<double> ty;
	extern std::vector<double> tz;

}
#endif /*LSF_RK4_H_*/