#ifndef LSF_H_
#define LSF_H_
/// Header file for LLG namespace
namespace LSF_arrays{

//==========================================================
// Namespace to store persistant LSF integration arrays
//==========================================================

   extern std::vector <double> x_lsf_array;
   extern std::vector <double> y_lsf_array;
   extern std::vector <double> z_lsf_array;

	extern std::vector <double> x_euler_array;
	extern std::vector <double> y_euler_array;
	extern std::vector <double> z_euler_array;

	extern std::vector <double> x_heun_array;
	extern std::vector <double> y_heun_array;
	extern std::vector <double> z_heun_array;

	extern std::vector <double> x_spin_storage_array;
	extern std::vector <double> y_spin_storage_array;
	extern std::vector <double> z_spin_storage_array;

	extern std::vector <double> x_initial_spin_array;
	extern std::vector <double> y_initial_spin_array;
	extern std::vector <double> z_initial_spin_array;

	extern bool LSF_set;

	extern std::vector <double> mod_S;

	extern std::vector <double> tx;
	extern std::vector <double> ty;
	extern std::vector <double> tz;

}
#endif /*LSF_H_*/
