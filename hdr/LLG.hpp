#ifndef LLG_H_
#define LLG_H_
// Header file for LLG namespace
namespace LLG_arrays{
	
//==========================================================
// Namespace to store persistant LLG integration arrays
//==========================================================

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

	extern bool LLG_set;

}
#endif /*LLG_H_*/
