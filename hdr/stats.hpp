#include <vector>

namespace stats
//==========================================================
// Namespace statistics
//==========================================================
{
	extern int num_atoms;				// Number of atoms for statistic purposes
	extern double inv_num_atoms;	//1.0/num_atoms
	extern double max_moment;		// Total Maximum moment

	extern double total_mag_actual[3];	///< Actual magnetisation components
	extern double total_mag_m_actual;	///< Actual magnitude of total magnetisation
	extern double total_mean_mag_m_actual;	///< Actual magnitude of total magnetisation
	
	extern double total_mag_norm[3];	///< Normalised magnetisation components
	extern double total_mag_m_norm;	///< Normalised magnitude of total magnetisation
	extern double total_mean_mag_m_norm;	///< Normalised magnitude of total magnetisation

	extern double data_counter;		// number of data points for averaging

	// Member Functions
	extern int mag_m();
	extern void mag_m_reset();
	extern double max_torque();
	
	extern std::vector <double> sublattice_mx_array;
	extern std::vector <double> sublattice_my_array;
	extern std::vector <double> sublattice_mz_array;
	extern std::vector <double> sublattice_magm_array;
	extern std::vector <double> sublattice_mean_magm_array;
	extern std::vector <double> sublattice_mom_array;
	extern std::vector <int> sublattice_nm_array;
	
	extern bool calculate_torque;
	extern double total_system_torque[3];
	extern double total_mean_system_torque[3];
	
	extern std::vector <double> sublattice_torque_x_array;
	extern std::vector <double> sublattice_torque_y_array;
	extern std::vector <double> sublattice_torque_z_array;

	extern double torque_data_counter;
}
