#ifndef VMPI_H_
#define VMPI_H_
//==========================================================
// Includes
//==========================================================
#include "create.hpp"
#include <cstdlib>
#include <string>
#include <vector>

#ifdef MPICF
	#include <mpi.h>
#endif

using std::string;
//using std::vector;

namespace vmpi{
	
	extern int my_rank; 					///< Local CPU ID
	extern int num_processors;			///< Total number of CPUs
	extern int mpi_mode; 				///< MPI Simulation Mode (0 = Geometric Decomposition, 1 = Replicated Data, 2 = Statistical Parallelism)
	extern int num_core_atoms;			///< Number of atoms on local CPU with no external communication
	extern int num_bdry_atoms;			///< Number of atoms on local CPU with external communication
	extern int num_halo_atoms;			///< Number of atoms on remote CPUs needed for boundary atom integration
	
	extern char hostname[20];			///< Hostname of local CPU
	extern double start_time;			///< Simulation start time on local CPU
	extern double end_time;				///< Simulation end time on local CPU
	extern double min_dimensions[3]; 	///< Minimum coordinates of system on local cpu
	extern double max_dimensions[3]; 	///< Maximum coordinates of system on local cpu
	
	extern std::vector<int> send_atom_translation_array;
	extern std::vector<int> send_start_index_array;
	extern std::vector<int> send_num_array;
	extern std::vector<double> send_spin_data_array;

	extern std::vector<int> recv_atom_translation_array;
	extern std::vector<int> recv_start_index_array;
	extern std::vector<int> recv_num_array;
	extern std::vector<double> recv_spin_data_array;
	
	#ifdef MPICF
		extern std::vector<MPI::Request> requests;
		extern std::vector<MPI::Status> stati;
	#endif

	//functions declarations
	extern int initialise();
	extern int hosts();
	extern int finalise();
	extern int geometric_decomposition(int, double []);
	extern int crystal_xyz(std::vector<cs::catom_t> &);
	extern int copy_halo_atoms(std::vector<cs::catom_t> &);
	extern int set_replicated_data(std::vector<cs::catom_t> &);
	extern int identify_boundary_atoms(std::vector<cs::catom_t> &, std::vector<std::vector <int> > &);
	extern int init_mpi_comms(std::vector<cs::catom_t> & catom_array);

}

#endif /*VMPI_H_*/
