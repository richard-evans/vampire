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
   extern unsigned int ppn;			///< Processors per node
	extern int num_core_atoms;			///< Number of atoms on local CPU with no external communication
	extern int num_bdry_atoms;			///< Number of atoms on local CPU with external communication
	extern int num_halo_atoms;			///< Number of atoms on remote CPUs needed for boundary atom integration

	extern int num_io_processors;		///< Total number of CPUs that perform IO
	extern int size_io_group;			///< Size of io mpi groups
	extern int my_io_rank;				///< Local CPU IO Comm Group Rank
	extern int my_io_group;				///< Local CPU IO Comm Group Rank
	extern int io_processor;			///< The group rank of processor who performs IO
#ifdef MPICF
	extern MPI_Comm io_comm;			///< MPI Communicator for IO
#endif


	extern bool replicated_data_staged; ///< Flag for staged system generation

	extern char hostname[20];			///< Hostname of local CPU
	extern double min_dimensions[3]; 	///< Minimum coordinates of system on local cpu
	extern double max_dimensions[3]; 	///< Maximum coordinates of system on local cpu

	// Timing variables
	extern double start_time;			///< Simulation start time on local CPU
	extern double end_time;				///< Simulation end time on local CPU
	extern double ComputeTime;			/// Temporary for storing time
	extern double WaitTime;				/// Temporary for storing time
	extern double TotalComputeTime;	/// Total time spent in computation
	extern double TotalWaitTime;		/// Total time spent waiting
	extern double AverageComputeTime;
	extern double AverageWaitTime;
	extern double MaximumComputeTime;
	extern double MaximumWaitTime;
	extern std::vector<double> ComputeTimeArray;
	extern std::vector<double> WaitTimeArray;
	extern bool DetailedMPITiming; /// flag to control logging of compute and wait times

	extern std::vector<int> send_atom_translation_array;
	extern std::vector<int> send_start_index_array;
	extern std::vector<int> send_num_array;
	extern std::vector<double> send_spin_data_array;

	extern std::vector<int> recv_atom_translation_array;
	extern std::vector<int> recv_start_index_array;
	extern std::vector<int> recv_num_array;
	extern std::vector<double> recv_spin_data_array;

	#ifdef MPICF
		extern std::vector<MPI_Request> requests;
		extern std::vector<MPI_Status> stati;
	#endif

	//functions declarations
	extern int initialise(int argc, char *argv[]);
	extern int hosts();
	extern int finalise();
   extern void geometric_decomposition(int, double []);
	extern int crystal_xyz(std::vector<cs::catom_t> &);
	extern int copy_halo_atoms(std::vector<cs::catom_t> &);
	extern int set_replicated_data(std::vector<cs::catom_t> &);
	extern int identify_boundary_atoms(std::vector<cs::catom_t> &, std::vector<std::vector <cs::neighbour_t> > &);
	extern int init_mpi_comms(std::vector<cs::catom_t> & catom_array);
	extern double SwapTimer(double, double&);

	// wrapper functions avoiding MPI library
	extern void barrier();

}

#endif /*VMPI_H_*/
