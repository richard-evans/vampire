//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
#ifndef VMPI_H_
#define VMPI_H_
//==========================================================
// Includes
//==========================================================
#include "create.hpp"
#include <cstdlib>
#include <string>
#include <cstdint>
#include <vector>

#ifdef MPICF
   #include <mpi.h>
#endif

namespace vmpi{

   extern bool master;              // boolean variable (only true on master process)
   extern int master_id;            // MPI process ID for master process
   extern int my_rank;              ///< Local CPU ID
   extern int num_processors;       ///< Total number of CPUs
   extern int mpi_mode;             ///< MPI Simulation Mode (0 = Geometric Decomposition, 1 = Replicated Data, 2 = Statistical Parallelism)
   extern unsigned int ppn;         ///< Processors per node
   extern int num_core_atoms;       ///< Number of atoms on local CPU with no external communication
   extern int num_bdry_atoms;       ///< Number of atoms on local CPU with external communication
   extern int num_halo_atoms;       ///< Number of atoms on remote CPUs needed for boundary atom integration
   extern int num_local_atoms;       // number of atoms on local processor

   extern int num_io_processors;    ///< Total number of CPUs that perform IO
   extern int size_io_group;        ///< Size of io mpi groups
   extern int my_io_rank;           ///< Local CPU IO Comm Group Rank
   extern int my_io_group;          ///< Local CPU IO Comm Group Rank
   extern int io_processor;         ///< The group rank of processor who performs IO
#ifdef MPICF
   extern MPI_Comm io_comm;         ///< MPI Communicator for IO
#endif


   extern bool replicated_data_staged; ///< Flag for staged system generation

   extern std::string hostname;        ///< Hostname of local CPU
   extern double min_dimensions[3];    ///< Minimum coordinates of system on local cpu
   extern double max_dimensions[3];    ///< Maximum coordinates of system on local cpu

   // Timing variables
   extern double start_time;        ///< Simulation start time on local CPU
   extern double end_time;          ///< Simulation end time on local CPU
   extern double ComputeTime;       /// Temporary for storing time
   extern double WaitTime;          /// Temporary for storing time
   extern double TotalComputeTime;  /// Total time spent in computation
   extern double TotalWaitTime;     /// Total time spent waiting
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
   extern std::vector<double> send_coord_data_array;

   extern std::vector<int> recv_atom_translation_array;
   extern std::vector<int> recv_start_index_array;
   extern std::vector<int> recv_num_array;
   extern std::vector<double> recv_spin_data_array;
   extern std::vector<double> recv_coord_data_array;

   #ifdef MPICF
      extern std::vector<MPI_Request> requests;
      extern std::vector<MPI_Status> stati;
   #endif

   //functions declarations
   extern void initialise(int argc, char *argv[]);
   extern int hosts();
   extern int finalise();
   extern void geometric_decomposition(int, double [3]);
   extern double SwapTimer(double, double&);

   // functions for sending/receiving halo data
   extern void mpi_init_halo_swap();
   extern void mpi_complete_halo_swap();

   // functions for sending/receiving halo data
   extern void mpi_init_halo_swap_coords();
   extern void mpi_complete_halo_swap_coords();

   // wrapper functions avoiding MPI library
   extern void barrier();
   extern std::uint64_t reduce_sum(std::uint64_t local);
   extern double reduce_sum(double local);
   extern double reduce_max(double local);
   extern std::uint64_t all_reduce_sum(std::uint64_t local);
   extern double all_reduce_sum(double local);
   extern void all_reduce_sum(std::vector<double>& array);
   extern void all_reduce_sum(std::vector<int>& array);

   extern void collate(std::vector<double>& input, std::vector<double>& output);
   extern void counts_and_displacements(std::vector<double>& input, std::vector<double>& output, std::vector<int>& counts, std::vector<int>& displacements);
   extern void fast_collate(std::vector<double>& input, std::vector<double>& output, std::vector<int>& counts, std::vector<int>& displacements);
   extern void broadcast(std::vector<char>& message, int source_rank);

   // function to seed random numbers in parallel
   uint32_t parallel_rng_seed(int seed);

   //---------------------------------------------------------------------------
   // Function to process input file parameters for mpi module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

}

#endif /*VMPI_H_*/
