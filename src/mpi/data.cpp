//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2018. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "vmpi.hpp"

// Internal vmpi header

namespace vmpi{

   bool master = false; // boolean variable (only true on master process)
   int master_id = 0; // MPI process ID for master process

   int mpi_mode=0;
   unsigned int ppn=1;  ///< Processors per node
   int my_rank=0;
   int num_processors=1;
   int num_core_atoms;
   int num_bdry_atoms;
   int num_halo_atoms;
   int num_local_atoms; // number of local atoms on processor

   bool replicated_data_staged=false;

   std::string hostname;

   // timing variables
   double start_time;
   double end_time;
   double ComputeTime;			/// Temporary for storing time
   double WaitTime;				/// Temporary for storing time
   double TotalComputeTime;	/// Total time spent in computation
   double TotalWaitTime;		/// Total time spent waiting
   double AverageComputeTime;
   double AverageWaitTime;
   double MaximumComputeTime;
   double MaximumWaitTime;
   bool DetailedMPITiming=false;
   std::vector<double> ComputeTimeArray(0);
   std::vector<double> WaitTimeArray(0);

   double min_dimensions[3]; ///< Minimum coordinates of system on local cpu
   double max_dimensions[3]; ///< Maximum coordinates of system on local cpu

   std::vector<int> send_atom_translation_array;
   std::vector<int> send_start_index_array;
   std::vector<int> send_num_array;
   std::vector<double> send_spin_data_array;

   std::vector<int> recv_atom_translation_array;
   std::vector<int> recv_start_index_array;
   std::vector<int> recv_num_array;
   std::vector<double> recv_spin_data_array;
   #ifdef MPICF
   std::vector<MPI_Request> requests(0);
   std::vector<MPI_Status> stati(0);
   #endif

}
