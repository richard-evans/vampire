//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------
Rory - this works and is nearly finished
// C++ standard library headers
#include <iomanip>
#include <fstream>
#include <sstream>

// Vampire headers
#include "config.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

// config headers
#include "internal.hpp"

namespace config{
   namespace internal{
      namespace mpi{

         int my_io_rank = 0; // my id in i/o communicator
         int my_io_master = 0; // my master id in i/o communicator who I send data to
         int num_io_processes_in_my_group = 1;

         // declare new split communicator for interprocess i/o
         #ifdef MPICF
            MPI_Comm MPI_COMM_IO;
         #endif

         std::vector<float> storage_buffer(0); // temporary buffer to store output in parallel version
         std::vector<int> buffer_offset(0); // offsets for mpi buffer locations
         std::vector<int> buffer_sizes(0); // number of data points from each process in MPI_COMM_IO

         //---------------------------------------------------------------------
         // Function to initialise MPI_COMM_IO communicator and data structures
         //---------------------------------------------------------------------
         void initialize(){

            #ifdef MPICF

            // check for post initialization of num_io_processes
            if(config::internal::set_num_io_nodes_to_ppn) config::internal::num_io_nodes = vmpi::ppn;

            // if using parallel i/o print informative message to log file
            if(config::internal::num_io_nodes > 1) zlog << zTs() << "Using parallel configuration output with " << config::internal::num_io_nodes << " output processes" << std::endl;

            // num_processes_per_group (num_processors/num_io_nodes) rounding up
            const int num_processes_per_group = (vmpi::num_processors - 1)/config::internal::num_io_nodes + 1;

            // determine my io group id
            const int my_group_id = vmpi::my_rank/num_processes_per_group;

            // split MPI_COMM_WORLD to form MPI_COMM_IO with same group ID
            MPI_Comm_split(MPI_COMM_WORLD, my_group_id, vmpi::my_rank, &MPI_COMM_IO);

            // get size of group and rank in group
            MPI_Comm_size(MPI_COMM_IO, &num_io_processes_in_my_group );
            MPI_Comm_rank(MPI_COMM_IO, &my_io_rank);

            // set last process in MPI_COMM_IO to output data
            my_io_master = num_io_processes_in_my_group - 1;

            // set io flag on master processes
            if(config::internal::mpi::my_io_rank == my_io_master){
               config::internal::io_master = true;
               // set buffer to store number of atoms/spins from each rank in i/o group
               buffer_sizes.resize(num_io_processes_in_my_group);
            }

            // determine number of atoms from each process
            int my_num_output_atoms = config::internal::local_output_atom_list.size();

            // gather number of atoms on master io process
            MPI_Gather(&my_num_output_atoms, 1, MPI_INT, &buffer_sizes[0], 1, MPI_INT, config::internal::mpi::my_io_master, MPI_COMM_IO);

            // Calculate total number of atoms to output
            int total_size = 0;
            if(my_io_rank == my_io_master){
               for(int i=0; i < num_io_processes_in_my_group; i++) total_size+=buffer_sizes[i];
               for(int i=0; i < num_io_processes_in_my_group; i++) std::cerr << "io\t" << vmpi::my_rank << "\t" << i << "\t" << buffer_sizes[i] << std::endl;
               std::cerr << total_size << std::endl;
            }

            std::cerr << vmpi::my_rank << "\t" << my_io_rank << "\t" << my_group_id << std::endl;

            MPI_Barrier(MPI_COMM_WORLD);

            err::vexit();

            config::internal::total_output_atoms = config::internal::local_output_atom_list.size();

            // Rory - got to here. Need to define output buffers here and collate data on master process

            // determine offset from each process


            // resize spin buffer on master process
            //if(config::internal::io_master){

               // calculate total number of output atoms
               //config::internal::total_output_atoms


            /* counts and displacement arrays are only required on the root */
   /*if(myid == mpi_root){
      counts=(int*)malloc(numnodes*sizeof(int));
      displacements=(int*)malloc(numnodes*sizeof(int));
   }*/
/* we gather the counts to the root */
   /*mpi_err = MPI_Gather((void*)myray,1,MPI_INT,
                    (void*)counts,  1,MPI_INT,
                    mpi_root,MPI_COMM_WORLD);*/
/* calculate displacements and the size of the recv array */
   /*if(myid == mpi_root){
      displacements[0]=0;
      for( i=1;i<numnodes;i++){
         displacements[i]=counts[i-1]+displacements[i-1];
      }
      size=0;
      for(i=0;i< numnodes;i++)
         size=size+counts[i];
      allray=(int*)malloc(size*sizeof(int));
   }*/
/* different amounts of data from each processor  */
/* is gathered to the root */
/*   mpi_err = MPI_Gatherv(myray, mysize,         MPI_INT,
                    allray,counts,displacements,MPI_INT,
                    mpi_root,
                    MPI_COMM_WORLD);*/


//determine if I am master process


  // set local buffers for all nodes

   //resize main buffer to zero on non-io nodes, to 3*tot spins on master node

            // Set CPUID on non-root process
         //   filename << std::setfill('0') << std::setw(5) << config vmpi::my_rank << "-";
// aggregate bandwidth
   #else
   // serial
      config::internal::total_output_atoms = config::internal::local_output_atom_list.size();
   #endif

}

//config::mpi::synchronize_buffer(){

  // copy data to local buffer

  // mpi reduce to buffer

  // swap buffers on master process


      } // end of namespace mpi
   } // end of internal namespace
} // end of config namespace
