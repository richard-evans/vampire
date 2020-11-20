//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Rory Pond and Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "config.hpp"
#include "vio.hpp"

// config module headers
#include "internal.hpp"

namespace config{

   namespace internal{

      //----------------------------------------------------------------------------
      // Function to initialize config module
      //----------------------------------------------------------------------------
      void initialize(){

         // Output informative message to log
         zlog << zTs() << "Initialising configuration output..." << std::flush;

         // Calculate total number of atoms for output
         #ifdef MPICF
            const uint64_t num_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
         #else
            const uint64_t num_atoms = atoms::num_atoms;
         #endif

         // get output bounds
         const double minB[3] = {atoms_output_min[0] * cs::system_dimensions[0],
                                 atoms_output_min[1] * cs::system_dimensions[1],
                                 atoms_output_min[2] * cs::system_dimensions[2]};

         const double maxB[3] = {atoms_output_max[0] * cs::system_dimensions[0],
                                 atoms_output_max[1] * cs::system_dimensions[1],
                                 atoms_output_max[2] * cs::system_dimensions[2]};

         // loop over all local atoms and determine atoms to be outputted
         for (uint64_t atom = 0; atom < num_atoms; atom++){

            const double cc[3] = {atoms::x_coord_array[atom], atoms::y_coord_array[atom], atoms::z_coord_array[atom]};

            // check atom within output bounds
            if ((cc[0] >= minB[0]) && (cc[0] <= maxB[0]))
            {
               if ((cc[1] >= minB[1]) && (cc[1] <= maxB[1]))
               {
                  if ((cc[2] >= minB[2]) && (cc[2] <= maxB[2]))
                  {
                     config::internal::local_output_atom_list.push_back(atom);
                  }
               }
            }

         }

         //------------------------------------------------------
         // calculate total atoms to output from all processors
         //------------------------------------------------------
         #ifdef MPICF
            uint64_t local_atoms = local_output_atom_list.size();
            uint64_t total_atoms;
            // add number of local output atoms on all processors
            MPI_Allreduce(&local_atoms, &total_atoms, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            config::internal::total_output_atoms = total_atoms;
         #else
            config::internal::total_output_atoms = local_output_atom_list.size();
         #endif

         // calculate total size of spin data in GB
         config::internal::io_data_size = 3.e-9 * double(sizeof(double)) * double(config::internal::total_output_atoms);

         // Resize local buffer
         config::internal::local_buffer.resize(3 * local_output_atom_list.size());

         // set number of io_processors to num processors for for file per process mode
         #ifdef MPICF
            if(config::internal::mode == fpprocess){
               config::internal::num_io_groups = vmpi::num_processors;
               config::internal::io_group_id = vmpi::my_rank;
            }
         #endif

         //------------------------------------------------------
         // For mpi-io output calculate data offsets
         //------------------------------------------------------
         #ifdef MPICF
            if(config::internal::mode == mpi_io){

               // array containing number of atoms per processor
               std::vector<uint64_t> atoms_per_processor(vmpi::num_processors, 0);

               // store number of local atoms in correct bin
               atoms_per_processor[vmpi::my_rank] = local_output_atom_list.size();

               // reduce on all processors
               MPI_Allreduce(MPI_IN_PLACE,&atoms_per_processor[0],vmpi::num_processors, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

               // calculate linear integer and 3 vector buffer offsets for my_rank
               uint64_t rank_offset = 0;
               for(unsigned int p=0; p < vmpi::my_rank; p++){
                  rank_offset += atoms_per_processor[p];
               }

               // set MPI offset (in bytes)
               config::internal::linear_offset = rank_offset * sizeof(int);
               config::internal::buffer_offset = rank_offset * 3 * sizeof(double);

            }
         #endif

         //------------------------------------------------------
         // For file per node output split MPI communicator
         //------------------------------------------------------
         #ifdef MPICF
            if(config::internal::mode == fpnode){

               // determine the group id to which I belong
               config::internal::io_group_id = vmpi::my_rank / ( 1 + (vmpi::num_processors - 1)/config::internal::num_io_groups);

               // Split communicator according to group id
               MPI_Comm_split(MPI_COMM_WORLD, config::internal::io_group_id, vmpi::my_rank, &config::internal::io_comm);

               // get my rank in group and group size
               MPI_Comm_rank(config::internal::io_comm, &config::internal::io_group_rank);
               MPI_Comm_size(config::internal::io_comm, &config::internal::io_group_size);

               // determine master process which actually outputs data from all processors in group
               config::internal::io_group_master_id  = config::internal::io_group_size - 1;

               // identify master process in io group
               if(config::internal::io_group_rank == config::internal::io_group_master_id) config::internal::io_group_master = true;

               // Reduce number of local atoms to be outputted on master process
               uint64_t local_size = local_output_atom_list.size();
               uint64_t collated_size = 0;

               MPI_Reduce(&local_size, &collated_size, 1, MPI_UINT64_T, MPI_SUM, config::internal::io_group_master_id, config::internal::io_comm);

               // resize output array on io master only
               if(config::internal::io_group_master) config::internal::collated_buffer.resize(3 * collated_size);

               // resize counts on master io process
               if(io_group_master) io_group_recv_counts.resize(io_group_size,0);
               if(io_group_master) io_group_displacements.resize(io_group_size,0);

               // get local data size
               int num_data = local_buffer.size();

               // gather number of data to be received from each processor in io group
               MPI_Gather(&num_data, 1, MPI_INT, &io_group_recv_counts[0], 1, MPI_INT, io_group_master_id, io_comm);

               // calculate displacements for gatherv and total number of points
               if(io_group_master){
                  int disp = 0; // initial displacement
                  for(int p=0; p<io_group_size; p++){
                     io_group_displacements[p] = disp;
                     disp += io_group_recv_counts[p]; // add number of counts to be recieved
                  }
               }

            }
         #endif

         // set initialised flag
         config::internal::initialised = true;

         // Output informative message to log
         zlog << " done!" << std::endl;

         return;

      }

   } // end of internal namespace

} // end of config namespace
