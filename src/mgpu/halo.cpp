//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) S Jenkins 2026. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

// Vampire headers
#include "errors.hpp"
#include "mgpu.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// mgpu module headers
#include "internal.hpp"

namespace mgpu{

namespace internal{

   int num_send = 0;
   int num_recv = 0;
   void* h_send_buffer = nullptr;
   void* h_recv_buffer = nullptr;
   int real_size = 8;

   uint64_t num_swaps = 0; // number of halo updates
   double mpi_time = 0.0;  // wall time in seconds spent in MPI, including waiting for other ranks

   //---------------------------------------------------------------------------
   // Function to exchange the packed halo spins with the other ranks (the spins
   // a rank sends to itself are copied on the device instead)
   //---------------------------------------------------------------------------
   void exchange_halo(){

      const double start = MPI_Wtime();
      const MPI_Datatype type = (real_size == sizeof(double)) ? MPI_DOUBLE : MPI_FLOAT;
      char* send = static_cast<char*>(h_send_buffer);
      char* recv = static_cast<char*>(h_recv_buffer);

      std::vector<MPI_Request> requests;
      for(int p = 0; p < vmpi::num_processors; p++){
         if(p == vmpi::my_rank) continue; // copied on the device (start_halo)
         if(peer_send_num[p] != 0){
            requests.push_back(MPI_REQUEST_NULL);
            MPI_Isend(send + 3 * peer_send_start[p] * real_size, 3 * peer_send_num[p], type, p, 48, MPI_COMM_WORLD, &requests.back());
         }
         if(peer_recv_num[p] != 0){
            requests.push_back(MPI_REQUEST_NULL);
            MPI_Irecv(recv + 3 * peer_recv_start[p] * real_size, 3 * peer_recv_num[p], type, p, 48, MPI_COMM_WORLD, &requests.back());
         }
      }
      MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
      mpi_time += MPI_Wtime() - start;
      num_swaps++;

   }

   //---------------------------------------------------------------------------
   // Function to update the spins of the halo atoms on the device
   //---------------------------------------------------------------------------
   void swap_halo(){
      if(num_send == 0 && num_recv == 0) return;
      start_halo();
      finish_halo();
   }

   //---------------------------------------------------------------------------
   // Function to write the atom counts and halo exchange times of every rank.
   //---------------------------------------------------------------------------
   void report_halo_times(){

      const int n = 8;
      const double mine[n] = {double(vmpi::num_core_atoms + vmpi::num_bdry_atoms), double(num_recv), double(num_send),
                              double(num_swaps), halo_copy_time(), mpi_time, halo_wait_time(), double(device_only_halo)};
      std::vector<double> all(n * vmpi::num_processors);
      MPI_Gather(mine, n, MPI_DOUBLE, all.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if(vmpi::my_rank != 0) return;

      std::ostringstream report;
      report << std::fixed << std::setprecision(3);
      for(int p = 0; p < vmpi::num_processors; p++){
         const double* r = &all[n * p];
         report << "Multi-GPU rank " << p << ": " << uint64_t(r[0]) << " local atoms, " << uint64_t(r[1]) << " halo atoms received, "
                << uint64_t(r[2]) << " sent; " << uint64_t(r[3]) << " halo updates: ";
         if(r[7] != 0.0) report << "copied on the device, not timed\n";
         else report << "GPU copies " << r[4] << " s, MPI " << r[5] << " s, GPU waiting " << r[6] << " s\n";
      }
      std::cout << report.str() << std::flush;
      zlog << zTs() << "Halo exchange times\n" << report.str() << std::flush;

   }

} // end of internal namespace

//------------------------------------------------------------------------------
// Function to set up halo exchange on the device and, if requested, to order
// the exchange neighbours of each atom by atom key
//------------------------------------------------------------------------------
void initialize_device(){

   if(err::check==true) std::cout << "mgpu::initialize_device has been called" << std::endl;

   if(decomposition_independent_order) internal::sort_neighbours_by_key();
   internal::initialize_noise();

   internal::num_send = vmpi::send_atom_translation_array.size();
   internal::num_recv = vmpi::recv_atom_translation_array.size();
   internal::allocate_device_halo(vmpi::send_atom_translation_array, vmpi::recv_atom_translation_array);

   zlog << zTs() << "Multi-GPU halo exchange: " << internal::num_send << " atoms sent and "
        << internal::num_recv << " atoms received per update" << std::endl;

}

//------------------------------------------------------------------------------
// Function to free buffers of the module
//------------------------------------------------------------------------------
void finalize(){
   if(err::check==true) std::cout << "mgpu::finalize has been called" << std::endl;
   internal::report_halo_times();
   internal::free_device_halo();
   internal::free_noise();
}

} // end of mgpu namespace
