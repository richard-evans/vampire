//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) S Jenkins 2026. All rights reserved.
//
//------------------------------------------------------------------------------
#ifndef MGPU_INTERNAL_H_
#define MGPU_INTERNAL_H_

// C++ standard library headers
#include <string>
#include <vector>

namespace mgpu{

   namespace internal{

      //------------------------------------------------------------------------
      // Halo exchange buffers: pinned host buffers of x,y,z triplets in the
      // precision of the CUDA build, holding the selected atoms of each rank in turn
      //------------------------------------------------------------------------
      extern int num_send;        // number of atoms sent
      extern int num_recv;        // number of atoms received
      extern void* h_send_buffer; // pinned host buffers (3 values per atom)
      extern void* h_recv_buffer;
      extern int real_size;       // size of one value in bytes (8 or 4)
      // first atom and number of atoms of each rank in the send and receive buffers
      extern std::vector<int> peer_send_start, peer_send_num, peer_recv_start, peer_recv_num;

      // device side (kernels.cu)
      int device_count();
      void device_pci_id(int device, char* id, int size);
      void allocate_device_halo(const std::vector<int>& send_atoms, const std::vector<int>& recv_atoms);
      void free_device_halo();
      // The halo update runs on its own stream: start_halo() packs the spins into
      // h_send_buffer, finish_halo() exchanges the buffers and returns the received spins
      // to the device; work launched in between overlaps. device_only_halo: all halo spins come from this rank.
      extern bool device_only_halo;
      void start_halo();
      void finish_halo();
      void select_halo_atoms(); // chooses the halo atoms to copy when the HAMR active window changes
      double halo_copy_time(); // GPU time in seconds spent packing, unpacking and copying halo spins
      double halo_wait_time(); // GPU time in seconds the default stream waited for halo spins

      // host side (halo.cpp): exchanges h_send_buffer and h_recv_buffer with the other ranks
      void exchange_halo();

      // full halo update of the spins on the device without overlap
      void swap_halo();

      // halo update and fields of the local atoms; external = true adds the
      // external (thermal, HAMR) fields of the step (kernels.cu)
      void update_fields(const bool external);

      // writes the halo exchange times of every rank to the log and screen (halo.cpp)
      void report_halo_times();

      // orders the neighbours of each atom in the device exchange list by atom key (kernels.cu)
      void sort_neighbours_by_key();

      // thermal noise from the atom keys (thermal.cu)
      void initialize_noise();
      void free_noise();
      void update_external_fields(); // external fields of one time step, incl. thermal and HAMR fields

      // stops with an error message for unsupported simulations (initialize.cpp)
      void unsupported(const std::string& feature);

   } // end of internal namespace

} // end of mgpu namespace

#endif //MGPU_INTERNAL_H_
