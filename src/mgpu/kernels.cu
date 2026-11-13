//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) S Jenkins 2026. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "errors.hpp"
#include "hamr.hpp"
#include "material.hpp"
#include "mgpu.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// CUDA module headers
#include "../cuda/cuda_utils.hpp"
#include "../cuda/data.hpp"
#include "../cuda/exchange_fields.hpp"
#include "../cuda/internal.hpp"
#include "../cuda/llg_heun.hpp"
#include "../cuda/typedefs.hpp"

// mgpu module headers
#include "internal.hpp"

namespace cu = ::vcuda::internal;

namespace mgpu {

namespace internal {

// Halo spins a rank sends to itself (periodic boundaries in a direction that
// is not divided between ranks) are copied on the device from their source
// atom (d_self_source) to the halo atom (d_self_halo); only the spins of the
// other ranks go through the host: the atoms sent (d_send_atoms) are packed
// into d_send_buffer, and the atoms received (d_recv_atoms) are unpacked from
// d_recv_buffer.
int *d_self_source = nullptr;
int *d_self_halo = nullptr;
int *d_send_atoms = nullptr;
int *d_recv_atoms = nullptr;
cu::cu_real_t *d_send_buffer = nullptr;
cu::cu_real_t *d_recv_buffer = nullptr;
int num_self = 0, num_remote_send = 0, num_remote_recv = 0; // atoms currently copied
bool device_only_halo = false;
std::vector<int> peer_send_start, peer_send_num, peer_recv_start, peer_recv_num;

// full halo lists, and for each entry whether its atom was integrated at the
// last choice of atoms (select_halo_atoms)
std::vector<int> all_send_atoms, all_recv_atoms;
std::vector<char> send_integrated, recv_integrated;

// Stream for the halo update.
cudaStream_t halo_stream;
cudaEvent_t step_done;

// Events timing the halo copies on the GPU and the waiting for the default
// stream to return the received spins; the times of one update are added after
// the next update has been packed
cudaEvent_t pack_start, pack_end, unpack_start, unpack_end, wait_start,
      wait_end, next_wait_start;
bool unpack_pending = false;
double copy_time = 0.0;
double wait_time = 0.0;

double elapsed(cudaEvent_t start, cudaEvent_t end) {
   float ms = 0.0f;
   cudaEventElapsedTime(&ms, start, end);
   return 1.0e-3 * ms;
}

//----------------------------------------------------------------
// Kernels to pack and unpack spins of listed atoms into buffers
//----------------------------------------------------------------
__global__ void pack_spins_kernel(const int *atoms, const int n,
                                                   const cu::cu_real_t *x,
                                                   const cu::cu_real_t *y,
                                                   const cu::cu_real_t *z,
                                                   cu::cu_real_t *buffer) {
   for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
          i += blockDim.x * gridDim.x) {
      const int atom = atoms[i];
      buffer[3 * i + 0] = x[atom];
      buffer[3 * i + 1] = y[atom];
      buffer[3 * i + 2] = z[atom];
   }
}

// copies the spins of the source atoms to their halo atoms
__global__ void copy_spins_kernel(const int *source, const int *halo,
                                                   const int n, cu::cu_real_t *x,
                                                   cu::cu_real_t *y, cu::cu_real_t *z) {
   for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
          i += blockDim.x * gridDim.x) {
      const int from = source[i];
      const int to = halo[i];
      x[to] = x[from];
      y[to] = y[from];
      z[to] = z[from];
   }
}

__global__ void unpack_spins_kernel(const int *atoms, const int n,
                                                      const cu::cu_real_t *buffer,
                                                      cu::cu_real_t *x, cu::cu_real_t *y,
                                                      cu::cu_real_t *z) {
   for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
          i += blockDim.x * gridDim.x) {
      const int atom = atoms[i];
      x[atom] = buffer[3 * i + 0];
      y[atom] = buffer[3 * i + 1];
      z[atom] = buffer[3 * i + 2];
   }
}

int device_count() {
   int n = 0;
   if (cudaGetDeviceCount(&n) != cudaSuccess)
      return 0;
   return n;
}

void device_pci_id(int device, char *id, int size) {
   if (cudaDeviceGetPCIBusId(id, size, device) != cudaSuccess)
      id[0] = '\0';
}

void allocate_device_halo(const std::vector<int> &send_atoms,
                                       const std::vector<int> &recv_atoms) {

   real_size = sizeof(cu::cu_real_t);
   cudaMalloc((void **)&d_self_source,
                   std::max<size_t>(send_atoms.size(), 1) * sizeof(int));
   cudaMalloc((void **)&d_self_halo,
                   std::max<size_t>(recv_atoms.size(), 1) * sizeof(int));
   cudaMalloc((void **)&d_send_atoms,
                   std::max<size_t>(send_atoms.size(), 1) * sizeof(int));
   cudaMalloc((void **)&d_recv_atoms,
                   std::max<size_t>(recv_atoms.size(), 1) * sizeof(int));
   cudaMalloc((void **)&d_send_buffer,
                   std::max<size_t>(3 * send_atoms.size(), 1) *
                         sizeof(cu::cu_real_t));
   cudaMalloc((void **)&d_recv_buffer,
                   std::max<size_t>(3 * recv_atoms.size(), 1) *
                         sizeof(cu::cu_real_t));
   cudaMallocHost(&h_send_buffer, std::max<size_t>(3 * send_atoms.size(), 1) *
                                                       sizeof(cu::cu_real_t));
   cudaMallocHost(&h_recv_buffer, std::max<size_t>(3 * recv_atoms.size(), 1) *
                                                       sizeof(cu::cu_real_t));
   const int me = vmpi::my_rank;
   if (vmpi::send_num_array[me] != vmpi::recv_num_array[me])
      err::fatal("Error: rank " + std::to_string(me) +
                      " sends a different number of halo atoms to itself than it "
                      "receives");
   device_only_halo =
         vmpi::send_num_array[me] == int(send_atoms.size()) &&
         vmpi::recv_num_array[me] == int(recv_atoms.size());
   all_send_atoms = send_atoms;
   all_recv_atoms = recv_atoms;
   send_integrated.assign(send_atoms.size(), 1);
   recv_integrated.assign(recv_atoms.size(), 1);
   select_halo_atoms();
   cudaStreamCreateWithFlags(&halo_stream, cudaStreamNonBlocking);
   cudaEventCreateWithFlags(&step_done, cudaEventDisableTiming);
   cudaEventCreate(&pack_start);
   cudaEventCreate(&pack_end);
   cudaEventCreate(&unpack_start);
   cudaEventCreate(&unpack_end);
   cudaEventCreate(&wait_start);
   cudaEventCreate(&wait_end);
   cudaEventCreate(&next_wait_start);

   if (cudaGetLastError() != cudaSuccess) {
      err::fatal("Error: unable to allocate multi-GPU halo buffers on rank " +
                      std::to_string(vmpi::my_rank));
   }
}

void free_device_halo() {
   cudaStreamDestroy(halo_stream);
   cudaEventDestroy(step_done);
   cudaEventDestroy(pack_start);
   cudaEventDestroy(pack_end);
   cudaEventDestroy(unpack_start);
   cudaEventDestroy(unpack_end);
   cudaEventDestroy(wait_start);
   cudaEventDestroy(wait_end);
   cudaEventDestroy(next_wait_start);
   cudaFree(d_self_source);
   cudaFree(d_self_halo);
   cudaFree(d_send_atoms);
   cudaFree(d_recv_atoms);
   cudaFree(d_send_buffer);
   cudaFree(d_recv_buffer);
   cudaFreeHost(h_send_buffer);
   cudaFreeHost(h_recv_buffer);
   d_self_source = d_self_halo = d_send_atoms = d_recv_atoms = nullptr;
   d_send_buffer = d_recv_buffer = nullptr;
   h_send_buffer = h_recv_buffer = nullptr;
}

//---------------------------------------------------------------------------
// Function to choose the halo atoms to copy. Without the HAMR active window
// every atom is copied. With it, a halo atom is copied while the grain of its
// atom is integrated, and at the step after the grain leaves the window, as
// its last corrector step changed the spin after the last copy; any other
// halo spin equals its source spin already. The grains are the same on every
// rank, so the sender and the receiver of each spin make the same choice.
//---------------------------------------------------------------------------
void select_halo_atoms() {

   const bool window = cu::d_active_atoms != nullptr;
   const int me = vmpi::my_rank;
   std::vector<int> self_source, self_halo, send, recv;
   peer_send_start.assign(vmpi::num_processors, 0);
   peer_send_num.assign(vmpi::num_processors, 0);
   peer_recv_start.assign(vmpi::num_processors, 0);
   peer_recv_num.assign(vmpi::num_processors, 0);

   for (int p = 0; p < vmpi::num_processors; p++) {
      peer_send_start[p] = send.size();
      const int first = vmpi::send_start_index_array[p];
      for (int i = first; i < first + vmpi::send_num_array[p]; i++) {
         const char integrated = !window || hamr::atom_in_window(all_send_atoms[i]);
         const bool copy = integrated || send_integrated[i];
         send_integrated[i] = integrated;
         if (!copy)
            continue;
         if (p == me) {
            self_source.push_back(all_send_atoms[i]);
            self_halo.push_back(
                  all_recv_atoms[vmpi::recv_start_index_array[me] + i - first]);
         } else
            send.push_back(all_send_atoms[i]);
      }
      peer_send_num[p] = send.size() - peer_send_start[p];

      peer_recv_start[p] = recv.size();
      if (p == me)
         continue;
      const int first_recv = vmpi::recv_start_index_array[p];
      for (int i = first_recv; i < first_recv + vmpi::recv_num_array[p]; i++) {
         const char integrated = !window || hamr::atom_in_window(all_recv_atoms[i]);
         const bool copy = integrated || recv_integrated[i];
         recv_integrated[i] = integrated;
         if (copy)
            recv.push_back(all_recv_atoms[i]);
      }
      peer_recv_num[p] = recv.size() - peer_recv_start[p];
   }

   // the last halo update has finished reading the lists: the default stream
   // waited for it, and the next update waits for this copy
   num_self = self_source.size();
   num_remote_send = send.size();
   num_remote_recv = recv.size();
   cudaMemcpy(d_self_source, self_source.data(), num_self * sizeof(int),
                   cudaMemcpyHostToDevice);
   cudaMemcpy(d_self_halo, self_halo.data(), num_self * sizeof(int),
                   cudaMemcpyHostToDevice);
   cudaMemcpy(d_send_atoms, send.data(), num_remote_send * sizeof(int),
                   cudaMemcpyHostToDevice);
   cudaMemcpy(d_recv_atoms, recv.data(), num_remote_recv * sizeof(int),
                   cudaMemcpyHostToDevice);
   check_cuda_errors(__FILE__, __LINE__);
}

// chooses the halo atoms again when the grains of the window have changed at
// this or the previous step (counted by the hamr module, so that every rank
// chooses again at the same step, also when its own atoms are unaffected)
void update_halo_selection() {
   static bool window_seen = false;
   static uint64_t changes_seen = 0;
   static bool changed_last_step = false;
   const bool window = cu::d_active_atoms != nullptr;
   const bool changed = window != window_seen ||
                                  (window && hamr::window_changes() != changes_seen);
   if (changed || changed_last_step)
      select_halo_atoms();
   window_seen = window;
   changes_seen = hamr::window_changes();
   changed_last_step = changed;
}

// adds the times of the last update once it has finished
void add_pending_times() {
   if (!unpack_pending)
      return;
   copy_time += elapsed(unpack_start, unpack_end);
   wait_time += elapsed(wait_start, wait_end);
   unpack_pending = false;
}

void start_halo() {
   // spins must not be packed before the launched integration step has finished
   cudaEventRecord(step_done, 0);
   cudaStreamWaitEvent(halo_stream, step_done, 0);
   if (num_send == 0)
      return;
   cudaEventRecord(pack_start, halo_stream);
   if (num_self != 0) {
      copy_spins_kernel<<<cu::grid_size, cu::block_size, 0, halo_stream>>>(
            d_self_source, d_self_halo, num_self, cu::atoms::d_x_spin,
            cu::atoms::d_y_spin, cu::atoms::d_z_spin);
      check_cuda_errors(__FILE__, __LINE__);
   }
   if (num_remote_send != 0) {
      pack_spins_kernel<<<cu::grid_size, cu::block_size, 0, halo_stream>>>(
            d_send_atoms, num_remote_send, cu::atoms::d_x_spin,
            cu::atoms::d_y_spin, cu::atoms::d_z_spin, d_send_buffer);
      check_cuda_errors(__FILE__, __LINE__);
      cudaMemcpyAsync(h_send_buffer, d_send_buffer,
                              3 * num_remote_send * sizeof(cu::cu_real_t),
                              cudaMemcpyDeviceToHost, halo_stream);
   }
   cudaEventRecord(pack_end, halo_stream);
}

void finish_halo() {
   if (device_only_halo) {
      // every halo spin comes from this rank and was copied on the device, so the
      // host need not wait (the updates are not timed, as their events are read
      // only after a wait)
      exchange_halo();
      cudaEventRecord(unpack_end, halo_stream);
      cudaStreamWaitEvent(0, unpack_end, 0);
      return;
   }
   // everything launched on the default stream so far does not need halo spins
   // (recorded in a second event, as wait_start still times the last update)
   cudaEventRecord(next_wait_start, 0);
   // the send buffer is complete and the last unpack has read the receive buffer
   cudaStreamSynchronize(halo_stream);
   if (num_send != 0)
      copy_time += elapsed(pack_start, pack_end);
   add_pending_times();
   std::swap(wait_start, next_wait_start);
   exchange_halo();
   if (num_remote_recv != 0) {
      cudaEventRecord(unpack_start, halo_stream);
      cudaMemcpyAsync(d_recv_buffer, h_recv_buffer,
                              3 * num_remote_recv * sizeof(cu::cu_real_t),
                              cudaMemcpyHostToDevice, halo_stream);
      unpack_spins_kernel<<<cu::grid_size, cu::block_size, 0, halo_stream>>>(
            d_recv_atoms, num_remote_recv, d_recv_buffer, cu::atoms::d_x_spin,
            cu::atoms::d_y_spin, cu::atoms::d_z_spin);
      check_cuda_errors(__FILE__, __LINE__);
   } else
      cudaEventRecord(unpack_start, halo_stream);
   cudaEventRecord(unpack_end, halo_stream);
   cudaStreamWaitEvent(0, unpack_end, 0);
   cudaEventRecord(wait_end, 0);
   unpack_pending = true;
}

double halo_copy_time() {
   cudaDeviceSynchronize();
   add_pending_times();
   return copy_time;
}

double halo_wait_time() {
   cudaDeviceSynchronize();
   add_pending_times();
   return wait_time;
}

//---------------------------------------------------------------------------
// Function to order the neighbours of each atom in the exchange list of the
// CUDA module (built in local atom order) by atom key, so that the exchange
// field of every atom is summed in the same order for any decomposition.
// The list is only reordered: neighbours and exchange constants are unchanged.
//---------------------------------------------------------------------------
void sort_neighbours_by_key() {

   if (!cu::exchange::use_isotropic_kernel) {
      zlog << zTs()
             << "Warning: exchange is not isotropic, so the order of exchange sums "
                  "depends on the decomposition"
             << std::endl;
      return;
   }

   const int num_atoms = ::atoms::num_atoms;
   std::vector<int> start(num_atoms + 1);
   cudaMemcpy(start.data(), cu::exchange::d_nbr_start,
                   start.size() * sizeof(int), cudaMemcpyDeviceToHost);
   const int num_nbrs = start[num_atoms];
   if (num_nbrs == 0)
      return;

   std::vector<int> nbr(num_nbrs);
   std::vector<cu::cu_real_t> J(num_nbrs);
   cudaMemcpy(nbr.data(), cu::exchange::d_nbr_atom, num_nbrs * sizeof(int),
                   cudaMemcpyDeviceToHost);
   cudaMemcpy(J.data(), cu::exchange::d_nbr_J, num_nbrs * sizeof(cu::cu_real_t),
                   cudaMemcpyDeviceToHost);

   std::vector<std::pair<uint64_t, int>> order;
   std::vector<int> sorted_nbr(num_nbrs);
   std::vector<cu::cu_real_t> sorted_J(num_nbrs);
   for (int atom = 0; atom < num_atoms; atom++) {
      order.clear();
      for (int i = start[atom]; i < start[atom + 1]; i++)
         order.push_back(std::make_pair(atom_key[nbr[i]], i));
      std::sort(order.begin(), order.end());
      for (size_t n = 0; n < order.size(); n++) {
         sorted_nbr[start[atom] + n] = nbr[order[n].second];
         sorted_J[start[atom] + n] = J[order[n].second];
      }
   }

   cudaMemcpy(cu::exchange::d_nbr_atom, sorted_nbr.data(),
                   num_nbrs * sizeof(int), cudaMemcpyHostToDevice);
   cudaMemcpy(cu::exchange::d_nbr_J, sorted_J.data(),
                   num_nbrs * sizeof(cu::cu_real_t), cudaMemcpyHostToDevice);
   zlog << zTs() << "Exchange neighbours of each atom ordered by atom key"
          << std::endl;
}
//---------------------------------------------------------------------------
// Function to update the halo spins and the fields of the local atoms. Core
// atoms have no halo neighbours, so their spin fields (and the external
// fields, which need no spins) are calculated while the halo is updated, as
// in the CPU-only MPI version; boundary atoms follow once the halo spins have
// arrived.
//---------------------------------------------------------------------------
void update_fields(const bool external) {

   const int num_local_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;

   if (!cu::exchange::use_isotropic_kernel || (num_send == 0 && num_recv == 0)) {
      swap_halo();
      cu::update_spin_fields(num_local_atoms);
      if (external)
         update_external_fields(); // thermal noise from the atom keys (thermal.cu)
      return;
   }

   start_halo();
   cu::update_spin_fields(vmpi::num_core_atoms);
   if (external)
      update_external_fields();
   finish_halo();
   cu::update_spin_fields(num_local_atoms, vmpi::num_core_atoms);
}

} // namespace internal

//------------------------------------------------------------------------------
// Function to perform a single LLG Heun step: Only the local (core and
// boundary) atoms, which come first, are integrated and get spin fields; halo
// atoms only change when they are received, so after a step they hold the
// predicted spins of the last update, as in the CPU-only MPI version.
//------------------------------------------------------------------------------
void llg_heun() {

   if (!cu::llg::initialized)
      cu::llg::__llg_init();

   const size_t num_local_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;

   // atoms to integrate (HAMR active window) and halo atoms to copy
   cu::update_active_atoms();
   internal::update_halo_selection();

   // fields at the start of the step
   internal::update_fields(true);

   cu::llg::llg_heun_predictor_step<<<cu::grid_size, cu::block_size>>>(
         cu::atoms::d_materials, cu::llg::d_heun_params, cu::atoms::d_x_spin,
         cu::atoms::d_y_spin, cu::atoms::d_z_spin, cu::d_x_spin_field,
         cu::d_y_spin_field, cu::d_z_spin_field, cu::d_x_external_field,
         cu::d_y_external_field, cu::d_z_external_field, cu::llg::d_ds_x,
         cu::llg::d_ds_y, cu::llg::d_ds_z, ::mp::dt,
         cu::num_integrated_atoms(num_local_atoms), cu::d_active_atoms);
   check_cuda_errors(__FILE__, __LINE__);

   // fields of the predicted spins
   internal::update_fields(false);

   cu::llg::llg_heun_corrector_step<<<cu::grid_size, cu::block_size>>>(
         cu::atoms::d_materials, cu::llg::d_heun_params, cu::atoms::d_x_spin,
         cu::atoms::d_y_spin, cu::atoms::d_z_spin, cu::d_x_spin_field,
         cu::d_y_spin_field, cu::d_z_spin_field, cu::d_x_external_field,
         cu::d_y_external_field, cu::d_z_external_field, cu::llg::d_x_spin_buffer,
         cu::llg::d_y_spin_buffer, cu::llg::d_z_spin_buffer, cu::llg::d_ds_x,
         cu::llg::d_ds_y, cu::llg::d_ds_z, ::mp::dt,
         cu::num_integrated_atoms(num_local_atoms), cu::d_active_atoms);
   check_cuda_errors(__FILE__, __LINE__);
}

} // namespace mgpu
