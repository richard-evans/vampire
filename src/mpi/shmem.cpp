//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <functional>
#include <string>
#include <utility>
#include <vector>

// Vampire headers
#include "vio.hpp"
#include "vmpi.hpp"

namespace vmpi{

// shared memory communicators
#ifdef MPICF
   MPI_Comm shared_comm = MPI_COMM_NULL;
   MPI_Comm leader_comm = MPI_COMM_NULL;
#endif
int my_shared_rank = 0;
int num_shared_ranks = 1;

namespace{
   bool shared_memory_initialised = false; // true once the shared memory communicators exist

   // arrays currently allocated, in allocation order (the same on every rank, since allocation is
   // collective), with a function freeing each one, so that finalise_shared_memory can free any left
   std::vector< std::pair<const void*, std::function<void()> > > allocated_arrays;

   void register_array(const void* array, std::function<void()> free_array){
      allocated_arrays.push_back(std::make_pair(array, free_array));
   }

   void unregister_array(const void* array){
      for(size_t i = 0; i < allocated_arrays.size(); i++){
         if(allocated_arrays[i].first == array){
            allocated_arrays.erase(allocated_arrays.begin() + i);
            return;
         }
      }
   }
}

//------------------------------------------------------------------------------
// Function to create the communicators for ranks sharing a region of memory
//
// shared_comm holds the ranks in the same NUMA region when the MPI library can
// identify it (MPI 4 hardware guided split or Open MPI), otherwise the ranks on
// the same node. leader_comm holds rank 0 of each shared_comm, for
// communication between regions, and is MPI_COMM_NULL on all other ranks.
// Collective over MPI_COMM_WORLD; calls after the first do nothing.
//------------------------------------------------------------------------------
void initialise_shared_memory(){

   if(shared_memory_initialised) return;

   #ifdef MPICF

      // split by NUMA region where the library supports it
      MPI_Comm numa_comm = MPI_COMM_NULL;
      #if MPI_VERSION >= 4
         MPI_Info info;
         MPI_Info_create(&info);
         MPI_Info_set(info, "mpi_hw_resource_type", "NUMANode");
         MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_HW_GUIDED, vmpi::my_rank, info, &numa_comm);
         MPI_Info_free(&info);
      #elif defined(OPEN_MPI)
         MPI_Comm_split_type(MPI_COMM_WORLD, OMPI_COMM_TYPE_NUMA, vmpi::my_rank, MPI_INFO_NULL, &numa_comm);
      #endif

      // ranks on the same node
      MPI_Comm node_comm = MPI_COMM_NULL;
      MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, vmpi::my_rank, MPI_INFO_NULL, &node_comm);

      // The NUMA split is not usable if any rank has no region, or if every region holds a single rank
      // while a node holds several, which is what Open MPI gives when ranks are not bound to cores. Then
      // memory is shared across the node instead, on every rank.
      int sizes[3] = {1, 0, 0}; // NUMA region found, largest NUMA region, largest node
      if(numa_comm == MPI_COMM_NULL) sizes[0] = 0;
      else MPI_Comm_size(numa_comm, &sizes[1]);
      MPI_Comm_size(node_comm, &sizes[2]);
      sizes[0] = -sizes[0]; // so that a max reduction gives the minimum
      MPI_Allreduce(MPI_IN_PLACE, sizes, 3, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      const bool use_numa = (sizes[0] == -1) && !(sizes[1] == 1 && sizes[2] > 1);

      std::string region = "NUMA region";
      if(use_numa){
         vmpi::shared_comm = numa_comm;
         MPI_Comm_free(&node_comm);
      }
      else{
         if(numa_comm != MPI_COMM_NULL) MPI_Comm_free(&numa_comm);
         vmpi::shared_comm = node_comm;
         region = "node";
      }
      MPI_Comm_rank(vmpi::shared_comm, &vmpi::my_shared_rank);
      MPI_Comm_size(vmpi::shared_comm, &vmpi::num_shared_ranks);

      // communicator of the first rank in each region
      const int leader = (vmpi::my_shared_rank == 0);
      MPI_Comm_split(MPI_COMM_WORLD, leader ? 0 : MPI_UNDEFINED, vmpi::my_rank, &vmpi::leader_comm);

      int num_regions = leader;
      MPI_Allreduce(MPI_IN_PLACE, &num_regions, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      zlog << zTs() << "Shared memory initialised with " << num_regions << " regions (one per " << region << "), "
           << vmpi::num_shared_ranks << " ranks share memory with rank " << vmpi::my_rank << std::endl;

   #endif

   shared_memory_initialised = true;

   return;

}

//------------------------------------------------------------------------------
// Function to free the shared memory communicators, called before MPI_Finalize.
// Shared arrays still allocated are freed first, newest first, and logged.
//------------------------------------------------------------------------------
void finalise_shared_memory(){

   if(!shared_memory_initialised) return;

   if(!allocated_arrays.empty()){
      zlog << zTs() << "Freeing " << allocated_arrays.size() << " shared memory arrays at finalisation" << std::endl;
   }
   while(!allocated_arrays.empty()){
      std::function<void()> free_array = allocated_arrays.back().second;
      free_array(); // removes the entry
   }

   #ifdef MPICF
      if(vmpi::leader_comm != MPI_COMM_NULL) MPI_Comm_free(&vmpi::leader_comm);
      MPI_Comm_free(&vmpi::shared_comm);
   #endif

   vmpi::my_shared_rank = 0;
   vmpi::num_shared_ranks = 1;
   shared_memory_initialised = false;

   return;

}

//------------------------------------------------------------------------------
// Function to allocate a contiguous array of num_elements shared by all ranks in
// shared_comm, so that it is stored once per region. Any previous array is freed.
// The array is set to zero.
//
// Ranks read and write the array directly, and call synchronise() between
// writing and reading data written by other ranks. Collective over
// MPI_COMM_WORLD on first use, over shared_comm afterwards.
//------------------------------------------------------------------------------
template <typename T>
void shared_array_t<T>::allocate(uint64_t n){

   free();
   vmpi::initialise_shared_memory();

   num_elements = n;

   #ifdef MPICF

      // allocate the whole array on the first rank of the region so that it is contiguous
      const MPI_Aint bytes = (vmpi::my_shared_rank == 0) ? n * sizeof(T) : 0;
      void* base = nullptr;
      MPI_Win_allocate_shared(bytes, sizeof(T), MPI_INFO_NULL, vmpi::shared_comm, &base, &window);

      // get the address of the array in this rank's address space
      MPI_Aint region_bytes = 0;
      int displacement_unit = 0;
      MPI_Win_shared_query(window, 0, &region_bytes, &displacement_unit, &base);
      ptr = static_cast<T*>(base);

      // passive target epoch for the lifetime of the array, so that ranks can use plain loads and stores
      MPI_Win_lock_all(MPI_MODE_NOCHECK, window);

      if(vmpi::my_shared_rank == 0) std::fill(ptr, ptr + n, T(0));
      synchronise();

   #else

      serial_data.assign(n, T(0));
      ptr = serial_data.data();

   #endif

   register_array(this, [this](){ this->free(); });

   return;

}

//------------------------------------------------------------------------------
// Function to make writes to the array by any rank in shared_comm visible to
// all ranks in shared_comm. Collective over shared_comm.
//------------------------------------------------------------------------------
template <typename T>
void shared_array_t<T>::synchronise(){

   #ifdef MPICF
      MPI_Win_sync(window); // complete this rank's writes
      MPI_Barrier(vmpi::shared_comm);
      MPI_Win_sync(window); // see the writes of the other ranks
   #endif

   return;

}

//------------------------------------------------------------------------------
// Function to free the array. Collective over shared_comm; does nothing if the
// array is not allocated.
//------------------------------------------------------------------------------
template <typename T>
void shared_array_t<T>::free(){

   #ifdef MPICF
      if(window != MPI_WIN_NULL){
         MPI_Win_unlock_all(window);
         MPI_Win_free(&window);
      }
   #else
      std::vector<T>().swap(serial_data);
   #endif

   unregister_array(this);
   ptr = nullptr;
   num_elements = 0;

   return;

}

//------------------------------------------------------------------------------
// Function to order this rank's loads and stores to the array with respect to
// later communication, without synchronising with other ranks.
//------------------------------------------------------------------------------
template <typename T>
void shared_array_t<T>::memory_barrier(){
   #ifdef MPICF
      MPI_Win_sync(window);
   #endif
   return;
}

// types available for shared arrays
template class shared_array_t<double>;
template class shared_array_t<int>;
template class shared_array_t<uint64_t>;

//------------------------------------------------------------------------------
// Function to set up a shared sum of num_components values for each of this
// rank's indices. The first rank of each region gathers the indices of every
// rank in the region. Collective over MPI_COMM_WORLD.
//------------------------------------------------------------------------------
void shared_sum_t::initialise(const std::vector<int>& indices, int components){

   vmpi::initialise_shared_memory();

   local_indices = indices;
   num_components = components;

   #ifdef MPICF
      const int num_indices = local_indices.size();
      const bool leader = (vmpi::my_shared_rank == 0);
      if(leader){
         index_counts.resize(vmpi::num_shared_ranks);
         index_displacements.resize(vmpi::num_shared_ranks);
         value_counts.resize(vmpi::num_shared_ranks);
         value_displacements.resize(vmpi::num_shared_ranks);
      }
      MPI_Gather(&num_indices, 1, MPI_INT, index_counts.data(), 1, MPI_INT, 0, vmpi::shared_comm);
      int total = 0;
      if(leader){
         for(int r = 0; r < vmpi::num_shared_ranks; r++){
            index_displacements[r] = total;
            value_counts[r] = num_components * index_counts[r];
            value_displacements[r] = num_components * total;
            total += index_counts[r];
         }
         region_indices.resize(total);
         region_values.resize(num_components * uint64_t(total));
      }
      MPI_Gatherv(local_indices.data(), num_indices, MPI_INT, region_indices.data(), index_counts.data(),
                  index_displacements.data(), MPI_INT, 0, vmpi::shared_comm);
   #endif

   return;

}

//------------------------------------------------------------------------------
// Function to sum values[c*n + i] (component c of this rank's index i, n
// indices) over all ranks into totals[c]. Collective over MPI_COMM_WORLD.
//
// Every rank calls this only after it has finished reading the totals of the
// previous sum, so once the gather completes the first rank of the region may
// overwrite them.
//------------------------------------------------------------------------------
void shared_sum_t::sum(const std::vector<double>& values, const std::vector< shared_array_t<double>* >& totals){

   const uint64_t num_indices = local_indices.size();

   #ifdef MPICF

      // complete this rank's reads of the previous totals before the first rank of the region writes them
      for(shared_array_t<double>* total : totals) total->memory_barrier();

      MPI_Gatherv(values.data(), num_components * num_indices, MPI_DOUBLE, region_values.data(), value_counts.data(),
                  value_displacements.data(), MPI_DOUBLE, 0, vmpi::shared_comm);

      // the first rank of each region adds the values of its ranks, then the regions are summed
      if(vmpi::my_shared_rank == 0){
         for(int c = 0; c < num_components; c++){
            shared_array_t<double>& total = *totals[c];
            std::fill(total.data(), total.data() + total.size(), 0.0);
            for(int r = 0; r < vmpi::num_shared_ranks; r++){
               const int n = index_counts[r];
               const int* index = &region_indices[index_displacements[r]];
               const double* value = &region_values[value_displacements[r] + uint64_t(c) * n];
               for(int i = 0; i < n; i++) total[index[i]] += value[i];
            }
         }
         int num_regions = 1;
         MPI_Comm_size(vmpi::leader_comm, &num_regions);
         if(num_regions > 1){
            for(shared_array_t<double>* total : totals){
               MPI_Allreduce(MPI_IN_PLACE, total->data(), total->size(), MPI_DOUBLE, MPI_SUM, vmpi::leader_comm);
            }
         }
      }

      // make the totals visible to every rank in the region (one barrier for all arrays)
      for(shared_array_t<double>* total : totals) total->memory_barrier();
      MPI_Barrier(vmpi::shared_comm);
      for(shared_array_t<double>* total : totals) total->memory_barrier();

   #else

      for(int c = 0; c < num_components; c++){
         shared_array_t<double>& total = *totals[c];
         std::fill(total.data(), total.data() + total.size(), 0.0);
         for(uint64_t i = 0; i < num_indices; i++) total[local_indices[i]] += values[c * num_indices + i];
      }

   #endif

   return;

}

} // end of namespace vmpi
