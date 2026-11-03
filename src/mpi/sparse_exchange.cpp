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
#include <string>
#include <vector>

// Vampire headers
#include "errors.hpp"
#include "vmpi.hpp"

namespace vmpi{

#ifdef MPICF

//------------------------------------------------------------------------------
// Function to set up an exchange of num_components values per index. indices[p]
// lists the indices this rank sends to (lists_t::sent) or receives from
// (lists_t::received) rank p of comm; the lists for the other side are found
// here. Persistent requests are set up for one message per peer. Collective
// over comm.
//------------------------------------------------------------------------------
void sparse_exchange_t::initialise(MPI_Comm user_comm, const std::vector< std::vector<int> >& indices, lists_t lists,
                                   int components){

   free();

   MPI_Comm_dup(user_comm, &comm);
   int num_ranks = 1;
   int rank = 0;
   MPI_Comm_size(comm, &num_ranks);
   MPI_Comm_rank(comm, &rank);
   num_components = components;

   if(int(indices.size()) != num_ranks){
      err::fatal("sparse exchange given " + std::to_string(indices.size()) + " index lists for a communicator of " +
                 std::to_string(num_ranks) + " ranks");
   }
   if(!indices[rank].empty()){
      err::fatal("sparse exchange given " + std::to_string(indices[rank].size()) + " indices to exchange with its own rank");
   }

   // this rank's side: the lists given, peer by peer in rank order
   std::vector<int>& my_peers   = (lists == lists_t::sent) ? send_peers   : recv_peers;
   std::vector<int>& my_offsets = (lists == lists_t::sent) ? send_offsets : recv_offsets;
   std::vector<int>& my_list    = (lists == lists_t::sent) ? send_list    : recv_list;
   std::vector<int>& peers      = (lists == lists_t::sent) ? recv_peers   : send_peers;
   std::vector<int>& offsets    = (lists == lists_t::sent) ? recv_offsets : send_offsets;
   std::vector<int>& list       = (lists == lists_t::sent) ? recv_list    : send_list;

   std::vector<int> my_counts(num_ranks);
   std::vector<int> my_displacements(num_ranks);
   my_offsets.assign(1, 0);
   for(int p = 0; p < num_ranks; p++){
      my_counts[p] = indices[p].size();
      my_displacements[p] = my_list.size();
      if(my_counts[p] > 0){
         my_peers.push_back(p);
         my_list.insert(my_list.end(), indices[p].begin(), indices[p].end());
         my_offsets.push_back(my_list.size());
      }
   }

   // the other side: every rank tells each peer how many indices, then which
   std::vector<int> counts(num_ranks);
   std::vector<int> displacements(num_ranks);
   MPI_Alltoall(my_counts.data(), 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
   offsets.assign(1, 0);
   for(int p = 0; p < num_ranks; p++){
      displacements[p] = offsets.back();
      if(counts[p] > 0){
         peers.push_back(p);
         offsets.push_back(offsets.back() + counts[p]);
      }
   }
   list.resize(offsets.back());
   MPI_Alltoallv(my_list.data(), my_counts.data(), my_displacements.data(), MPI_INT,
                 list.data(), counts.data(), displacements.data(), MPI_INT, comm);

   // one persistent receive, then one persistent send, per peer
   send_buffer.resize(uint64_t(num_components) * send_list.size());
   recv_buffer.resize(uint64_t(num_components) * recv_list.size());
   requests.resize(recv_peers.size() + send_peers.size());
   const int tag = 0;
   for(size_t p = 0; p < recv_peers.size(); p++){
      const int count = num_components * (recv_offsets[p+1] - recv_offsets[p]);
      MPI_Recv_init(&recv_buffer[uint64_t(num_components) * recv_offsets[p]], count, MPI_DOUBLE, recv_peers[p], tag,
                    comm, &requests[p]);
   }
   for(size_t p = 0; p < send_peers.size(); p++){
      const int count = num_components * (send_offsets[p+1] - send_offsets[p]);
      MPI_Send_init(&send_buffer[uint64_t(num_components) * send_offsets[p]], count, MPI_DOUBLE, send_peers[p], tag,
                    comm, &requests[recv_peers.size() + p]);
   }

   return;

}

//------------------------------------------------------------------------------
// Function to send arrays[c][i] for each index i sent, and to assign or add
// each value received into arrays[c] at its index. All values are packed
// before any received value is combined, and received values are combined in
// order of source rank. Collective over the communicator.
//------------------------------------------------------------------------------
void sparse_exchange_t::exchange(const std::vector<double*>& arrays, combine_t combine){

   if(int(arrays.size()) != num_components){
      err::fatal("sparse exchange given " + std::to_string(arrays.size()) + " arrays for " +
                 std::to_string(num_components) + " components");
   }

   const int num_recvs = recv_peers.size();
   if(num_recvs > 0) MPI_Startall(num_recvs, requests.data());

   const uint64_t num_send = send_list.size();
   for(uint64_t k = 0; k < num_send; k++){
      for(int c = 0; c < num_components; c++) send_buffer[k * num_components + c] = arrays[c][send_list[k]];
   }

   const int num_sends = send_peers.size();
   if(num_sends > 0) MPI_Startall(num_sends, requests.data() + num_recvs);
   if(!requests.empty()) MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

   const uint64_t num_recv = recv_list.size();
   switch(combine){
      case combine_t::assign:
         for(uint64_t k = 0; k < num_recv; k++){
            for(int c = 0; c < num_components; c++) arrays[c][recv_list[k]] = recv_buffer[k * num_components + c];
         }
         break;
      case combine_t::add:
         for(uint64_t k = 0; k < num_recv; k++){
            for(int c = 0; c < num_components; c++) arrays[c][recv_list[k]] += recv_buffer[k * num_components + c];
         }
         break;
   }

   return;

}

//------------------------------------------------------------------------------
// Function to free the requests and the communicator. Collective over the
// communicator; does nothing if the exchange is not set up, or after
// MPI_Finalize (when the object outlives MPI).
//------------------------------------------------------------------------------
void sparse_exchange_t::free(){

   int finalised = 0;
   MPI_Finalized(&finalised);
   if(comm != MPI_COMM_NULL && !finalised){
      for(MPI_Request& request : requests) MPI_Request_free(&request);
      MPI_Comm_free(&comm);
   }
   comm = MPI_COMM_NULL;
   num_components = 0;
   send_peers.clear();
   send_offsets.clear();
   send_list.clear();
   recv_peers.clear();
   recv_offsets.clear();
   recv_list.clear();
   send_buffer.clear();
   recv_buffer.clear();
   requests.clear();

   return;

}

#endif

} // end of namespace vmpi
