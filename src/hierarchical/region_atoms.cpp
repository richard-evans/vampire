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
//   Table of the atoms of each level 0 cell, held once per shared memory region,
//   for the atomistic (near field) part of the hierarchical dipole tensors
//   (calc_inter and calc_intra). It replaces, for the hierarchical solver, the
//   per-cell atom lists of the cells module (an inner vector for every cell in
//   the system on every process) and the all-process gather of
//   dipole::internal::initialise_atomistic_cell_data.
//
//   A process needs every atom of each cell within the atomistic cutoff of one
//   of its own cells with atoms (calc_inter). The table of a region holds:
//
//      - the magnetic atoms of all processes in the region, gathered to the
//        first process of the region;
//      - the magnetic atoms of other regions in cells whose centre lies within
//        the cutoff of the box bounding the region's cell centres, sent point to
//        point between region leaders. Regions further apart exchange nothing.
//
//   The table is sorted by cell, then by source process, then by the atom's
//   index on that process, which is the order in which initialise_atomistic_
//   cell_data delivers them, so that the tensor sums are unchanged.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "cells.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"

// hierarchical module headers
#include "internal.hpp"

namespace hierarchical{
namespace internal{

region_atoms_t region_atoms; // the table of the region of this process

namespace{

   // one record per atom: cell, source process, index on that process, x, y, z
   constexpr int record_size = 6;

   // box bounding a set of points, empty if min > max
   struct box_t{
      double min[3] = { std::numeric_limits<double>::max(),  std::numeric_limits<double>::max(),  std::numeric_limits<double>::max()};
      double max[3] = {-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};
      bool empty() const { return min[0] > max[0]; }
      void add(const double* p){ for(int i = 0; i < 3; i++){ min[i] = std::min(min[i], p[i]); max[i] = std::max(max[i], p[i]); } }
   };

   // square of the distance from point p to box b (zero inside)
   double distance_sq(const double* p, const box_t& b){
      double d2 = 0.0;
      for(int i = 0; i < 3; i++){
         const double d = std::max({b.min[i] - p[i], 0.0, p[i] - b.max[i]});
         d2 += d * d;
      }
      return d2;
   }

   // square of the distance between boxes a and b (zero if they overlap)
   double distance_sq(const box_t& a, const box_t& b){
      double d2 = 0.0;
      for(int i = 0; i < 3; i++){
         const double d = std::max({a.min[i] - b.max[i], 0.0, b.min[i] - a.max[i]});
         d2 += d * d;
      }
      return d2;
   }

}

//------------------------------------------------------------------------------
// Builds the table of the region (collective over all processes). cutoff is the
// atomistic cutoff (Angstroms) between cell centres in cells::pos_and_mom_array.
//------------------------------------------------------------------------------
void region_atoms_t::build(const double cutoff,
                           const std::vector<double>& atom_coords_x,
                           const std::vector<double>& atom_coords_y,
                           const std::vector<double>& atom_coords_z){

   vutil::vtimer_t timer;
   timer.start();

   vmpi::initialise_shared_memory(); // region communicators (no shared array may have been allocated yet)
   free();

   // records of this process's magnetic atoms (the atoms of the cells module's per-cell lists), in index order
   std::vector<double> records;
   for(int atom = 0; atom < vmpi::num_local_atoms; atom++){
      if(mp::material[atoms::type_array[atom]].non_magnetic != 0) continue;
      records.insert(records.end(), {double(cells::atom_cell_id_array[atom]), double(vmpi::my_rank), double(atom),
                                     atom_coords_x[atom], atom_coords_y[atom], atom_coords_z[atom]});
   }

   // statistics for the log (first process of each region)
   int num_region_atoms = 0;
   int num_halo_atoms = 0;
   int num_neighbours = 0;

   #ifdef MPICF

      //---------------------------------------------------------------------------
      // gather the records of the region to its first process
      //---------------------------------------------------------------------------
      const bool leader = (vmpi::my_shared_rank == 0);
      int num_values = records.size();
      std::vector<int> counts(leader ? vmpi::num_shared_ranks : 0);
      std::vector<int> displacements(leader ? vmpi::num_shared_ranks : 0);
      MPI_Gather(&num_values, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, vmpi::shared_comm);
      std::vector<double> region_records;
      if(leader){
         std::partial_sum(counts.begin(), counts.end() - 1, displacements.begin() + 1);
         region_records.resize(displacements.back() + counts.back());
      }
      MPI_Gatherv(records.data(), num_values, MPI_DOUBLE, region_records.data(), counts.data(), displacements.data(),
                  MPI_DOUBLE, 0, vmpi::shared_comm);
      std::vector<double>().swap(records);

      //---------------------------------------------------------------------------
      // exchange the atoms near each other region between region leaders
      //---------------------------------------------------------------------------
      if(leader){

         num_region_atoms = region_records.size() / record_size;

         // box bounding the centres of the region's cells with atoms
         box_t box;
         for(int r = 0; r < num_region_atoms; r++) box.add(&cells::pos_and_mom_array[4 * int(region_records[record_size * r])]);

         int num_leaders = 1;
         int my_leader = 0;
         MPI_Comm_size(vmpi::leader_comm, &num_leaders);
         MPI_Comm_rank(vmpi::leader_comm, &my_leader);
         std::vector<box_t> boxes(num_leaders);
         MPI_Allgather(&box, sizeof(box_t), MPI_BYTE, boxes.data(), sizeof(box_t), MPI_BYTE, vmpi::leader_comm);

         // regions within the cutoff of this one (a symmetric relation, so each side expects a message from the
         // other, possibly empty), with a margin of 1 A so that a cell exactly at the cutoff is never missed
         const double range = cutoff + 1.0;
         const double range_sq = range * range;
         std::vector<int> neighbours;
         for(int q = 0; q < num_leaders; q++){
            if(q != my_leader && !box.empty() && !boxes[q].empty() && distance_sq(box, boxes[q]) <= range_sq) neighbours.push_back(q);
         }
         num_neighbours = neighbours.size();

         // records of this region's atoms in cells within range of each neighbour's box
         std::vector< std::vector<double> > send(num_neighbours);
         for(int n = 0; n < num_neighbours; n++){
            for(int r = 0; r < num_region_atoms; r++){
               const double* centre = &cells::pos_and_mom_array[4 * int(region_records[record_size * r])];
               if(distance_sq(centre, boxes[neighbours[n]]) <= range_sq){
                  send[n].insert(send[n].end(), &region_records[record_size * r], &region_records[record_size * (r + 1)]);
               }
            }
         }

         // sizes, then records
         std::vector<int> send_sizes(num_neighbours), recv_sizes(num_neighbours);
         std::vector<MPI_Request> requests(2 * num_neighbours);
         for(int n = 0; n < num_neighbours; n++){
            send_sizes[n] = send[n].size();
            MPI_Irecv(&recv_sizes[n], 1, MPI_INT, neighbours[n], 0, vmpi::leader_comm, &requests[n]);
            MPI_Isend(&send_sizes[n], 1, MPI_INT, neighbours[n], 0, vmpi::leader_comm, &requests[num_neighbours + n]);
         }
         MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

         std::vector< std::vector<double> > recv(num_neighbours);
         for(int n = 0; n < num_neighbours; n++){
            recv[n].resize(recv_sizes[n]);
            MPI_Irecv(recv[n].data(), recv_sizes[n], MPI_DOUBLE, neighbours[n], 1, vmpi::leader_comm, &requests[n]);
            MPI_Isend(send[n].data(), send_sizes[n], MPI_DOUBLE, neighbours[n], 1, vmpi::leader_comm, &requests[num_neighbours + n]);
         }
         MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
         std::vector< std::vector<double> >().swap(send);

         for(int n = 0; n < num_neighbours; n++){
            region_records.insert(region_records.end(), recv[n].begin(), recv[n].end());
            num_halo_atoms += recv[n].size() / record_size;
         }

      }

      records.swap(region_records);

   #else
      num_region_atoms = records.size() / record_size;
   #endif

   //---------------------------------------------------------------------------
   // sort by cell, source process and index, and fill the shared table
   //---------------------------------------------------------------------------
   uint64_t sizes[2] = {0, 0}; // atoms, cells
   std::vector<int> order;
   std::vector<int> present_cells;
   std::vector<int> present_starts;
   if(vmpi::my_shared_rank == 0){

      const int num_records = records.size() / record_size;
      order.resize(num_records);
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&records](const int a, const int b){
         const double* ra = &records[record_size * a];
         const double* rb = &records[record_size * b];
         return std::lexicographical_compare(ra, ra + 3, rb, rb + 3);
      });

      for(int i = 0; i < num_records; i++){
         const int cell = records[record_size * order[i]];
         if(present_cells.empty() || present_cells.back() != cell){
            present_cells.push_back(cell);
            present_starts.push_back(i);
         }
      }
      present_starts.push_back(num_records);
      sizes[0] = num_records;
      sizes[1] = present_cells.size();

   }
   #ifdef MPICF
      MPI_Bcast(sizes, 2, MPI_UINT64_T, 0, vmpi::shared_comm);
   #endif

   coords.allocate(3 * sizes[0]);
   cell_ids.allocate(sizes[1]);
   cell_starts.allocate(sizes[1] + 1);
   if(coords.is_writer()){
      for(uint64_t i = 0; i < sizes[0]; i++){
         for(int c = 0; c < 3; c++) coords[3 * i + c] = records[record_size * order[i] + 3 + c];
      }
      for(uint64_t k = 0; k < sizes[1]; k++){
         cell_ids[k] = present_cells[k];
         cell_starts[k] = present_starts[k];
      }
      cell_starts[sizes[1]] = present_starts[sizes[1]];
   }
   coords.synchronise();
   cell_ids.synchronise();
   cell_starts.synchronise();

   timer.stop();

   const double mb = double(coords.size() * sizeof(double) + (cell_ids.size() + cell_starts.size()) * sizeof(int)) / 1.0e6;
   const double max_atoms = vmpi::reduce_max(double(num_region_atoms));
   const double max_halo = vmpi::reduce_max(double(num_halo_atoms));
   const double max_neighbours = vmpi::reduce_max(double(num_neighbours));
   const double max_mb = vmpi::reduce_max(mb);
   zlog << zTs() << "\tHierarchical atom table per region (max over regions): " << max_atoms << " atoms of the region, "
        << max_halo << " from " << max_neighbours << " neighbouring regions, " << max_mb << " MB shared; built in "
        << timer.elapsed_time() << " s" << std::endl;

   return;

}

//------------------------------------------------------------------------------
// Index of cell in the table, or -1 if the table does not hold it
//------------------------------------------------------------------------------
int region_atoms_t::find(const int cell) const{
   const int* begin = cell_ids.data();
   const int* end = begin + cell_ids.size();
   const int* it = std::lower_bound(begin, end, cell);
   return (it != end && *it == cell) ? int(it - begin) : -1;
}

//------------------------------------------------------------------------------
// Frees the table (collective)
//------------------------------------------------------------------------------
void region_atoms_t::free(){
   coords.free();
   cell_ids.free();
   cell_starts.free();
}

} // end of internal namespace
} // end of hierarchical namespace
