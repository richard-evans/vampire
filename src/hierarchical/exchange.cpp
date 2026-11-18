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
//   Point-to-point exchange of the cell moments between shared memory regions,
//   chosen with hierarchical:moment-exchange = point-to-point. It replaces the
//   all-reduce of every cell (vmpi::shared_sum_t, mag.cpp) by:
//
//      - near cells (levels up to hierarchical:moment-exchange-level): each
//        region leader sends its partial moments to the region owning the
//        cell (the region with most of its magnetic atoms), and the owner
//        sends the totals to every region that reads the cell;
//      - far cells (the higher levels, one contiguous range): an all-reduce
//        between the region leaders.
//
//   A region needs correct totals only for the cells its processes read: the
//   cells in their interaction lists and their own level 0 cells.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "cells.hpp"
#include "dipole.hpp"
#include "errors.hpp"
#include "micromagnetic.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"

// hierarchical module headers
#include "internal.hpp"

// alias internal hierarchical namespace for brevity
namespace ha = hierarchical::internal;

namespace hierarchical{
namespace internal{

namespace{

   int near_cells_end = 0; // cells [0, near_cells_end) are exchanged point to point, the rest all-reduced

   // phases of the last exchange on this process (add, contribute, far and serve: region leaders only)
   enum exchange_phase_t { gather = 0, add = 1, contribute = 2, far = 3, serve = 4, barrier = 5, num_phases = 6 };
   double exchange_time[num_phases] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

   #ifdef MPICF

   const double bytes_per_cell = 3.0 * sizeof(double); // x, y and z moments

   vmpi::sparse_exchange_t contribute_exchange; // partial moments to the owner of each near cell (region leaders)
   vmpi::sparse_exchange_t serve_exchange;      // totals from the owner to each region reading the cell

   // gather of the contributions of the processes of a region to its leader, as in vmpi::shared_sum_t
   std::vector<int> region_indices;      // leader: contributing cells of every process, process by process
   std::vector<int> index_counts;        // leader: number of contributing cells of each process
   std::vector<int> index_displacements;
   std::vector<int> value_counts;        // leader: number of values (3 per cell) of each process
   std::vector<int> value_displacements;
   std::vector<double> region_values;    // leader: gathered values, x then y then z for each process

   std::vector<int> region_read_cells;   // leader: cells (all levels) read by the region, for the cross-check
   bool cross_checked = false;           // the first exchange is compared with the all-reduce

   // minimum, mean and maximum of a value over the processes of a communicator as min/mean/max
   std::string spread_over(const double value, MPI_Comm comm){
      int size = 1;
      MPI_Comm_size(comm, &size);
      double min = value, sum = value, max = value;
      MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, comm);
      MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
      MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, comm);
      std::ostringstream text;
      text << min << "/" << sum / double(size) << "/" << max;
      return text.str();
   }

   //---------------------------------------------------------------------------
   // Function to check both stages through the exchanges themselves, with
   // test values in place of moments: after the contribute stage the owner of
   // each near cell must hold one contribution from every region with atoms in
   // it, and after the serve stage every region must hold the owner's value of
   // each near cell it reads. Returns the number of wrong values. Region
   // leaders only.
   //---------------------------------------------------------------------------
   uint64_t check_exchange(const std::vector<int>& region_atoms, const std::vector<int>& region_reads,
                           const std::vector<int>& owner, const std::vector<int>& max_atoms,
                           const std::vector<int>& contributors, const int my_region){

      uint64_t wrong = 0;
      std::vector<double> x(near_cells_end, 0.0);
      std::vector<double> y(near_cells_end, 0.0);
      std::vector<double> z(near_cells_end, 0.0);
      const std::vector<double*> arrays = {x.data(), y.data(), z.data()};

      // contribute: one per region with atoms, so the owner's total is the number of contributing regions
      for(int cell = 0; cell < near_cells_end; cell++){
         const double value = (region_atoms[cell] > 0) ? 1.0 : 0.0;
         x[cell] = value;
         y[cell] = -value;
         z[cell] = 2.0 * value;
      }
      contribute_exchange.exchange(arrays, vmpi::sparse_exchange_t::combine_t::add);
      for(int cell = 0; cell < near_cells_end; cell++){
         if(owner[cell] != my_region || max_atoms[cell] == 0) continue;
         const double expected = contributors[cell];
         if(x[cell] != expected || y[cell] != -expected || z[cell] != 2.0 * expected) wrong++;
      }

      // serve: owners hold a value unique to the cell, which every reader must then hold
      for(int cell = 0; cell < near_cells_end; cell++){
         const double value = (owner[cell] == my_region) ? double(cell) + 1.0 : 0.0;
         x[cell] = value;
         y[cell] = -value;
         z[cell] = 2.0 * value;
      }
      serve_exchange.exchange(arrays, vmpi::sparse_exchange_t::combine_t::assign);
      for(int cell = 0; cell < near_cells_end; cell++){
         if(!region_reads[cell] || max_atoms[cell] == 0) continue;
         const double expected = double(cell) + 1.0;
         if(x[cell] != expected || y[cell] != -expected || z[cell] != 2.0 * expected) wrong++;
      }

      return wrong;

   }

   #endif

}

//------------------------------------------------------------------------------
// Function to set up the point-to-point exchange of cell moments between shared
// memory regions. The owner of each near cell is the region with most of its
// magnetic atoms (ties go to the lower region). Falls back to the all-reduce,
// with a log line, in the serial build, with a single region, and for
// micromagnetic or multiscale discretisations. Called at the end of
// initialisation, after the interaction lists are built; collective over all
// processes.
//------------------------------------------------------------------------------
void initialise_moment_exchange(){

   vutil::vtimer_t timer;
   timer.start();

   // cells up to the exchange level are near
   const int level = std::min(ha::moment_exchange_level, ha::num_levels - 1);
   near_cells_end = ha::cells_level_end_index[level];

   #ifdef MPICF

      int num_regions = (vmpi::my_shared_rank == 0) ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &num_regions, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      std::string fallback = "";
      if(num_regions == 1) fallback = "a single shared memory region";
      if(micromagnetic::discretisation_type != 0) fallback = "a micromagnetic or multiscale discretisation";
      if(!fallback.empty()){
         zlog << zTs() << "\tHierarchical moment exchange: point-to-point is not available for " << fallback
              << ", using all-reduce" << std::endl;
         ha::moment_exchange = ha::moment_exchange_t::all_reduce;
         return;
      }

      // magnetic atoms of this process in its contributing cells (its level 0 cells and their coarse cells, the
      // only cells where it has atoms), indexed like ha::contributing_cells, summed up the tree
      std::vector<int> atoms_in_cell(ha::contributing_cells.size(), 0);
      for(int atom = 0; atom < vmpi::num_local_atoms; ++atom){
         if(atoms::magnetic[atom]) atoms_in_cell[ha::atom_contributing_index[atom]]++;
      }
      for(size_t i = 0; i < ha::upward_pass_index.size(); i++){
         for(int c = ha::upward_children_start[i]; c < ha::upward_children_start[i+1]; c++){
            atoms_in_cell[ha::upward_pass_index[i]] += atoms_in_cell[ha::upward_children[c]];
         }
      }

      // the near cells among them, as (cell, atoms) pairs
      std::vector<int> atom_pairs;
      for(size_t i = 0; i < ha::contributing_cells.size(); i++){
         if(ha::contributing_cells[i] < near_cells_end && atoms_in_cell[i] > 0){
            atom_pairs.push_back(ha::contributing_cells[i]);
            atom_pairs.push_back(atoms_in_cell[i]);
         }
      }
      atoms_in_cell.clear();
      atoms_in_cell.shrink_to_fit();

      // cells read by this process: its interaction lists and its own level 0 cells, one bit per cell
      const int num_read_bytes = (ha::total_num_cells + 7) / 8;
      std::vector<unsigned char> reads(num_read_bytes, 0);
      auto mark_read = [&reads](const int cell){ reads[cell >> 3] |= static_cast<unsigned char>(1u << (cell & 7)); };
      for(const int cell : ha::interaction_list) mark_read(cell);
      for(int lc = 0; lc < dipole::cells_num_local_cells; lc++) mark_read(cells::cell_id_array[lc]);

      // the same for the region, on its leader (atoms in the near cells, reads in all cells)
      const bool leader = (vmpi::my_shared_rank == 0);
      std::vector<int> region_atoms(leader ? near_cells_end : 0, 0);
      std::vector<int> region_reads(leader ? ha::total_num_cells : 0, 0);
      {
         // atoms: gather the pairs of every process of the region and add them (integers, so in any order)
         const int num_pair_values = atom_pairs.size();
         std::vector<int> pair_counts(leader ? vmpi::num_shared_ranks : 0, 0);
         std::vector<int> pair_displacements(leader ? vmpi::num_shared_ranks : 0, 0);
         MPI_Gather(&num_pair_values, 1, MPI_INT, pair_counts.data(), 1, MPI_INT, 0, vmpi::shared_comm);
         int total_pair_values = 0;
         for(int r = 0; r < int(pair_counts.size()); r++){
            pair_displacements[r] = total_pair_values;
            total_pair_values += pair_counts[r];
         }
         std::vector<int> region_pairs(total_pair_values, 0);
         MPI_Gatherv(atom_pairs.data(), num_pair_values, MPI_INT, region_pairs.data(), pair_counts.data(),
                     pair_displacements.data(), MPI_INT, 0, vmpi::shared_comm);
         for(int i = 0; i < total_pair_values; i += 2) region_atoms[region_pairs[i]] += region_pairs[i+1];

         // reads: bitwise or of the bits of every process of the region
         std::vector<unsigned char> region_read_bits(leader ? num_read_bytes : 0, 0);
         MPI_Reduce(reads.data(), region_read_bits.data(), num_read_bytes, MPI_UNSIGNED_CHAR, MPI_BOR, 0, vmpi::shared_comm);
         for(int cell = 0; cell < int(region_reads.size()); cell++) region_reads[cell] = (region_read_bits[cell >> 3] >> (cell & 7)) & 1;
      }
      atom_pairs.clear();
      atom_pairs.shrink_to_fit();
      reads.clear();
      reads.shrink_to_fit();
      region_read_cells.clear();
      for(int cell = 0; cell < int(region_reads.size()); cell++) if(region_reads[cell]) region_read_cells.push_back(cell);

      // gather of the contributing cells of every process of the region to its leader
      const int num_indices = ha::contributing_cells.size();
      index_counts.assign(leader ? vmpi::num_shared_ranks : 0, 0);
      index_displacements.assign(leader ? vmpi::num_shared_ranks : 0, 0);
      value_counts.assign(leader ? vmpi::num_shared_ranks : 0, 0);
      value_displacements.assign(leader ? vmpi::num_shared_ranks : 0, 0);
      MPI_Gather(&num_indices, 1, MPI_INT, index_counts.data(), 1, MPI_INT, 0, vmpi::shared_comm);
      int total_indices = 0;
      for(int r = 0; r < int(index_counts.size()); r++){
         index_displacements[r] = total_indices;
         value_counts[r] = 3 * index_counts[r];
         value_displacements[r] = 3 * total_indices;
         total_indices += index_counts[r];
      }
      region_indices.assign(total_indices, 0);
      region_values.assign(3 * uint64_t(total_indices), 0.0);
      MPI_Gatherv(ha::contributing_cells.data(), num_indices, MPI_INT, region_indices.data(), index_counts.data(),
                  index_displacements.data(), MPI_INT, 0, vmpi::shared_comm);

      uint64_t wrong = 0;
      if(leader){

         int my_region = 0;
         MPI_Comm_rank(vmpi::leader_comm, &my_region);

         // owner of each near cell, its number of atoms, and the number of regions with atoms in it
         std::vector<int> owner(near_cells_end);
         std::vector<int> max_atoms(near_cells_end);
         std::vector<int> contributors(near_cells_end);
         {
            struct owner_t{ int atoms; int region; };
            std::vector<owner_t> pairs(near_cells_end);
            for(int cell = 0; cell < near_cells_end; cell++){
               pairs[cell] = {region_atoms[cell], my_region};
               contributors[cell] = (region_atoms[cell] > 0);
            }
            MPI_Allreduce(MPI_IN_PLACE, pairs.data(), near_cells_end, MPI_2INT, MPI_MAXLOC, vmpi::leader_comm);
            MPI_Allreduce(MPI_IN_PLACE, contributors.data(), near_cells_end, MPI_INT, MPI_SUM, vmpi::leader_comm);
            for(int cell = 0; cell < near_cells_end; cell++){
               owner[cell] = pairs[cell].region;
               max_atoms[cell] = pairs[cell].atoms;
            }
         }

         // cells sent to each owner, and cells received from each owner; cells with no atoms stay zero
         std::vector< std::vector<int> > to_owner(num_regions);
         std::vector< std::vector<int> > from_owner(num_regions);
         for(int cell = 0; cell < near_cells_end; cell++){
            if(owner[cell] == my_region || max_atoms[cell] == 0) continue;
            if(region_atoms[cell] > 0) to_owner[owner[cell]].push_back(cell);
            if(region_reads[cell]) from_owner[owner[cell]].push_back(cell);
         }
         contribute_exchange.initialise(vmpi::leader_comm, to_owner, vmpi::sparse_exchange_t::lists_t::sent, 3);
         serve_exchange.initialise(vmpi::leader_comm, from_owner, vmpi::sparse_exchange_t::lists_t::received, 3);

         wrong = check_exchange(region_atoms, region_reads, owner, max_atoms, contributors, my_region);

      }

      MPI_Allreduce(MPI_IN_PLACE, &wrong, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
      if(wrong > 0){
         err::fatal("Programmer error: the point-to-point exchange of hierarchical cell moments gives " +
                    std::to_string(wrong) + " wrong test values at initialisation");
      }

      // log the size of each stage, min/mean/max over region leaders
      if(leader){
         const vmpi::sparse_exchange_t& c = contribute_exchange;
         const vmpi::sparse_exchange_t& s = serve_exchange;
         const int far_cells = ha::total_num_cells - near_cells_end;
         zlog << zTs() << "\tHierarchical moment exchange: point-to-point between " << num_regions << " regions for levels 0-"
              << level << " (" << near_cells_end << " cells), all-reduce of " << far_cells << " cells ("
              << double(far_cells) * bytes_per_cell / 1.0e6 << " MB) above; test values checked" << std::endl;
         zlog << zTs() << "\tHierarchical moment exchange per region leader (min/mean/max):"
              << " contribute_peers=" << spread_over(c.num_send_peers() + c.num_recv_peers(), vmpi::leader_comm)
              << " contribute_mb=" << spread_over(double(c.num_send_indices() + c.num_recv_indices()) * bytes_per_cell / 1.0e6, vmpi::leader_comm)
              << " serve_peers=" << spread_over(s.num_send_peers() + s.num_recv_peers(), vmpi::leader_comm)
              << " serve_mb=" << spread_over(double(s.num_send_indices() + s.num_recv_indices()) * bytes_per_cell / 1.0e6, vmpi::leader_comm)
              << std::endl;
      }

   #else

      zlog << zTs() << "\tHierarchical moment exchange: point-to-point is not available in the serial build, using all-reduce" << std::endl;
      ha::moment_exchange = ha::moment_exchange_t::all_reduce;
      return;

   #endif

   timer.stop();
   zlog << zTs() << "\tHierarchical moment exchange set up in " << timer.elapsed_time() << " s" << std::endl;

   return;

}

//------------------------------------------------------------------------------
// Function to sum the contributions of all processes (ha::mag_reduction_values,
// x then y then z for each of ha::contributing_cells) into the shared cell
// moments, for the cells each region reads. Replaces ha::mag_reduction.sum when
// hierarchical:moment-exchange = point-to-point. The first call also runs the
// all-reduce beforehand and logs the largest difference over the cells read.
// Collective over all processes.
//
// Every process calls this only after it has finished reading the moments of
// the previous update, so once the gather completes the leader may overwrite them.
//------------------------------------------------------------------------------
void exchange_moments(){

   #ifdef MPICF

      vmpi::shared_array_t<double>* totals[3] = {&ha::mag_array_x, &ha::mag_array_y, &ha::mag_array_z};
      const int num_indices = ha::contributing_cells.size();
      const bool leader = (vmpi::my_shared_rank == 0);

      // once: the all-reduce first, as the reference for the check at the end. The moments are then zeroed, so
      // that a cell the exchange fails to fill shows as a difference (cells with no atoms are zero in both).
      const bool check = !cross_checked;
      std::vector<double> reference;
      if(check){
         cross_checked = true;
         ha::mag_reduction.sum(ha::mag_reduction_values, {&ha::mag_array_x, &ha::mag_array_y, &ha::mag_array_z});
         if(leader){
            reference.reserve(3 * region_read_cells.size());
            for(vmpi::shared_array_t<double>* total : totals){
               for(const int cell : region_read_cells) reference.push_back((*total)[cell]);
               std::fill(total->data(), total->data() + total->size(), 0.0);
            }
         }
      }

      // complete this process's reads of the previous moments before the leader writes them
      for(vmpi::shared_array_t<double>* total : totals) total->memory_barrier();

      double time = MPI_Wtime();
      MPI_Gatherv(ha::mag_reduction_values.data(), 3 * num_indices, MPI_DOUBLE, region_values.data(), value_counts.data(),
                  value_displacements.data(), MPI_DOUBLE, 0, vmpi::shared_comm);
      exchange_time[gather] = MPI_Wtime() - time;
      for(int phase = add; phase < barrier; phase++) exchange_time[phase] = 0.0;

      if(leader){

         // the region's partial moments: zero the cells it contributes to and the far cells, then add
         time = MPI_Wtime();
         for(int c = 0; c < 3; c++){
            double* total = totals[c]->data();
            for(const int cell : region_indices) total[cell] = 0.0;
            std::fill(total + near_cells_end, total + ha::total_num_cells, 0.0);
            for(int r = 0; r < int(index_counts.size()); r++){
               const int n = index_counts[r];
               const int* index = &region_indices[index_displacements[r]];
               const double* value = &region_values[value_displacements[r] + uint64_t(c) * n];
               for(int i = 0; i < n; i++) total[index[i]] += value[i];
            }
         }
         exchange_time[add] = MPI_Wtime() - time;

         const std::vector<double*> arrays = {totals[0]->data(), totals[1]->data(), totals[2]->data()};

         // near cells: partial moments to their owners, who add them in order of region
         time = MPI_Wtime();
         contribute_exchange.exchange(arrays, vmpi::sparse_exchange_t::combine_t::add);
         exchange_time[contribute] = MPI_Wtime() - time;

         // far cells: summed over all regions
         time = MPI_Wtime();
         const int far_cells = ha::total_num_cells - near_cells_end;
         if(far_cells > 0){
            for(double* total : arrays){
               MPI_Allreduce(MPI_IN_PLACE, total + near_cells_end, far_cells, MPI_DOUBLE, MPI_SUM, vmpi::leader_comm);
            }
         }
         exchange_time[far] = MPI_Wtime() - time;

         // near cells: totals from their owners to the regions that read them
         time = MPI_Wtime();
         serve_exchange.exchange(arrays, vmpi::sparse_exchange_t::combine_t::assign);
         exchange_time[serve] = MPI_Wtime() - time;

      }

      // make the moments visible to every process in the region (one barrier for all arrays)
      time = MPI_Wtime();
      for(vmpi::shared_array_t<double>* total : totals) total->memory_barrier();
      MPI_Barrier(vmpi::shared_comm);
      for(vmpi::shared_array_t<double>* total : totals) total->memory_barrier();
      exchange_time[barrier] = MPI_Wtime() - time;

      // once: compare the cells each region reads with the all-reduce
      if(check){
         double difference = 0.0;
         double largest = 0.0;
         if(leader){
            uint64_t i = 0;
            for(vmpi::shared_array_t<double>* total : totals){
               for(const int cell : region_read_cells){
                  difference = std::max(difference, std::abs((*total)[cell] - reference[i]));
                  largest = std::max(largest, std::abs(reference[i]));
                  i++;
               }
            }
         }
         difference = vmpi::reduce_max(difference);
         largest = vmpi::reduce_max(largest);
         const double relative = largest > 0.0 ? difference / largest : difference;
         zlog << zTs() << "\tHierarchical moment exchange check: largest difference from all-reduce " << difference
              << ", largest moment " << largest << ", relative difference " << relative << std::endl;
         if(!(relative <= 1.0e-10)){ // also catches NaN
            err::warning("Warning - the point-to-point exchange of hierarchical cell moments differs from the all-reduce by " +
                         std::to_string(relative) + " of the largest moment.");
         }
      }

   #endif

   return;

}

//------------------------------------------------------------------------------
// Function returning the time of each phase of the last exchange for the
// DIPOLE_SCALING line, with the keys of the all-reduce where they match (the
// far all-reduce as reduce_allreduce). The max over processes includes the
// wait for the slowest; the min, over region leaders where only they take part,
// is the cost without waiting. Collective over all processes.
//------------------------------------------------------------------------------
std::string moment_exchange_phase_times(){

   const bool leader = vmpi::my_shared_rank == 0;
   const double not_leader = 1.0e30; // excluded from the minimum
   auto leader_min = [&](const double value){ return vmpi::reduce_min(leader ? value : not_leader); };

   std::ostringstream times;
   times << " reduce_gather_s="         << vmpi::reduce_max(exchange_time[gather])
         << " reduce_gather_min_s="     << vmpi::reduce_min(exchange_time[gather])
         << " reduce_add_s="            << vmpi::reduce_max(exchange_time[add])
         << " reduce_contribute_s="     << vmpi::reduce_max(exchange_time[contribute])
         << " reduce_contribute_min_s=" << leader_min(exchange_time[contribute])
         << " reduce_allreduce_s="      << vmpi::reduce_max(exchange_time[far])
         << " reduce_allreduce_min_s="  << leader_min(exchange_time[far])
         << " reduce_serve_s="          << vmpi::reduce_max(exchange_time[serve])
         << " reduce_serve_min_s="      << leader_min(exchange_time[serve])
         << " reduce_barrier_s="        << vmpi::reduce_max(exchange_time[barrier])
         << " reduce_barrier_min_s="    << vmpi::reduce_min(exchange_time[barrier]);

   return times.str();

}

} // end of internal namespace
} // end of hierarchical namespace
