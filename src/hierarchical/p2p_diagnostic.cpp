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
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "cells.hpp"
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

   // minimum, mean and maximum of a value over processes as min/mean/max
   std::string spread(const double value){
      const double min = vmpi::reduce_min(value);
      const double mean = vmpi::reduce_sum(value) / double(vmpi::num_processors);
      const double max = vmpi::reduce_max(value);
      std::ostringstream text;
      text << min << "/" << mean << "/" << max;
      return text.str();
   }

   // number of nonzero entries in a list with one entry per process
   double count_nonzero(const std::vector<int>& per_process){
      return double(per_process.size() - std::count(per_process.begin(), per_process.end(), 0));
   }

   // sum of a list with one entry per process
   double sum_of(const std::vector<int>& per_process){
      double total = 0.0;
      for(const int value : per_process) total += value;
      return total;
   }

   // sum over processes of the entry that each process holds for this process
   int sum_entries_for_this_process(const std::vector<int>& per_process){
      int total = 0;
      #ifdef MPICF
         MPI_Reduce_scatter_block(per_process.data(), &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      #endif
      return total;
   }

}

//------------------------------------------------------------------------------
// Function to write to the log what a point-to-point exchange of cell moments
// would cost, compared with the reduction of all cells. Nothing here changes the
// calculation or is kept. Each cell is owned by the process with most of its
// atoms. An update would then send the partial moment of each contributed cell
// to its owner, and the total of each cell in an interaction list from its owner.
// One line is written for each level L: levels up to L use the exchange, and
// higher levels stay in an all-reduce. Cumulative values are min/mean/max over
// processes, and cells of other processes only are counted. Called during
// initialisation, collective over all processes.
//------------------------------------------------------------------------------
void write_p2p_diagnostic(){

   vutil::vtimer_t timer;
   timer.start();

   // magnetic atoms of this process in each cell at all levels
   std::vector<int> atoms_in_cell(ha::total_num_cells, 0);
   if(micromagnetic::discretisation_type != 1){
      for(int atom = 0; atom < vmpi::num_local_atoms; ++atom){
         if(atoms::magnetic[atom]) atoms_in_cell[cells::atom_cell_id_array[atom]]++;
      }
   }
   else{
      for(int lc = 0; lc < micromagnetic::number_of_micromagnetic_cells; lc++) atoms_in_cell[micromagnetic::list_of_micromagnetic_cells[lc]] = 1;
   }
   for(const int cell : ha::upward_pass_cells){
      for(int i = ha::cells_in_cells_start_index[cell]; i < ha::cells_in_cells_end_index[cell]; i++){
         atoms_in_cell[cell] += atoms_in_cell[ha::cells_in_cells[i]];
      }
   }

   // owner of each cell, and the largest number of atoms any process has in it (zero for an empty cell)
   struct owner_t{ int atoms; int rank; };
   std::vector<owner_t> owner(ha::total_num_cells);
   for(int cell = 0; cell < ha::total_num_cells; cell++) owner[cell] = {atoms_in_cell[cell], vmpi::my_rank};
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, owner.data(), ha::total_num_cells, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);
   #endif

   // cells read in the field loop, in order
   std::vector<int> read_list(ha::interaction_list);
   std::sort(read_list.begin(), read_list.end());
   read_list.erase(std::unique(read_list.begin(), read_list.end()), read_list.end());

   // cells with atoms above each level, which would stay in the all-reduce
   std::vector<uint64_t> cells_above_level(ha::num_levels, 0);
   uint64_t cells_with_atoms = 0;
   for(int level = ha::num_levels - 1; level >= 0; level--){
      uint64_t cells_in_level = 0;
      for(int cell = ha::cells_level_start_index[level]; cell < ha::cells_level_end_index[level]; cell++) if(owner[cell].atoms > 0) cells_in_level++;
      cells_above_level[level] = cells_with_atoms;
      cells_with_atoms += cells_in_level;
   }

   const double bytes_per_cell = 3.0 * sizeof(double);
   zlog << zTs() << "\tHierarchical P2P diagnostic: " << vmpi::num_processors << " processes, " << cells_with_atoms
        << " cells with atoms, all-reduce of all cells " << double(cells_with_atoms) * bytes_per_cell / 1.0e6 << " MB" << std::endl;

   // cells of levels up to L exchanged with each process, added to level by level
   std::vector<int> cells_read_from(vmpi::num_processors, 0);
   std::vector<int> cells_contributed_to(vmpi::num_processors, 0);
   size_t next_read = 0;
   size_t next_contributed = 0;

   for(int level = 0; level < ha::num_levels; level++){

      const int level_end = ha::cells_level_end_index[level];
      for(; next_read < read_list.size() && read_list[next_read] < level_end; next_read++){
         const owner_t& cell_owner = owner[read_list[next_read]];
         if(cell_owner.atoms > 0 && cell_owner.rank != vmpi::my_rank) cells_read_from[cell_owner.rank]++;
      }
      for(; next_contributed < ha::contributing_cells.size() && ha::contributing_cells[next_contributed] < level_end; next_contributed++){
         const int cell = ha::contributing_cells[next_contributed];
         if(atoms_in_cell[cell] > 0 && owner[cell].rank != vmpi::my_rank) cells_contributed_to[owner[cell].rank]++;
      }

      // processes this process hears from and sends to in each stage
      std::vector<int> reads_from_rank(vmpi::num_processors, 0);
      std::vector<int> contributes_to_rank(vmpi::num_processors, 0);
      for(int rank = 0; rank < vmpi::num_processors; rank++){
         reads_from_rank[rank] = (cells_read_from[rank] > 0);
         contributes_to_rank[rank] = (cells_contributed_to[rank] > 0);
      }

      const double read_cells = sum_of(cells_read_from);
      const double read_sources = count_nonzero(cells_read_from);
      const double served_cells = sum_entries_for_this_process(cells_read_from);
      const double served_readers = sum_entries_for_this_process(reads_from_rank);
      const double contributed_cells = sum_of(cells_contributed_to);
      const double contributed_owners = count_nonzero(cells_contributed_to);
      const double received_cells = sum_entries_for_this_process(cells_contributed_to);
      const double received_contributors = sum_entries_for_this_process(contributes_to_rank);

      const double messages = read_sources + served_readers + contributed_owners + received_contributors;
      const double megabytes = (read_cells + served_cells + contributed_cells + received_cells) * bytes_per_cell / 1.0e6;

      zlog << zTs() << "DIPOLE_P2P level=" << level
           << " allreduce_cells=" << cells_above_level[level]
           << " read_cells=" << spread(read_cells)
           << " read_sources=" << spread(read_sources)
           << " serve_cells=" << spread(served_cells)
           << " serve_readers=" << spread(served_readers)
           << " contribute_cells=" << spread(contributed_cells)
           << " contribute_owners=" << spread(contributed_owners)
           << " receive_cells=" << spread(received_cells)
           << " receive_contributors=" << spread(received_contributors)
           << " messages=" << spread(messages)
           << " mb=" << spread(megabytes) << std::endl;

   }

   timer.stop();
   zlog << zTs() << "\tHierarchical P2P diagnostic took " << timer.elapsed_time() << " s" << std::endl;

   return;

}

} // end of internal namespace
} // end of hierarchical namespace
