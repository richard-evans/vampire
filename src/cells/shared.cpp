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
//   Cell arrays held once per shared memory region (cells::shared_*).
//
//   With cells::full_cell_arrays the cells module builds its full size arrays
//   on every process as before (an MPI_Allreduce over every cell) and copies
//   them into the shared arrays. Otherwise (hierarchical dipole solver with no
//   other module reading the cell arrays) each process sums only its own
//   cells, and the sums are reduced into the shared arrays by
//   vmpi::shared_sum_t, so no process holds an array over every cell.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "cells.hpp"
#include "material.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// cells module headers
#include "internal.hpp"

namespace cells{
namespace internal{

//------------------------------------------------------------------------------
// Builds cells::num_atoms_in_cell (this process), cells::cell_id_array, cells::atomic_volume and the shared arrays
// from the local cells of each process (collective). The sums over each cell's atoms are taken in atom order, as
// in the full size path, so cells held by one or two processes have the same values.
//------------------------------------------------------------------------------
void initialise_shared_arrays_compact(const int num_local_atoms,
                                      const std::vector<double>& atom_coords_x,
                                      const std::vector<double>& atom_coords_y,
                                      const std::vector<double>& atom_coords_z,
                                      const int total_atoms_non_filler,
                                      const double unit_cell_size_x,
                                      const double unit_cell_size_y,
                                      const double unit_cell_size_z){

   const int num_cells = cells::num_cells;

   // magnetic atoms of this process in each cell, and the cells holding them
   cells::num_atoms_in_cell.assign(num_cells, 0);
   int num_atoms_magnetic = 0;
   for(int atom = 0; atom < num_local_atoms; atom++){
      if(mp::material[atoms::type_array[atom]].non_magnetic == 0){
         cells::num_atoms_in_cell[cells::atom_cell_id_array[atom]]++;
         num_atoms_magnetic++;
      }
   }
   for(int cell = 0; cell < num_cells; cell++){
      if(cells::num_atoms_in_cell[cell] > 0) cells::cell_id_array.push_back(cell);
   }

   // sums over the magnetic atoms of each local cell: x mu, y mu, z mu, mu and the number of atoms
   const int n = cells::cell_id_array.size();
   std::vector<double> values(5 * n, 0.0);
   for(int atom = 0; atom < num_local_atoms; atom++){
      const int type = atoms::type_array[atom];
      if(mp::material[type].non_magnetic != 0) continue;
      const double mus = mp::material[type].mu_s_SI;
      const int cell = cells::atom_cell_id_array[atom];
      const int k = std::lower_bound(cells::cell_id_array.begin(), cells::cell_id_array.end(), cell) - cells::cell_id_array.begin();
      values[        k] += atom_coords_x[atom] * mus;
      values[    n + k] += atom_coords_y[atom] * mus;
      values[2 * n + k] += atom_coords_z[atom] * mus;
      values[3 * n + k] += mus;
      values[4 * n + k] += 1.0;
   }

   // sums over all processes, one copy per region
   std::vector< vmpi::shared_array_t<double> > totals(5);
   for(auto& total : totals) total.allocate(num_cells);
   vmpi::shared_sum_t sum;
   sum.initialise(cells::cell_id_array, 5);
   sum.sum(values, {&totals[0], &totals[1], &totals[2], &totals[3], &totals[4]});
   std::vector<double>().swap(values);

   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &num_atoms_magnetic, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif

   // magnetic atomic volume, as in the full size path
   const double factor_for_volume = double(total_atoms_non_filler)/double(num_atoms_magnetic);
   const double atomic_volume =  factor_for_volume * unit_cell_size_x*unit_cell_size_y*unit_cell_size_z/double(cells::num_atoms_in_unit_cell);
   cells::atomic_volume = atomic_volume;
   zlog << zTs() << "Calculated (magnetic) atomic volume: " << cells::atomic_volume << " cubic Angstroms" << std::endl;

   // magnetic centre and moment, volume and number of magnetic atoms of each cell
   cells::shared_pos_and_mom.allocate(4 * uint64_t(num_cells));
   cells::shared_volume.allocate(num_cells);
   cells::shared_num_atoms_global.allocate(num_cells);
   if(cells::shared_pos_and_mom.is_writer()){
      for(int cell = 0; cell < num_cells; cell++){
         const int count = int(totals[4][cell]);
         const double moment = totals[3][cell];
         cells::shared_num_atoms_global[cell] = count;
         cells::shared_pos_and_mom[4*cell+3] = moment;
         for(int i = 0; i < 3; i++) cells::shared_pos_and_mom[4*cell+i] = count > 0 ? totals[i][cell] / moment : totals[i][cell];
         cells::shared_volume[cell] = double(count) * atomic_volume;
      }
   }
   cells::shared_pos_and_mom.synchronise();
   cells::shared_volume.synchronise();
   cells::shared_num_atoms_global.synchronise();

   for(auto& total : totals) total.free();

   zlog << zTs() << "Cell centres, moments and volumes held once per shared memory region (" << vmpi::num_shared_ranks
        << " processes share a copy); full size cell arrays not built" << std::endl;

   return;

}

//------------------------------------------------------------------------------
// Copies the full size cell arrays into the shared arrays (collective)
//------------------------------------------------------------------------------
void copy_full_arrays_to_shared(){

   const int num_cells = cells::num_cells;
   cells::shared_pos_and_mom.allocate(4 * uint64_t(num_cells));
   cells::shared_volume.allocate(num_cells);
   cells::shared_num_atoms_global.allocate(num_cells);
   if(cells::shared_pos_and_mom.is_writer()){
      std::copy(cells::pos_and_mom_array.begin(), cells::pos_and_mom_array.begin() + 4 * num_cells, cells::shared_pos_and_mom.data());
      std::copy(cells::volume_array.begin(), cells::volume_array.begin() + num_cells, cells::shared_volume.data());
      std::copy(cells::num_atoms_in_cell_global.begin(), cells::num_atoms_in_cell_global.begin() + num_cells, cells::shared_num_atoms_global.data());
   }
   cells::shared_pos_and_mom.synchronise();
   cells::shared_volume.synchronise();
   cells::shared_num_atoms_global.synchronise();

   return;

}

} // end of internal namespace
} // end of cells namespace
