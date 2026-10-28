//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2018. All rights reserved.
//
//   Email: sarah.jenkins@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <cmath>
//#include <cstdlib>
#include <iostream>
#include <string>

// Vampire headers
#include "cells.hpp" // needed for cells::cell_id_array but to be removed
#include "dipole.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"
#include "atoms.hpp"
#include "hierarchical.hpp"
#include "micromagnetic.hpp"

// hierarchical module headers
#include "internal.hpp"
#include "../dipole/internal.hpp"

// alias internal hierarchical namespace for brevity
namespace ha = hierarchical::internal;

//----------------------------------------------------------------------------------
// Function to return the index 0=x, 1=y, 2=z of the maximum value of three numbers
//----------------------------------------------------------------------------------
int largest(int x, int y, int z){
   int largest = x;
   int N = 0;
   if (y > largest) {
      largest = y;
      N = 1;
   }
   if (z > largest) {
      largest = z;
      N = 2;
   }
   return N;
}

//------------------------------------------------------------------------------
// Function to return the minimum value of three numbers
//------------------------------------------------------------------------------
int min(int x, int y, int z){
   int min = x;
   if (y < min) min = y;
   if (z > min) min = z;
   return min;
}


namespace hierarchical{

   // Forward declaration of cells_in_cells function
   void build_cells_in_cells(const std::vector<int>& cell_grid_i,
                             const std::vector<int>& cell_grid_j,
                             const std::vector<int>& cell_grid_k,
                             const std::vector<int>& level_ncy,
                             const std::vector<int>& level_ncz);

//----------------------------------------------------------------------------
// Function to initialize hierarchical module
//----------------------------------------------------------------------------
//
// Cell magnetization data is stored in a single 1D array for all levels
//
// |          < Level 1 cells >            | < Level 2 cells > |  <L3>   |<L4>|
// | L1 | L1 | L1 | L1 | L1 | L1 | L1 | L1 | L2 | L2 | L2 | L2 | L3 | L3 | L4 |
//
// This allows a parallel reduction of moments and distributed accumulation
// of cell moments and a single all_reduce for each spatial dimension
//
// The interaction matrix is distributed over all processors for memory
// efficiency with the hierarchical data is replicated (but even for millions
// of cells is only a few MB per process)
//
//----------------------------------------------------------------------------
void initialize(const double system_dimensions_x,
                const double system_dimensions_y,
                const double system_dimensions_z,
                std::vector<double>& atom_coords_x, //atomic coordinates
                std::vector<double>& atom_coords_y,
                std::vector<double>& atom_coords_z,
                int num_atoms){

   // store segment timings for optional printing
   std::vector <double> segment_time(9,0.0);

   // instantiate timer
   vutil::vtimer_t timer;

   //---------------------------------------------------------------------------
   // Segment 1
   //---------------------------------------------------------------------------
   // Setup up temporary data for computing hierarchical data
   //---------------------------------------------------------------------------

   //  start timer
   timer.start();

   // create the shared memory communicators (logs the region chosen: NUMA region or node)
   vmpi::initialise_shared_memory();

   // calculate the number of cells in x, y and z rounding up to powers of 2 (RE-?)
   int A = ceil(std::log(2.0 * system_dimensions_x / cells::macro_cell_size_x ) / std::log(2.0));
   int B = ceil(std::log(2.0 * system_dimensions_y / cells::macro_cell_size_y ) / std::log(2.0));
   int C = ceil(std::log(2.0 * system_dimensions_z / cells::macro_cell_size_z ) / std::log(2.0));
   // std::cout << dipole::cutoff*cells::macro_cell_size <<std::endl;

   //calculate the number of hierarchical levels necessary for the
   // calculation based on largest requirement for x,y,z
   int N = largest(A, B, C);

   if (N == 0) ha::num_levels = A;
   if (N == 1) ha::num_levels = B;
   if (N == 2) ha::num_levels = C;

   // always have a minumum of 1 level
   if (ha::num_levels < 1) ha::num_levels = 1;
   zlog << zTs() << "\tNumber of hierarchical levels: "  << ha::num_levels << std::endl; //"\t" << system_dimensions_x << '\t' << cells::macro_cell_size_x << "\t" << C << "\t" << N << "\t" << std::endl;

   // resize data to store level data
   ha::cells_level_start_index.resize(ha::num_levels,0.0);
   ha::cells_level_end_index.resize(ha::num_levels,0.0);
   ha::interaction_range.resize(ha::num_levels,0.0);

   // cell data from the cells module, read only (initialise_atomistic_cell_data takes const references)
   const std::vector < std::vector < int > >& cells_index_atoms_array = cells::index_atoms_array;
   const std::vector < int >& cells_num_atoms_in_cell        = cells::num_atoms_in_cell;        // this process's atoms only
   const std::vector < int >& cells_num_atoms_in_cell_global = cells::num_atoms_in_cell_global;
   const std::vector < int >& cells_local_cell_array         = cells::local_cell_array;
   const std::vector < double >& cells_pos_and_mom_array     = cells::pos_and_mom_array;
   int cells_num_cells = cells::num_cells;
   int cells_num_local_cells = cells::num_local_cells;

   // stop the timer for segment 1 and save it
   timer.stop();
   segment_time[0] = timer.elapsed_time();

   //------------------------------------------------
   // end of data copying
   //------------------------------------------------

   //---------------------------------------------------------------------------
   // Segment 2
   //---------------------------------------------------------------------------
   // For parallel execution only - rationalise data on all CPUs
   //---------------------------------------------------------------------------

   // restart the timer
   timer.start();

      // all processors wait here before starting the timer
      vmpi::barrier();

      // get realspace cutoff distance
      const double real_cutoff = dipole::cutoff * cells::macro_cell_size;

      std::vector< std::vector<double> > atoms_in_cells_array; // 2D list of [cell][atom] for local cells needed for computing dipole tensor
      std::vector<int> list_of_atoms_with_cells; // list of cell IDs to enable parsing of atomistic data

      // distribute atomistic data to enable hierarchical dipole calculation
      dipole::internal::initialise_atomistic_cell_data(cells_num_cells,
                                                       cells_num_local_cells,
                                                       real_cutoff,                     // cutoff range for dipole tensor construction (Angstroms)
                                                       cells_num_atoms_in_cell,         // number of atoms in each cell (local CPU)
                                                       cells_local_cell_array,          // numerical list of cells containing atoms on local processor
                                                       cells_num_atoms_in_cell_global,  // number of atoms in each cell
                                                       cells_pos_and_mom_array,         // array of positions and cell moments
                                                       cells_index_atoms_array,         // 2D array of [cells][atomID]
                                                       atom_coords_x,                   // input arrays of atom coordinates
                                                       atom_coords_y,                   //
                                                       atom_coords_z,                   //
                                                       atoms::m_spin_array,             // input array of atom moments (Bohr magnetons)
                                                       list_of_atoms_with_cells,        // 2D list of [cell][atom] for local cells needed for computing dipole tensor (output)
                                                       atoms_in_cells_array             // list of cell IDs to enable parsing of atomistic data (output)
                                                      );

   // After transferring the data across cores, assign value cells_num_atoms_in_cell[] from cells_num_atoms_in_cell_global[]
   for(unsigned int i=0; i<cells_num_atoms_in_cell_global.size(); i++){
      //if(cells_num_atoms_in_cell_global[i]>0 && cells_num_atoms_in_cell[i]==0){
         dipole::internal::cells_num_atoms_in_cell[i] = cells_num_atoms_in_cell_global[i];
      //}
   }

   // stop the timer for segment 2 and save it
   timer.stop();
   segment_time[1] = timer.elapsed_time();

   //------------------------------------------------
   // end of parallel data rationalisation
   //------------------------------------------------

   //---------------------------------------------------------------------------
   // Segment 3
   //---------------------------------------------------------------------------
   // loop over all levels to calculate the positions and sizes of the cells in
   // the levels
   //---------------------------------------------------------------------------

   // restart the timer
   timer.start();

   // determine average cell size
   internal::av_cell_size =  (cells::macro_cell_size_x +  cells::macro_cell_size_y +  cells::macro_cell_size_z)/3.0;

   // incrementing 1D index counter separating cells from different levels
   int index = 0;

   // per-level grid dimensions, used to fill the cell positions and by Segment 4 to locate each cell's
   // parent by direct index arithmetic instead of an O(N^2) bounding-box search
   std::vector<int> level_ncx(ha::num_levels,0);
   std::vector<int> level_ncy(ha::num_levels,0);
   std::vector<int> level_ncz(ha::num_levels,0);

   //loop over all levels to calculate the number of cells and interaction range of each level
   for (int level = 0; level < ha::num_levels; level ++){

      // determine cell size at each level 2^L
      double cell_size_x = pow(2,level) * cells::macro_cell_size_x;
      double cell_size_y = pow(2,level) * cells::macro_cell_size_y;
      double cell_size_z = pow(2,level) * cells::macro_cell_size_z;

      // at the lowest level use the standard dipole-dipole cutoff range * 2 (*3?)
      if (level == 0) ha::interaction_range[0] = dipole::cutoff * internal::av_cell_size * 2.0 * 1.5;
      else  ha::interaction_range[level] = ha::interaction_range[level - 1]*2.0;

      // for detailed info output level cutoffs
      //if(ha::info)
      zlog << zTs() << "\tHierarchical level " << level << "\t range: " << ha::interaction_range[level] << "\t" << dipole::cutoff*cells::macro_cell_size_x<< "\t" << internal::av_cell_size  <<'\t' << cells::macro_cell_size_x << "\t" << cells::macro_cell_size_y << "\t" << cells::macro_cell_size_z << "\t" <<  std::endl;

      // For last level set range to infinity (RE-should be 1e9?)
      if (level == ha::num_levels - 1 ) ha::interaction_range[level] = internal::av_cell_size * dipole::cutoff * 10000;

      // Calculate number of macrocells at each level
      // determine number of cells in x and y (global)
      level_ncx[level] = static_cast<unsigned int>(ceil((system_dimensions_x+0.01)/cell_size_x));
      level_ncy[level] = static_cast<unsigned int>(ceil((system_dimensions_y+0.01)/cell_size_y));
      level_ncz[level] = static_cast<unsigned int>(ceil((system_dimensions_z+0.01)/cell_size_z));

      const int temp_num_cells = level_ncx[level]*level_ncy[level]*level_ncz[level];

      //set the start and end index for the level
      ha::cells_level_start_index[level] = index;
      index = index + temp_num_cells;
      ha::cells_level_end_index[level] = index;

   } //end of loop over levels

   // total number of cells is last value of the index
   ha::total_num_cells = index;

   // allocate the hierarchical cell data shared by all processes in a shared memory region (collective).
   // Every cell below the top level is in exactly one cell of the level above.
   const int num_cells_in_cells = (ha::num_levels > 1) ? ha::cells_level_end_index[ha::num_levels - 2] : 0;
   ha::cell_positions.allocate(3 * uint64_t(ha::total_num_cells));
   ha::cell_dimensions.allocate(3 * uint64_t(ha::total_num_cells));
   ha::cell_positions_mom.allocate(4 * uint64_t(ha::total_num_cells));
   ha::cells_in_cells.allocate(num_cells_in_cells);
   ha::cells_in_cells_start_index.allocate(ha::total_num_cells);
   ha::cells_in_cells_end_index.allocate(ha::total_num_cells);

   // per-cell grid indices (i,j,k) within their level, index-aligned with ha::cell_positions/ha::cell_dimensions,
   // used by Segment 4 below (on the process filling the shared arrays only)
   std::vector<int> cell_grid_i;
   std::vector<int> cell_grid_j;
   std::vector<int> cell_grid_k;

   //loop over all levels to calculate the positions and sizes of the cells in the levels (one process per region)
   if (ha::cell_positions.is_writer()){

      cell_grid_i.resize(ha::total_num_cells);
      cell_grid_j.resize(ha::total_num_cells);
      cell_grid_k.resize(ha::total_num_cells);

      for (int level = 0; level < ha::num_levels; level ++){

         // determine cell size at each level 2^L
         double cell_size_x = pow(2,level) * cells::macro_cell_size_x;
         double cell_size_y = pow(2,level) * cells::macro_cell_size_y;
         double cell_size_z = pow(2,level) * cells::macro_cell_size_z;

         const int ncx = level_ncx[level];
         const int ncy = level_ncy[level];
         const int ncz = level_ncz[level];

         int cell = ha::cells_level_start_index[level];
         double size_x,size_y,size_z;
         double x,y,z;
         // loop over all cells in level in x,y,z (RE-complicated!)
         for(int i=0;i<ncx;++i){
            if (i < ncx -1)  size_x = cell_size_x;
            else             size_x = system_dimensions_x - (ncx-1)*cell_size_x;
            x = i*cell_size_x + size_x/2.0;
            for(int j=0;j<ncy;++j){
               if (j < ncy -1) size_y = cell_size_y;
               else            size_y = system_dimensions_y - (ncy-1)*cell_size_y;
               y = j*cell_size_y + size_y/2.0;
               // store cell coordinates
               for(int k=0; k<ncz; ++k){
                  if (k < ncz -1) size_z = cell_size_z;
                  else            size_z = system_dimensions_z - (ncz-1)*cell_size_z;
                  z = k*cell_size_z + size_z/2.0;

                  // now store cell positions and sizes
                  ha::cell_positions[3*cell + 0] = x;
                  ha::cell_positions[3*cell + 1] = y;
                  ha::cell_positions[3*cell + 2] = z;
                  ha::cell_dimensions[3*cell + 0] = size_x;
                  ha::cell_dimensions[3*cell + 1] = size_y;
                  ha::cell_dimensions[3*cell + 2] = size_z;
                  cell_grid_i[cell] = i;
                  cell_grid_j[cell] = j;
                  cell_grid_k[cell] = k;
                  cell++;
               }
            }
         }
      } //end of loop over levels

   }

   // stop the timer for segment 3 and save it
   timer.stop();
   segment_time[2] = timer.elapsed_time();

   //------------------------------------------------
   // end of hierarchical cell position calculation
   //------------------------------------------------

   //---------------------------------------------------------------------------
   // Segment 4
   //---------------------------------------------------------------------------
   // loop over all levels to determine a list of cells in level L that reside
   // in a cell at L+1
   //---------------------------------------------------------------------------

   // restart the timer
   timer.start();

   if (ha::cells_in_cells.is_writer()) build_cells_in_cells(cell_grid_i, cell_grid_j, cell_grid_k, level_ncy, level_ncz);

   // stop the timer for segment 4 and save it
   timer.stop();
   segment_time[3] = timer.elapsed_time();

   //------------------------------------------------
   // end of determination of cells in cells
   //------------------------------------------------

   //---------------------------------------------------------------------------
   // Segment 5
   //---------------------------------------------------------------------------
   // populate zero level hierarchical cells data
   //---------------------------------------------------------------------------

   // restart the timer
   timer.start();

   // determine number of zero level cells
   ha::num_zero_level_cells = ha::cells_level_end_index[0];
   //std::cout << ha::num_zero_level_cells << std::endl;

   // resize data structure to store
   ha::interaction_list_start_index.resize(ha::num_zero_level_cells,false);
   ha::interaction_list_end_index.resize(ha::num_zero_level_cells,false);
   // cell magnetisations, shared by the processes of a region (collective), and this process's contributions
   ha::mag_array_x.allocate(ha::total_num_cells);
   ha::mag_array_y.allocate(ha::total_num_cells);
   ha::mag_array_z.allocate(ha::total_num_cells);
   ha::partial_mag_array_x.assign(ha::total_num_cells,0.0);
   ha::partial_mag_array_y.assign(ha::total_num_cells,0.0);
   ha::partial_mag_array_z.assign(ha::total_num_cells,0.0);

   // check that the number of zero levels cells is correct (must be the same for hierachical and cells structures)
   if (ha::num_zero_level_cells != cells_num_cells) std::cout << "ERROR ZERO CELLS != NUM CELLS" << std::endl;
   //    std::cout << ha::num_zero_level_cells << '\t' << cells_num_cells << '\t' << ha::cell_positions_mom.size()/4. << '\t' << cells_pos_and_mom_array.size()/4. <<"\t" <<  ha::num_atoms_in_cell.size() << '\t' << cells_num_atoms_in_cell.size() <<std::endl;

   // cell positions and moments are the same on all processes (reduced in the cells module)
   if (ha::cell_positions_mom.is_writer()){
      for (int cell = 0; cell < cells_num_cells; cell ++){

         ha::cell_positions_mom[4*cell+0] = cells_pos_and_mom_array[4*cell+0];
         ha::cell_positions_mom[4*cell+1] = cells_pos_and_mom_array[4*cell+1];
         ha::cell_positions_mom[4*cell+2] = cells_pos_and_mom_array[4*cell+2];
         ha::cell_positions_mom[4*cell+3] = cells_pos_and_mom_array[4*cell+3];

      }
   }

   // stop the timer for segment 5 and save it
   timer.stop();
   segment_time[4] = timer.elapsed_time();

   //------------------------------------------------------
   // end of zero level cell data population
   //------------------------------------------------------

   //---------------------------------------------------------------------------
   // Segment 6
   //---------------------------------------------------------------------------
   // populate other level hierarchical cells data
   //---------------------------------------------------------------------------

   // restart the timer
   timer.start();

   for (int level = 1; level < ha::num_levels && ha::cell_positions_mom.is_writer(); level ++ ){
      int start = ha::cells_level_start_index[level];
      int end   = ha::cells_level_end_index[level];
      //   std::cout << "A" << level << '\t' << start << "\t" << end <<std::endl;
      for (int cell = start; cell < end; cell++){
         int start_cell_in_cell = ha::cells_in_cells_start_index[cell];
         int end_cell_in_cell = ha::cells_in_cells_end_index[cell];
         // centre of moment of the subcells (empty subcells have zero moment and position), and total moment
         for (int interaction = start_cell_in_cell; interaction < end_cell_in_cell; interaction ++ ){
            int cell_in_cell = ha::cells_in_cells[interaction];
            const double mom = ha::cell_positions_mom[4*cell_in_cell+3];
            ha::cell_positions_mom[4*cell+0] += ha::cell_positions_mom[4*cell_in_cell+0] * mom;
            ha::cell_positions_mom[4*cell+1] += ha::cell_positions_mom[4*cell_in_cell+1] * mom;
            ha::cell_positions_mom[4*cell+2] += ha::cell_positions_mom[4*cell_in_cell+2] * mom;
            ha::cell_positions_mom[4*cell+3] += mom;

         }
         const double total_mom = ha::cell_positions_mom[4*cell+3];
         for (int i = 0; i < 3; i++){
            // cells with no moment use their geometric centre
            if (total_mom > 0.0) ha::cell_positions_mom[4*cell+i] = ha::cell_positions_mom[4*cell+i] / total_mom;
            else                 ha::cell_positions_mom[4*cell+i] = ha::cell_positions[3*cell+i];
         }
      }
   }

   // make the hierarchical cell data visible to all processes in the region
   ha::cell_positions.synchronise();
   ha::cell_dimensions.synchronise();
   ha::cell_positions_mom.synchronise();
   ha::cells_in_cells.synchronise();
   ha::cells_in_cells_start_index.synchronise();
   ha::cells_in_cells_end_index.synchronise();

   // list the coarse cells holding level 0 cells whose moments this process accumulates (its local cells, or its
   // micromagnetic cells), level by level, so that the upward pass skips cells where its partial sums are zero
   std::vector<bool> has_local_cells(ha::total_num_cells, false);
   for (int lc = 0; lc < cells_num_local_cells; lc++) has_local_cells[cells::cell_id_array[lc]] = true;
   for (int lc = 0; lc < micromagnetic::number_of_micromagnetic_cells; lc++) has_local_cells[micromagnetic::list_of_micromagnetic_cells[lc]] = true;
   ha::upward_pass_cells.clear();
   for (int level = 1; level < ha::num_levels; level ++ ){
      for (int cell = ha::cells_level_start_index[level]; cell < ha::cells_level_end_index[level]; cell++){
         for (int interaction = ha::cells_in_cells_start_index[cell]; interaction < ha::cells_in_cells_end_index[cell]; interaction ++ ){
            if (has_local_cells[ha::cells_in_cells[interaction]]) has_local_cells[cell] = true;
         }
         if (has_local_cells[cell]) ha::upward_pass_cells.push_back(cell);
      }
   }

   // cells this process contributes moments to (its level 0 cells and their coarse cells, in order), and the sum
   // of the contributions over processes. Every magnetic atom of this process must be in one of these cells.
   ha::contributing_cells.clear();
   for (int cell = 0; cell < ha::total_num_cells; cell++) if (has_local_cells[cell]) ha::contributing_cells.push_back(cell);
   if( micromagnetic::discretisation_type != 1 ){
      for(int atom = 0; atom < vmpi::num_local_atoms; ++atom){
         if( atoms::magnetic[atom] && !has_local_cells[cells::atom_cell_id_array[atom]] ){
            err::fatal("Programmer error: magnetic atom " + std::to_string(atom) + " is in cell " + std::to_string(cells::atom_cell_id_array[atom]) + ", which is not a local cell of the hierarchical dipole solver");
         }
      }
   }
   ha::mag_reduction_values.assign(3 * ha::contributing_cells.size(), 0.0);
   ha::mag_reduction.initialise(ha::contributing_cells, 3);

   // stop the timer for segment 6 and save it
   timer.stop();
   segment_time[5] = timer.elapsed_time();

   //------------------------------------------------------
   // end of other level cell data population
   //------------------------------------------------------

   //---------------------------------------------------------------------------
   // Segment 7
   //---------------------------------------------------------------------------
   // Calculate interaction list for all cells (looping over LOCAL cells only)
   // this interaction list is therefore distributed, saving memory
   //---------------------------------------------------------------------------

   // restart the timer
   timer.start();

   // level of each cell, and the roots of the hierarchy (cells that are not a subcell of another cell)
   std::vector<int> cell_level(ha::total_num_cells, 0);
   std::vector<bool> is_root(ha::total_num_cells, true);
   for (int level = 0; level < ha::num_levels; level ++ ){
      for (int cell = ha::cells_level_start_index[level]; cell < ha::cells_level_end_index[level]; cell++){
         cell_level[cell] = level;
         if (level > 0){
            for (int interaction = ha::cells_in_cells_start_index[cell]; interaction < ha::cells_in_cells_end_index[cell]; interaction ++ ){
               is_root[ha::cells_in_cells[interaction]] = false;
            }
         }
      }
   }
   std::vector<int> roots;
   for (int cell = 0; cell < ha::total_num_cells; cell++) if (is_root[cell]) roots.push_back(cell);

   // determine the interaction list for local cells by descending the hierarchy from the roots: a cell at level
   // L > 0 interacts as a whole if all its corners are beyond the level L interaction range, otherwise its subcells
   // are considered; level 0 cells always interact. The list is ordered by level (highest first) and by cell id
   // within a level.
   std::vector < std::vector <int> > level_interactions_lc(ha::num_levels);
   std::vector <int> stack;
   int interaction_num = 0;
   for (int lc = 0; lc < cells_num_local_cells; lc ++){

      ha::interaction_list_start_index[lc] = interaction_num;
      int cell_i = cells::cell_id_array[lc];
      double xi = ha::cell_positions[cell_i*3 + 0];
      double yi = ha::cell_positions[cell_i*3 + 1];
      double zi = ha::cell_positions[cell_i*3 + 2];

      for (int level = 0; level < ha::num_levels; level ++ ) level_interactions_lc[level].clear();
      stack = roots;
      while (!stack.empty()){

         const int cell_j = stack.back();
         stack.pop_back();
         const int level = cell_level[cell_j];

         if (level == 0){
            level_interactions_lc[0].push_back(cell_j);
            continue;
         }

         // count the corners of cell j beyond the interaction range (corners as in ha::calculate_corners)
         const double min_x = ha::cell_positions[cell_j*3 + 0] - ha::cell_dimensions[cell_j*3 + 0]/2.0;
         const double min_y = ha::cell_positions[cell_j*3 + 1] - ha::cell_dimensions[cell_j*3 + 1]/2.0;
         const double min_z = ha::cell_positions[cell_j*3 + 2] - ha::cell_dimensions[cell_j*3 + 2]/2.0;
         const double max_x = ha::cell_positions[cell_j*3 + 0] + ha::cell_dimensions[cell_j*3 + 0]/2.0;
         const double max_y = ha::cell_positions[cell_j*3 + 1] + ha::cell_dimensions[cell_j*3 + 1]/2.0;
         const double max_z = ha::cell_positions[cell_j*3 + 2] + ha::cell_dimensions[cell_j*3 + 2]/2.0;
         const double corner_x[2] = {min_x, max_x};
         const double corner_y[2] = {min_y, max_y};
         const double corner_z[2] = {min_z, max_z};

         int is_corner_outside_range = 0;
         for (int corner = 0; corner < 8; corner ++){
            const double cx = corner_x[(corner >> 2) & 1];
            const double cy = corner_y[(corner >> 1) & 1];
            const double cz = corner_z[corner & 1];
            double dx_2 = (xi - cx)*(xi - cx);
            double dy_2 = (yi - cy)*(yi - cy);
            double dz_2 = (zi - cz)*(zi - cz);
            double r = sqrt(dx_2 + dy_2 + dz_2);
            if (r > ha::interaction_range[level]) is_corner_outside_range ++;
         }

         if (is_corner_outside_range == 8) level_interactions_lc[level].push_back(cell_j);
         else {
            for (int interaction = ha::cells_in_cells_start_index[cell_j]; interaction < ha::cells_in_cells_end_index[cell_j]; interaction ++ ){
               stack.push_back(ha::cells_in_cells[interaction]);
            }
         }

      }

      for (int level = ha::num_levels - 1; level > -1; level -- ){
         std::sort(level_interactions_lc[level].begin(), level_interactions_lc[level].end());
         for (const int cell_j : level_interactions_lc[level]) ha::interaction_list.push_back(cell_j);
         interaction_num += level_interactions_lc[level].size();
      }
      ha::interaction_list_end_index[lc] = interaction_num;

   }

   // output the number of interactions with cells of each level (local cells on this process)
   std::vector<uint64_t> level_interactions(ha::num_levels, 0);
   for (size_t i = 0; i < ha::interaction_list.size(); i++){
      for (int level = 0; level < ha::num_levels; level++){
         if (ha::interaction_list[i] >= ha::cells_level_start_index[level] && ha::interaction_list[i] < ha::cells_level_end_index[level]) level_interactions[level]++;
      }
   }
   for (int level = 0; level < ha::num_levels; level++){
      zlog << zTs() << "\tHierarchical level " << level << " interactions (local cells): " << level_interactions[level] << std::endl;
   }

   // spread over processes of the field loop work: local cells, interactions and distinct cells read from ha::mag_array
   {
      std::vector<int> distinct_cells(ha::interaction_list);
      std::sort(distinct_cells.begin(), distinct_cells.end());
      distinct_cells.erase(std::unique(distinct_cells.begin(), distinct_cells.end()), distinct_cells.end());
      const double counts[3] = { double(dipole::cells_num_local_cells), double(ha::interaction_list.size()), double(distinct_cells.size()) };
      const char* names[3] = { "local cells", "interactions", "distinct interacting cells" };
      for(int i = 0; i < 3; i++){
         const double min = vmpi::reduce_min(counts[i]);
         const double mean = vmpi::reduce_sum(counts[i]) / double(vmpi::num_processors);
         const double max = vmpi::reduce_max(counts[i]);
         zlog << zTs() << "\tHierarchical " << names[i] << " per process (min, mean, max): " << min << ", " << mean << ", " << max << std::endl;
      }
   }

   // stop the timer for segment 7 and save it
   timer.stop();
   segment_time[6] = timer.elapsed_time();

   //------------------------------------------------------
   // end of calculation of the interaction list
   //------------------------------------------------------

   //---------------------------------------------------------------------
   // Check memory requirements and print to screen
   //---------------------------------------------------------------------

   const double total_memory = vmpi::reduce_sum( double(ha::interaction_list.size() ) * 6.0 * 8.0 / 1.0e6 ); // in Megabytes

   zlog << zTs() << "\tTotal memory for hierarchical dipole calculation (all CPUs): " << total_memory << " MB of RAM" << std::endl;
   std::cout << "Total memory for hierarchical dipole calculation (all CPUs): " << total_memory << " MB of RAM" << std::endl;

   zlog << zTs() << "\tNumber of local cells for hierarchical dipole calculation = " << cells_num_local_cells << std::endl;
   zlog << zTs() << "\tNumber of total cells for hierarchical dipole calculation = " << ha::total_num_cells << std::endl;

   // memory held by each process for the arrays replicated on every process (O(total cells)) and for the
   // distributed interaction arrays, the candidates for shared memory buffers (dp_plan step 4)
   auto mb = [](const auto& v){ return double(v.size() * sizeof(v[0])) / 1.0e6; };
   const double shared_ranks = double(vmpi::num_shared_ranks); // shared arrays counted as a share of the region's copy
   const double geometry_mb = (mb(ha::cell_positions) + mb(ha::cell_positions_mom) + mb(ha::cell_dimensions)) / shared_ranks;
   const double topology_mb = (mb(ha::cells_in_cells) + mb(ha::cells_in_cells_start_index) + mb(ha::cells_in_cells_end_index)) / shared_ranks;
   const double moment_mb   = (mb(ha::mag_array_x) + mb(ha::mag_array_y) + mb(ha::mag_array_z)) / shared_ranks
                            + mb(ha::partial_mag_array_x) + mb(ha::partial_mag_array_y) + mb(ha::partial_mag_array_z);
   const double cells_mb    = mb(cells::pos_and_mom_array) + mb(cells::pos_array) + mb(cells::num_atoms_in_cell) + mb(cells::num_atoms_in_cell_global)
                            + mb(cells::mag_array_x) + mb(cells::mag_array_y) + mb(cells::mag_array_z) + mb(cells::volume_array)
                            + 4.0 * double(cells::index_atoms_array.size() * sizeof(std::vector<int>)) / 1.0e6; // empty inner vectors of index_atoms_array and atom_in_cell_coords_array_*
   const double dipole_mb   = mb(dipole::cells_field_array_x) + mb(dipole::cells_field_array_y) + mb(dipole::cells_field_array_z)
                            + mb(dipole::cells_volume_array) + mb(dipole::internal::cells_num_atoms_in_cell)
                            + mb(dipole::internal::cells_pos_and_mom_array);
   const double interaction_mb = mb(ha::interaction_list) + mb(ha::interaction_list_start_index) + mb(ha::interaction_list_end_index)
                               + 6.0 * double(ha::interaction_list.size() * sizeof(double)) / 1.0e6; // tensors are resized below
   zlog << zTs() << "\tHierarchical memory per process (MB): replicated geometry " << geometry_mb << ", topology " << topology_mb
        << ", moments " << moment_mb << ", cells module " << cells_mb << ", dipole module " << dipole_mb
        << "; distributed interactions " << vmpi::reduce_max(interaction_mb) << " (largest process)" << std::endl;

   //---------------------------------------------------------------------------
   // Segment 8
   //---------------------------------------------------------------------------
   // Resize arrays to store hierarchical interaction tensors
   //---------------------------------------------------------------------------

   // restart the timer
   timer.start();

   //  std::cout << interaction_num << "\t" << cells_num_local_cells << '\t' << cells_num_local_cells*cells_num_cells << std::endl;
   ha::rij_tensor_xx.resize(interaction_num,0.0);
   ha::rij_tensor_xy.resize(interaction_num,0.0);
   ha::rij_tensor_xz.resize(interaction_num,0.0);
   ha::rij_tensor_yy.resize(interaction_num,0.0);
   ha::rij_tensor_yz.resize(interaction_num,0.0);
   ha::rij_tensor_zz.resize(interaction_num,0.0);

   // resize B-field cells array
   dipole::cells_field_array_x.resize(cells_num_cells,0.0);
   dipole::cells_field_array_y.resize(cells_num_cells,0.0);
   dipole::cells_field_array_z.resize(cells_num_cells,0.0);

   // stop the timer for segment 8 and save it
   timer.stop();
   segment_time[7] = timer.elapsed_time();


   //---------------------------------------------------------------------------
   // Segment 9
   //---------------------------------------------------------------------------
   // Calculate dipole tensor between cells
   // for nearby cells this is done atomistically
   // for cells further away the dipole-dipole approximation is used
   //---------------------------------------------------------------------------

   // restart the timer
   timer.start();

   //std::cout << cells_num_local_cells << std::endl;

   // loop over all local cells which need to be calculated
   for(int lc=0; lc<cells_num_local_cells; lc++){

      // get the global cell ID
      int cell_i = cells::cell_id_array[lc];

      // only calculate interaction for local cells with atoms
      if (cells_num_atoms_in_cell[cell_i] != 0){

         // for each cell loop over cel-cell interactions for all levels, getting frirst and last interactions
         const int start = ha::interaction_list_start_index[lc];
         const int end = ha::interaction_list_end_index[lc];
         // std::cerr << cell_i << '\t' << "interaction" << "\t" << end - start << "\t" << cells_num_cells << "\t" << cells_num_atoms_in_cell[cell_i] << std::endl;

         // now loop over all interactions
         for(int j = start; j < end; j++){

            // get global cell_j ID
            int cell_j = ha::interaction_list[j];

            // coarse cells (level >= 1) always use the bare dipole tensor: they are accepted only beyond the level
            // interaction range, which exceeds the atomistic cutoff
            const bool coarse = cell_j >= ha::num_zero_level_cells;

            // only calculate interaction if there are atoms in a remote level 0 cell
            if ( coarse || cells_num_atoms_in_cell_global[cell_j] != 0 ){

               //--------------------------------------------------------------
               // Calculation of inter part of dipolar tensor
               //--------------------------------------------------------------
               // check that remote cell is actually magnetic
               //if ( cells_num_atoms_in_cell_global[cell_j] > 0 ){

                  // check if cell is not the same
                  if( coarse || cell_i != cell_j ){
                     // calculate inter term of dipolar tensor
                     internal::calc_inter(cell_i,
                                          cell_j,
                                          j,
                                          real_cutoff,
                                          cells_num_atoms_in_cell_global,
                                          cells_pos_and_mom_array,
                                          list_of_atoms_with_cells,
                                          atoms_in_cells_array);

                  } // End of Inter part

                  //--------------------------------------------------------------
                  // Calculation of intra part of dipolar tensor
                  //--------------------------------------------------------------
                  else if( cell_i == cell_j ){

                     //Compute inter component of dipolar tensor
                     internal::calc_intra(cell_i,
                                          cell_j,
                                          j,
                                          cells_num_atoms_in_cell_global,
                                          list_of_atoms_with_cells,
                                          atoms_in_cells_array);

                  }
                  // End of Intra part
               //}
               //    // check for close to zero value tensors and round down to zero
               //
               // //   if ( cells_num_atoms_in_cell[cell_j] == 0 && cell_i == 0) std::cout << "HHHHM" << cell_j <<  "\t" << cells_num_atoms_in_cell[cell_j] << std::endl;
               // // if cell is not the same, compute the inter tensor else compute the intra tensor
               // if(cell_i!=cell_j) {
               //    ha::calc_inter(cell_i,cell_j,j, cells_atom_in_cell_coords_array_x, cells_atom_in_cell_coords_array_y, cells_atom_in_cell_coords_array_z);
               // }
               // else if (cell_i == cell_j){
               //    ha::calc_intra(cell_i,cell_j,j, cells_atom_in_cell_coords_array_x, cells_atom_in_cell_coords_array_y, cells_atom_in_cell_coords_array_z);
               // }

               // zero excessively small components (avoids round off errors comparing different calculation methods)
               if (ha::rij_tensor_xx[j]*ha::rij_tensor_xx[j] < 1e-15) ha::rij_tensor_xx[j] =0;
               if (ha::rij_tensor_xy[j]*ha::rij_tensor_xy[j] < 1e-15) ha::rij_tensor_xy[j] =0;
               if (ha::rij_tensor_xz[j]*ha::rij_tensor_xz[j] < 1e-15) ha::rij_tensor_xz[j] =0;
               if (ha::rij_tensor_yy[j]*ha::rij_tensor_yy[j] < 1e-15) ha::rij_tensor_yy[j] =0;
               if (ha::rij_tensor_yz[j]*ha::rij_tensor_yz[j] < 1e-15) ha::rij_tensor_yz[j] =0;
               if (ha::rij_tensor_zz[j]*ha::rij_tensor_zz[j] < 1e-15) ha::rij_tensor_zz[j] =0;
               //            std::cout  << cell_i << '\t' << cell_j << "\t" << ha::rij_tensor_xx[j] << '\t' << ha::rij_tensor_xy[j] << '\t' << ha::rij_tensor_xz[j] << '\t' << ha::rij_tensor_yy[j] << '\t' << ha::rij_tensor_yz[j] <<'\t' << ha::rij_tensor_zz[j] << std::endl;
               //             //std::cout << cells_num_atoms_in_cell[cell_i] << '\t' << cells_num_atoms_in_cell[cell_j] << '\t' <<  cell_i << '\t' << cell_j << "\t" <<  ha::rij_tensor_xx[j] << '\t' << ha::rij_tensor_xy[j] << '\t' << ha::rij_tensor_xz[j] << '\t' << ha::rij_tensor_yy[j] << '\t' << ha::rij_tensor_yz[j] <<'\t' << ha::rij_tensor_zz[j] << std::endl;
            }
         }
      }
   }

   // stop the timer for segment 9 and save it
   timer.stop();
   segment_time[8] = timer.elapsed_time();

   //---------------------------------------------------------------------------
   // Output detailed timings for the hierarchical initialisation
   //---------------------------------------------------------------------------

   // cell geometry is needed only to build the interactions
   ha::cell_positions.free();
   ha::cell_dimensions.free();
   ha::cell_positions_mom.free();

   // print out detailed hierarchical timings
   zlog << zTs() << "\tSegment times for hierarchical initialisation:" << std::endl;
   for(int i=0; i<9; i++) zlog << zTs() << "\t\t Segment time " << i+1 << ": " << segment_time[i] << " s" << std::endl;

   // Woohoo we made it!
   return;

}

//----------------------------------------------------------------------------
// Segment 4: determine, for every cell, the list of cells in the level below
// that it contains, with O(N) index arithmetic (the original O(N_level *
// N_sublevel) bounding-box search, build_cells_in_cells_original, was removed
// with the move to shared memory arrays; see git history).
//
// Valid because macrocell size doubles exactly each level (Segment 3), so a
// sublevel cell's grid index (i,j,k) always maps to parent grid index
// (i/2, j/2, k/2) -- including at the truncated last cell in each
// dimension, whose position (only its size) is unaffected by the
// truncation.
//
// Implemented as a three-pass counting sort per level pair: (1) find each
// sublevel cell's parent and tally children per parent, (2) prefix-sum
// those counts into contiguous slots of ha::cells_in_cells (allocated by the
// caller with one slot per cell below the top level), (3) scatter each
// sublevel cell into its parent's slot. Downstream readers (Segment 6/7 here,
// hierarchical::mag.cpp) only ever sum/mark over a parent's full child list.
// Called on the process filling the shared arrays only.
//----------------------------------------------------------------------------
void build_cells_in_cells(const std::vector<int>& cell_grid_i,
                           const std::vector<int>& cell_grid_j,
                           const std::vector<int>& cell_grid_k,
                           const std::vector<int>& level_ncy,
                           const std::vector<int>& level_ncz){

   // incrementing 1D index counter separating cells from different levels
   int index2 = 0;

   for (int level = 1; level < ha::num_levels; level ++){

      int sublevel         = level - 1;
      int start_level      = ha::cells_level_start_index[level];
      int start_sublevel   = ha::cells_level_start_index[sublevel];
      int end_level        = ha::cells_level_end_index[level];
      int end_sublevel     = ha::cells_level_end_index[sublevel];

      const int n_parents  = end_level - start_level;
      const int n_children = end_sublevel - start_sublevel;

      const int pncy = level_ncy[level];
      const int pncz = level_ncz[level];

      // pass 1: locate each child's parent (local index within `level`)
      // and count children per parent
      std::vector<int> parent_local(n_children);
      std::vector<int> child_count(n_parents, 0);

      for (int c = 0; c < n_children; c++){
         const int cell_sl = start_sublevel + c;
         const int ii = cell_grid_i[cell_sl] / 2;
         const int jj = cell_grid_j[cell_sl] / 2;
         const int kk = cell_grid_k[cell_sl] / 2;
         const int local_parent = (ii * pncy + jj) * pncz + kk;
         parent_local[c] = local_parent;
         child_count[local_parent]++;
      }

      // pass 2: prefix sum to assign each parent a contiguous slot
      std::vector<int> parent_offset(n_parents);
      int running = index2;
      for (int p = 0; p < n_parents; p++){
         parent_offset[p] = running;
         ha::cells_in_cells_start_index[start_level + p] = running;
         running += child_count[p];
         ha::cells_in_cells_end_index[start_level + p] = running;
      }
      index2 = running;

      // pass 3: scatter children into their parent's slot
      std::vector<int> fill_cursor = parent_offset;
      for (int c = 0; c < n_children; c++){
         const int cell_sl = start_sublevel + c;
         const int local_parent = parent_local[c];
         ha::cells_in_cells[fill_cursor[local_parent]++] = cell_sl;
      }
   }

}

} // end of hierarchical namespace
