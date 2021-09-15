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
//#include <string>
#include <cmath>
//#include <cstdlib>
#include <iostream>

// Vampire headers
#include "cells.hpp" // needed for cells::cell_id_array but to be removed
#include "dipole.hpp"
#include "vio.hpp"
#include "vutil.hpp"
#include "atoms.hpp"
#include "hierarchical.hpp"

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

   std::vector < std::vector < int > > cells_index_atoms_array;

   std::vector < int > cells_num_atoms_in_cell;
   std::vector < int > cells_num_atoms_in_cell_global;
   std::vector < int > cells_local_cell_array;
   std::vector <  double > cells_pos_and_mom_array;
   int cells_num_cells = cells::num_cells;
   int cells_num_local_cells = cells::num_local_cells;

   // allocate data for 2D arrays
   cells_index_atoms_array.resize(cells::num_cells);

   cells_num_atoms_in_cell.resize(cells::num_cells);
   cells_num_atoms_in_cell_global.resize(cells::num_cells);
   cells_local_cell_array.resize(cells::num_local_cells);
   cells_pos_and_mom_array.resize(cells::num_cells*4);

   // consider all cells
   for(int i = 0; i < cells_num_cells; i++){

      // copy number of atoms in each cell
      cells_num_atoms_in_cell[i]        = cells::num_atoms_in_cell[i];
      cells_num_atoms_in_cell_global[i] = cells::num_atoms_in_cell_global[i];


      // populate cell and moment array
      cells_pos_and_mom_array[4*i + 0] = cells::pos_and_mom_array[4*i + 0];
      cells_pos_and_mom_array[4*i + 1] = cells::pos_and_mom_array[4*i + 1];
      cells_pos_and_mom_array[4*i + 2] = cells::pos_and_mom_array[4*i + 2];
      cells_pos_and_mom_array[4*i + 3] = cells::pos_and_mom_array[4*i + 3];

      // resize second dimension of arrays
      cells_index_atoms_array[i].resize(cells::num_atoms_in_cell[i]);

      // copy a list of atoms in each cell
      for (int atom = 0; atom <cells_num_atoms_in_cell[i]; atom ++ ){
         cells_index_atoms_array[i][atom] = cells::index_atoms_array[i][atom];
      }

   }

   //
   for (int lc = 0; lc < cells_num_local_cells; lc ++){
      cells_local_cell_array[lc] = cells::local_cell_array[lc];
   }

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

   // Assign updated value of cells_num_atoms_in_cell to dipole::dipole_cells_num_atoms_in_cell. It is needed to print the config file. The actual value cells::num_atoms_in_cell is not changed instead
   dipole::dipole_cells_num_atoms_in_cell = cells_num_atoms_in_cell;

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

   //loop over all levels to calculate the positions and sizes of the cells in the levels.
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
      int ncx = static_cast<unsigned int>(ceil((system_dimensions_x+0.01)/cell_size_x));
      int ncy = static_cast<unsigned int>(ceil((system_dimensions_y+0.01)/cell_size_y));
      int ncz = static_cast<unsigned int>(ceil((system_dimensions_z+0.01)/cell_size_z));

      const int temp_num_cells = ncx*ncy*ncz;

      //set the start and end index for the level
      ha::cells_level_start_index[level] = index;
      index = index + temp_num_cells;
      ha::cells_level_end_index[level] = index;
      // std::cout << "l:\t" << level << '\t' <<  << "\t ha:\t" >> ha::cells_level_start_index[level] << '\t' << ha::cells_level_end_index[level] << '\t' << temp_num_cells << "\t" << cells_num_cells << "\t" << ncx << '\t' << ncy << '\t' << ncz << std::endl;
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
               //    if (level == 0) {
               ha::cell_positions.push_back(x);
               ha::cell_positions.push_back(y);
               ha::cell_positions.push_back(z);
               ha::cell_dimensions.push_back(size_x);
               ha::cell_dimensions.push_back(size_y);
               ha::cell_dimensions.push_back(size_z);
               //    index++; - this part to delete
               //    ha::cells_level_end_index[level] = index;
               // }
               // else if (size_x != 0 && size_y != 0 && size_z != 0) {
               //    ha::cell_positions.push_back(x);
               //    ha::cell_positions.push_back(y);
               //    ha::cell_positions.push_back(z);
               //    ha::cell_dimensions.push_back(size_x);
               //    ha::cell_dimensions.push_back(size_y);
               //    ha::cell_dimensions.push_back(size_z);
               // }
            }
         }
      }
   } //end of loop over levels

   // total number of cells is last value of the index
   ha::total_num_cells = index;

   // std::cout << "total\t" << ha::total_num_cells <<std::endl;
   // std::cout << ha::cells_level_start_index[0] << '\t'  << ha::cells_level_end_index[0] << std::endl;
   // std::cout << "A" << '\t' << ha::num_levels << '\t' << ha::cells_level_end_index[0] << '\t' << ha::cells_level_start_index[0] << '\t' << ha::total_num_cells << '\t' << std::endl;

   // resize arrays storing hierarchical cell data to correct size
   ha::cells_in_cells_start_index.resize(ha::total_num_cells,0.0);
   ha::cells_in_cells_end_index.resize(ha::total_num_cells,0.0);
   ha::cell_positions_mom.resize(ha::total_num_cells*4,0.0);

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

   // temporary data structure for debugging (RE- to delete?)
   //std::vector <int > cells_help(ha::total_num_cells, 0);

   // incrementing 1D index counter separating cells from different levels
   int index2 = 0;

   // loop over all levels
   for (int level = 1; level < ha::num_levels; level ++){

      int sublevel         = level - 1;
      int start_level      = ha::cells_level_start_index[level];
      int start_sublevel   = ha::cells_level_start_index[sublevel];
      int end_level        = ha::cells_level_end_index[level];
      int end_sublevel     = ha::cells_level_end_index[sublevel];

      // loop over all cells in level
      for (int cell_l = start_level; cell_l < end_level; cell_l++){
         ha::cells_in_cells_start_index[cell_l] = index2;
         //  std::cout << ha::cells_in_cells_start_index[cell_l] <<std::endl;
         double x = ha::cell_positions[cell_l*3 + 0];
         double y = ha::cell_positions[cell_l*3 + 1];
         double z = ha::cell_positions[cell_l*3 + 2];
         double size_x = ha::cell_dimensions[cell_l*3 + 0];
         double size_y = ha::cell_dimensions[cell_l*3 + 1];
         double size_z = ha::cell_dimensions[cell_l*3 + 2];

         double min_x = x - size_x/2.0;
         double min_y = y - size_y/2.0;
         double min_z = z - size_z/2.0;

         double max_x = x + size_x/2.0;
         double max_y = y + size_y/2.0;
         double max_z = z + size_z/2.0;

         // loop over all cells in the level below
         for (int cell_sl = start_sublevel; cell_sl < end_sublevel; cell_sl++){

            double sc_x = ha::cell_positions[cell_sl*3 + 0];
            double sc_y = ha::cell_positions[cell_sl*3 + 1];
            double sc_z = ha::cell_positions[cell_sl*3 + 2];

            // if cell is in cell at higher level, then add to the list
            if ((sc_x >= min_x) && (sc_x <= max_x) && (sc_y >= min_y) && (sc_y <= max_y) && (sc_z >= min_z) && (sc_z <= max_z)){
               //std::cout << "A" << index2 << '\t' << ha::cells_in_cells.size() <<std::endl;
               ha::cells_in_cells.push_back(cell_sl);
               index2 ++;
               ha::cells_in_cells_end_index[cell_l] = index2;
               // debugging code (RE- to remove?)
               //cells_help[cell_sl] ++;
               //if (level == 1 ){
                  //     std::cout << "cell\t" << cell_l << "\t" << ha::cell_positions[cell_l*3 + 0] << "\t" << ha::cell_positions[cell_l*3 + 1] << "\t" << ha::cell_positions[cell_l*3 + 2] <<  "\t" <<ha::cell_dimensions[cell_l*3 + 0] << "\t" <<ha::cell_dimensions[cell_l*3 + 1] <<"\t" <<ha::cell_dimensions[cell_l*3 + 2] << std::endl;
                  //           std::cout  <<"subcell" << '\t' <<  index2 << '\t' << ha::cells_in_cells_start_index[cell_l] << '\t' <<ha::cells_in_cells_end_index[cell_l] << '\t' << cell_sl <<  "\t" << ha::cell_positions[cell_sl*3 + 0] << "\t" << ha::cell_positions[cell_sl*3 + 1] << "\t" << ha::cell_positions[cell_sl*3 + 2] <<  "\t" <<ha::cell_dimensions[cell_sl*3 + 0] << "\t" <<ha::cell_dimensions[cell_sl*3 + 1] <<"\t" <<ha::cell_dimensions[cell_sl*3 + 2] << std::endl;
                  //std::cout << cell_l << "\t" << ha::cell_positions[cell_l*3 + 0] << "\t" << ha::cell_positions[cell_l*3 + 1] << "\t" << ha::cell_positions[cell_l*3 + 2] <<  "\t" <<ha::cell_dimensions[cell_l*3 + 0] << "\t" <<ha::cell_dimensions[cell_l*3 + 1] <<"\t" <<ha::cell_dimensions[cell_l*3 + 2] << std::endl;
               //}

            }

         }
      }
   }

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
   ha::num_atoms_in_cell.resize(ha::total_num_cells,0);
   ha::interaction_list_start_index.resize(ha::num_zero_level_cells,false);
   ha::interaction_list_end_index.resize(ha::num_zero_level_cells,false);
   ha::mag_array_x.resize(ha::total_num_cells,0);
   ha::mag_array_y.resize(ha::total_num_cells,0);
   ha::mag_array_z.resize(ha::total_num_cells,0);

   // check that the number of zero levels cells is correct (must be the same for hierachical and cells structures)
   if (ha::num_zero_level_cells != cells_num_cells) std::cout << "ERROR ZERO CELLS != NUM CELLS" << std::endl;
   //    std::cout << ha::num_zero_level_cells << '\t' << cells_num_cells << '\t' << ha::cell_positions_mom.size()/4. << '\t' << cells_pos_and_mom_array.size()/4. <<"\t" <<  ha::num_atoms_in_cell.size() << '\t' << cells_num_atoms_in_cell.size() <<std::endl;

   // incrementing counter determining number of atoms in cell (RE- to remove?)
   int n_cells = 0;

   for (int cell = 0; cell < cells_num_cells; cell ++){

      ha::cell_positions_mom[4*cell+0] = cells_pos_and_mom_array[4*cell+0];
      ha::cell_positions_mom[4*cell+1] = cells_pos_and_mom_array[4*cell+1];
      ha::cell_positions_mom[4*cell+2] = cells_pos_and_mom_array[4*cell+2];
      ha::cell_positions_mom[4*cell+3] = cells_pos_and_mom_array[4*cell+3];

      //    std::cout << ha::cell_positions_mom[4*cell+2] << '\t' << ha::cell_positions[3*cell+2] << std::endl;

      // determine number of atoms in each cell
      ha::num_atoms_in_cell[cell] = cells_num_atoms_in_cell[cell];

      //std::cout << ha::num_atoms_in_cell[cell] <<std::endl;
      n_cells = n_cells+ ha::num_atoms_in_cell[cell];

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

   for (int level = 1; level < ha::num_levels; level ++ ){
      n_cells = 0;
      int start = ha::cells_level_start_index[level];
      int end   = ha::cells_level_end_index[level];
      //   std::cout << "A" << level << '\t' << start << "\t" << end <<std::endl;
      for (int cell = start; cell < end; cell++){
         int start_cell_in_cell = ha::cells_in_cells_start_index[cell];
         int end_cell_in_cell = ha::cells_in_cells_end_index[cell];
         int N_cells = end_cell_in_cell - start_cell_in_cell;
         for (int interaction = start_cell_in_cell; interaction < end_cell_in_cell; interaction ++ ){
            int cell_in_cell = ha::cells_in_cells[interaction];
            // std::cout << "B" << level << '\t' << cell << '\t' << cell_in_cell << "\t" << ha::cell_positions[cell*3 + 2] <<  "\t" <<ha::cell_dimensions[cell*3 + 2] << "\t" << ha::cell_positions[cell_in_cell*3 + 2] <<  "\t" <<ha::cell_dimensions[cell_in_cell*3 + 2] << "\t" << start_cell_in_cell << '\t' << end_cell_in_cell <<std::endl;
            ha::cell_positions_mom[4*cell+0] += ha::cell_positions_mom[4*cell_in_cell+0];
            ha::cell_positions_mom[4*cell+1] += ha::cell_positions_mom[4*cell_in_cell+1];
            ha::cell_positions_mom[4*cell+2] += ha::cell_positions_mom[4*cell_in_cell+2];
            ha::cell_positions_mom[4*cell+3] += ha::cell_positions_mom[4*cell_in_cell+3];
            ha::num_atoms_in_cell[cell] += ha::num_atoms_in_cell[cell_in_cell];

         }
         ha::cell_positions_mom[4*cell+0] = ha::cell_positions_mom[4*cell+0]/N_cells;
         ha::cell_positions_mom[4*cell+1] = ha::cell_positions_mom[4*cell+1]/N_cells;
         ha::cell_positions_mom[4*cell+2] = ha::cell_positions_mom[4*cell+2]/N_cells;
         ha::cell_positions_mom[4*cell+3] = ha::cell_positions_mom[4*cell+3]/N_cells;
         // std::cout << ha::cell_positions_mom[4*cell+0] << '\t' << ha::cell_positions[3*cell+0] << "\t" << N_cells <<  std::endl;

         n_cells = n_cells + ha::num_atoms_in_cell[cell]; // this is never used (RE- to delete?)
      }
   }

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

   // data sturcture to store corner positions for each cell
   std::vector < std::vector < double> >corners;
   corners.resize(8);
   for (int i = 0; i < 8; ++i) corners[i].resize(3,0.0);

   // determine the interaction list for local cells
   int interaction_num = 0;
   for (int lc = 0; lc < cells_num_local_cells; lc ++){

      std::vector < bool > interacted(ha::total_num_cells,false);
      ha::interaction_list_start_index[lc] = interaction_num;
      int cell_i = cells::cell_id_array[lc];
      double xi = ha::cell_positions[cell_i*3 + 0];
      double yi = ha::cell_positions[cell_i*3 + 1];
      double zi = ha::cell_positions[cell_i*3 + 2];
      for (int level = ha::num_levels - 1; level > -1; level -- ){

         int start = ha::cells_level_start_index[level];
         int end   = ha::cells_level_end_index[level];
         for (int cell_j = start; cell_j < end; cell_j++){

            if (interacted[cell_j] == false && level != 0){

               double xj = ha::cell_positions[cell_j*3 + 0];
               double yj = ha::cell_positions[cell_j*3 + 1];
               double zj = ha::cell_positions[cell_j*3 + 2];
               double size_xj = ha::cell_dimensions[cell_j*3 + 0];
               double size_yj = ha::cell_dimensions[cell_j*3 + 1];
               double size_zj = ha::cell_dimensions[cell_j*3 + 2];

               corners = ha::calculate_corners(xj,yj,zj,size_xj,size_yj,size_zj);

               int is_corner_outside_range = 0;
               for (int corner = 0; corner < 8; corner ++){

                  double dx_2 = (xi - corners[corner][0])*(xi-corners[corner][0]);
                  double dy_2 = (yi - corners[corner][1])*(yi-corners[corner][1]);
                  double dz_2 = (zi - corners[corner][2])*(zi-corners[corner][2]);

                  double r = sqrt(dx_2 + dy_2 + dz_2);
                  if (r > ha::interaction_range[level]) is_corner_outside_range ++;
               }

               if (is_corner_outside_range == 8){
                  interaction_num ++;
                  ha::interaction_list_end_index[lc] = interaction_num;
                  interacted[cell_j] = true;
                  int start_cell_in_cell = ha::cells_in_cells_start_index[cell_j];
                  int end_cell_in_cell = ha::cells_in_cells_end_index[cell_j];
                  ha::interaction_list.push_back(cell_j);

                  //if (cell_i == 0 ha::cell_positions[cell_j*3 + 0]  == 115)  std::cout <<"interacted" << std::endl;//std::cout << cell_j <<  "\t" << ha::cell_positions[cell_j*3 + 0] << "\t" << ha::cell_positions[cell_j*3 + 1] << "\t" << ha::cell_positions[cell_j*3 + 2] <<  "\t" <<ha::cell_dimensions[cell_j*3 + 0] << "\t" <<ha::cell_dimensions[cell_j*3 + 1] <<"\t" <<ha::cell_dimensions[cell_j*3 + 2] << std::endl;
                  ha::interaction_list_end_index[lc] = interaction_num;
                  // if (cell_i == 0) std::cout << "CELL POSITION\t" << start_cell_in_cell << '\t' << end_cell_in_cell << '\t' << cell_j << '\t' << ha::cell_positions[cell_j*3 + 0] << "\t" << ha::cell_positions[cell_j*3 + 1] << "\t" << ha::cell_positions[cell_j*3 + 2] <<  "\t" <<ha::cell_dimensions[cell_j*3 + 0] << "\t" <<ha::cell_dimensions[cell_j*3 + 1] <<"\t" <<ha::cell_dimensions[cell_j*3 + 2] << std::endl;
                  for (int interaction = start_cell_in_cell; interaction < end_cell_in_cell; interaction ++ ){
                     int sub_cell = ha::cells_in_cells[interaction];
                     interacted[sub_cell] = true;
                     //    if (cell_i == 0)     std::cout << "SUBCELL POSITION\t" << sub_cell <<  "\t" << ha::cell_positions[sub_cell*3 + 0] << "\t" << ha::cell_positions[sub_cell*3 + 1] << "\t" << ha::cell_positions[sub_cell*3 + 2] <<  "\t" <<ha::cell_dimensions[sub_cell*3 + 0] << "\t" <<ha::cell_dimensions[sub_cell*3 + 1] <<"\t" <<ha::cell_dimensions[sub_cell*3 + 2] << std::endl;
                  }
               }
               //if (cell_i == 0)  std::cout << '\t' << std::endl;
            }
            else if (level != 0){
               int start_cell_in_cell = ha::cells_in_cells_start_index[cell_j];
               int end_cell_in_cell = ha::cells_in_cells_end_index[cell_j];

               for (int interaction = start_cell_in_cell; interaction < end_cell_in_cell; interaction ++ ){
                  int sub_cell = ha::cells_in_cells[interaction];
                  interacted[sub_cell] = true;
               }
            }
            else if (level == 0 && interacted[cell_j] == false){
               ha::interaction_list.push_back(cell_j);
               interaction_num ++;
               ha::interaction_list_end_index[lc] = interaction_num;
            }
         }

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

   // resize mu_0*Hd-field cells array
   dipole::cells_mu0Hd_field_array_x.resize(cells_num_cells,0.0);
   dipole::cells_mu0Hd_field_array_y.resize(cells_num_cells,0.0);
   dipole::cells_mu0Hd_field_array_z.resize(cells_num_cells,0.0);

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

            // only calculate interaction if there are atoms in remote cell - I dont understand either of these conditions!
            // 1) why does the local number of atoms in the cell matter - surely the global number is more meaningful?
            // 2) why would I only want to calculate interactions with zero-level cells??? surely we want all interactions?
            if ( cells_num_atoms_in_cell_global[cell_j] != 0 && cell_j < ha::num_zero_level_cells){

               //--------------------------------------------------------------
               // Calculation of inter part of dipolar tensor
               //--------------------------------------------------------------
               // check that remote cell is actually magnetic
               //if ( cells_num_atoms_in_cell_global[cell_j] > 0 ){

                  // check if cell is not the same
                  if( cell_i != cell_j ){
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

   // print out detailed hierarchical timings
   zlog << zTs() << "\tSegment times for hierarchical initialisation:" << std::endl;
   for(int i=0; i<9; i++) zlog << zTs() << "\t\t Segment time " << i+1 << ": " << segment_time[i] << " s" << std::endl;

   // Woohoo we made it!
   return;

}

} // end of hierarchical namespace
