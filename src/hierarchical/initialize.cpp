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
// C++ standard library headers
#include <string>
#include <cmath>
#include <cstdlib>
#include <iostream>
// C library headers
#include <fenv.h>
#include <signal.h>
#include <math.h>
// Vampire headers
#include "cells.hpp" // needed for cells::cell_id_array but to be removed
#include "dipole.hpp"
#include "vio.hpp"
#include "vutil.hpp"
#include "atoms.hpp"
// Vampire headers
#include "hierarchical.hpp"


// hierarchical module headers
#include "internal.hpp"
#include "../dipole/internal.hpp"

namespace ha = hierarchical::internal;

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

int min(int x, int y, int z){
    int min = x;
    if (y < min)
      min = y;
    if (z > min)
      min = z;

   return min;
}


namespace hierarchical{

   //----------------------------------------------------------------------------
   // Function to initialize hierarchical module
   //----------------------------------------------------------------------------
   void initialize(double system_dimensions_x, double system_dimensions_y, double system_dimensions_z){

      // instantiate timer
      vutil::vtimer_t timer;


      //  start timer
       timer.start();

       int A = ceil(std::log(2*system_dimensions_x/(dipole::cutoff*2*cells::macro_cell_size_x))/std::log(2.0));
       int B = ceil(std::log(2*system_dimensions_y/(dipole::cutoff*2*cells::macro_cell_size_y))/std::log(2.0));
       int C = ceil(std::log(2*system_dimensions_z/(dipole::cutoff*2*cells::macro_cell_size_z))/std::log(2.0));
   //   std::cout << dipole::cutoff*cells::macro_cell_size <<std::endl;
        //caclualte the number of levels necessary for the calculation
        int N = largest(A, B,C);
        int cell_size = 0;
        if (N == 0)
          ha::num_levels = A;
        if (N == 1)
          ha::num_levels = B;
        if (N == 2)
          ha::num_levels = C;

        ha::cells_level_start_index.resize(ha::num_levels,0.0);
        ha::cells_level_end_index.resize(ha::num_levels,0.0);
        ha::interaction_range.resize(ha::num_levels,0.0);

        // Call parallelisation function
        // Exchange cells data

       std::vector < std::vector < double > > cells_atom_in_cell_coords_array_x;
       std::vector < std::vector < double > > cells_atom_in_cell_coords_array_y;
       std::vector < std::vector < double > > cells_atom_in_cell_coords_array_z;
       std::vector < std::vector < int > > cells_index_atoms_array;

       std::vector < int > cells_num_atoms_in_cell;
       std::vector < int > cells_local_cell_array;
       std::vector <  double > cells_pos_and_mom_array;
       int cells_num_cells = cells::num_cells;
       int cells_num_local_cells = cells::num_local_cells;

       cells_atom_in_cell_coords_array_x.resize(cells::num_cells);
       cells_atom_in_cell_coords_array_y.resize(cells::num_cells);
       cells_atom_in_cell_coords_array_z.resize(cells::num_cells);
       cells_index_atoms_array.resize(cells::num_cells);

       cells_num_atoms_in_cell.resize(cells::num_cells);
       cells_local_cell_array.resize(cells::num_local_cells);
       cells_pos_and_mom_array.resize(cells::num_cells*4);




          for(int i=0; i<cells_num_cells; i++){
             // resize arrays
             cells_num_atoms_in_cell[i] = cells::num_atoms_in_cell[i];
             cells_pos_and_mom_array[4*i + 0] = cells::pos_and_mom_array[4*i + 0];
             cells_pos_and_mom_array[4*i + 1] = cells::pos_and_mom_array[4*i + 1];
             cells_pos_and_mom_array[4*i + 2] = cells::pos_and_mom_array[4*i + 2];
             cells_pos_and_mom_array[4*i + 3] = cells::pos_and_mom_array[4*i + 3];
             cells_atom_in_cell_coords_array_x[i].resize(cells::num_atoms_in_cell[i]);
             cells_atom_in_cell_coords_array_y[i].resize(cells::num_atoms_in_cell[i]);
             cells_atom_in_cell_coords_array_z[i].resize(cells::num_atoms_in_cell[i]);
             cells_index_atoms_array[i].resize(cells::num_atoms_in_cell[i]);
             //std::cout << i << "\t" << cells_num_atoms_in_cell[i] << "\t" <<  std::endl;
             for (int atom = 0; atom <cells_num_atoms_in_cell[i]; atom ++ ){
                cells_index_atoms_array[i][atom] = cells::index_atoms_array[i][atom];
                cells_atom_in_cell_coords_array_x[i][atom] = cells::atom_in_cell_coords_array_x[i][atom];
                cells_atom_in_cell_coords_array_y[i][atom] = cells::atom_in_cell_coords_array_y[i][atom];
                cells_atom_in_cell_coords_array_z[i][atom] = cells::atom_in_cell_coords_array_z[i][atom];
             }

          }

          for (int lc = 0; lc < cells_num_local_cells; lc ++){
             cells_local_cell_array[lc] = cells::local_cell_array[lc];
          }
         timer.stop();
         std::cout << "\t1 [ " << timer.elapsed_time() << " s ]" << std::endl;
         zlog << zTs() <<  "\tDIPOLE UPDATE1. Time taken: " << timer.elapsed_time() << " s"<< std::endl;


                timer.start();

                #ifdef MPICF

                const int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;

                std::vector<double> atom_pos_x(num_local_atoms,0.0);
                std::vector<double> atom_pos_y(num_local_atoms,0.0);
                std::vector<double> atom_pos_z(num_local_atoms,0.0);

                for(int atom=0; atom<num_local_atoms; atom++){
                   atom_pos_x[atom]=atoms::x_coord_array[atom];
                   atom_pos_y[atom]=atoms::y_coord_array[atom];
                   atom_pos_z[atom]=atoms::z_coord_array[atom];

                   //if (vmpi::my_rank == 0) std::cerr << atom << '\t' << num_local_atoms << '\t' <<atom_pos_x[atom] << '\t' << atom_pos_y[atom] <<std::endl;

                }


              dipole::internal::send_recv_cells_data(dipole::internal::proc_cell_index_array1D,
                                                     cells_atom_in_cell_coords_array_x,
                                                     cells_atom_in_cell_coords_array_y,
                                                     cells_atom_in_cell_coords_array_z,
                                                     cells_index_atoms_array,
                                                     cells_pos_and_mom_array,
                                                     cells_num_atoms_in_cell,
                                                     cells::cell_id_array,
                                                     cells_local_cell_array,
                                                     cells_num_local_cells,
                                                     cells_num_cells);

             // Exchange atoms data
             dipole::internal::send_recv_atoms_data(dipole::internal::proc_cell_index_array1D,
                                                    cells::cell_id_array,
                                                    cells_local_cell_array,
                                                    atom_pos_x,
                                                    atom_pos_y,
                                                    atom_pos_z,
                                                    dipole::internal::atom_type_array, // atomic moments (from dipole;:internal::atom_type_array)
                                                    cells_atom_in_cell_coords_array_x,
                                                    cells_atom_in_cell_coords_array_y,
                                                    cells_atom_in_cell_coords_array_z,
                                                    cells_index_atoms_array,
                                                    cells_pos_and_mom_array,
                                                    cells_num_atoms_in_cell,
                                                    cells_num_local_cells,
                                                    cells_num_cells,
                                                    cells::macro_cell_size);


                dipole::internal::sort_data(dipole::internal::proc_cell_index_array1D,
                                            cells::cell_id_array,
                                            cells_atom_in_cell_coords_array_x,
                                            cells_atom_in_cell_coords_array_y,
                                            cells_atom_in_cell_coords_array_z,
                                            cells_index_atoms_array,
                                            cells_pos_and_mom_array,
                                            cells_num_atoms_in_cell,
                                            cells_num_local_cells,
                                            cells_num_cells);

                        //                    std::cout <<"CELLS" <<  cells::num_atoms_in_cell_global.size() <<std::endl;

                //After transferring the data across cores, assign value cells_num_atoms_in_cell[] from cells_num_atoms_in_cell_global[]
                for(unsigned int i=0; i<cells::num_atoms_in_cell_global.size(); i++){
                   if(cells::num_atoms_in_cell_global[i]>0 && cells_num_atoms_in_cell[i]==0){
                      cells_num_atoms_in_cell[i] = cells::num_atoms_in_cell_global[i];
                     // std::cout <<"CELLS" <<  cells_atom_in_cell_coords_array_x.size() <<std::endl;
                   }
                }

                // Clear atom_pos_x,y,z
                atom_pos_x.clear();
                atom_pos_y.clear();
                atom_pos_z.clear();

              #endif

             timer.stop();
             std::cout << "\t2 [ " << timer.elapsed_time() << " s ]" << std::endl;
             zlog << zTs() <<  "\tDIPOLE UPDATE2. Time taken: " << timer.elapsed_time() << " s"<< std::endl;


            timer.start();

              int index = 0;
              internal::av_cell_size = cells::macro_cell_size_x*cells::macro_cell_size_x + cells::macro_cell_size_y*cells::macro_cell_size_y + cells::macro_cell_size_z*cells::macro_cell_size_z;
              internal::av_cell_size = sqrt(internal::av_cell_size);
              //loop over all levels to calvculate the positions and sizes of the cells in the levels.
              for (int level = 0; level < ha::num_levels; level ++){

                 double cell_size_x = pow(2,level)*cells::macro_cell_size_x;
                 double cell_size_y = pow(2,level)*cells::macro_cell_size_y;
                 double cell_size_z = pow(2,level)*cells::macro_cell_size_z;

                 ha::interaction_range[level] = internal::av_cell_size*level*dipole::cutoff + dipole::cutoff;
                 if (level == ha::num_levels - 1 ) ha::interaction_range[level] = cell_size*dipole::cutoff*10000;
                 // Calculate number of microcells
                 // determine number of cells in x and y (global)
                 int ncx = static_cast<unsigned int>(ceil((system_dimensions_x+0.01)/cell_size_x));
                 int ncy = static_cast<unsigned int>(ceil((system_dimensions_y+0.01)/cell_size_y));
                 int ncz = static_cast<unsigned int>(ceil((system_dimensions_z+0.01)/cell_size_z));

                 int temp_num_cells = ncx*ncy*ncz;
                 //set the start and end index for the level
                 ha::cells_level_start_index[level] = index;
                 index = index + temp_num_cells;
                 ha::cells_level_end_index[level] = index;
                 // std::cout << "l:\t" << level << '\t' <<  << "\t ha:\t" >> ha::cells_level_start_index[level] << '\t' << ha::cells_level_end_index[level] << '\t' << temp_num_cells << "\t" << cells_num_cells << "\t" << ncx << '\t' << ncy << '\t' << ncz << std::endl;
                 double size_x,size_y,size_z;
                 double x,y,z;
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
                     //    if (level == 0) {
                            ha::cell_positions.push_back(x);
                            ha::cell_positions.push_back(y);
                            ha::cell_positions.push_back(z);
                            ha::cell_dimensions.push_back(size_x);
                            ha::cell_dimensions.push_back(size_y);
                            ha::cell_dimensions.push_back(size_z);
                        //    index++;
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
                } //loop over levels
                ha::total_num_cells = index;
                int index2 = 0;



         //   std::cout << "total\t" << ha::total_num_cells <<std::endl;
            //  std::cout << ha::cells_level_start_index[0] << '\t'  << ha::cells_level_end_index[0] << std::endl;

                ha::cells_in_cells_start_index.resize(ha::total_num_cells,0.0);
                ha::cells_in_cells_end_index.resize(ha::total_num_cells,0.0);
                ha::cell_positions_mom.resize(ha::total_num_cells*4,0.0);
            //   std::cout <<"total:\t" <<  ha::total_num_cells << '\t' << " zero:" << '\t' << cells_num_cells << std::endl;
               timer.stop();
               std::cout << "\t3 [ " << timer.elapsed_time() << " s ]" << std::endl;
               zlog << zTs() <<  "\tDIPOLE UPDATE3. Time taken: " << timer.elapsed_time() << " s"<< std::endl;


            //   std::cout << "A" << '\t' << ha::num_levels << '\t' << ha::cells_level_end_index[0] << '\t' << ha::cells_level_start_index[0] << '\t' << ha::total_num_cells << '\t' << std::endl;
              timer.start();

              std::vector <int > cells_help(ha::total_num_cells, 0);

                 for (int level = 1; level < ha::num_levels; level ++){

                  int sublevel         = level - 1;
                  int start_level      = ha::cells_level_start_index[level];
                  int start_sublevel   = ha::cells_level_start_index[sublevel];
                  int end_level        = ha::cells_level_end_index[level];
                  int end_sublevel     = ha::cells_level_end_index[sublevel];

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

                     for (int cell_sl = start_sublevel; cell_sl < end_sublevel; cell_sl++){

                      double sc_x = ha::cell_positions[cell_sl*3 + 0];
                      double sc_y = ha::cell_positions[cell_sl*3 + 1];
                      double sc_z = ha::cell_positions[cell_sl*3 + 2];

                      if ((sc_x >= min_x) && (sc_x <= max_x) && (sc_y >= min_y) && (sc_y <= max_y) && (sc_z >= min_z) && (sc_z <= max_z)){
                       //std::cout << "A" << index2 << '\t' << ha::cells_in_cells.size() <<std::endl;
                       ha::cells_in_cells.push_back(cell_sl);
                       index2 ++;
                       ha::cells_in_cells_end_index[cell_l] = index2;

                       //cells_help[cell_sl] ++;
                       if (level == 1 ){
                  //     std::cout << "cell\t" << cell_l << "\t" << ha::cell_positions[cell_l*3 + 0] << "\t" << ha::cell_positions[cell_l*3 + 1] << "\t" << ha::cell_positions[cell_l*3 + 2] <<  "\t" <<ha::cell_dimensions[cell_l*3 + 0] << "\t" <<ha::cell_dimensions[cell_l*3 + 1] <<"\t" <<ha::cell_dimensions[cell_l*3 + 2] << std::endl;


               //           std::cout  <<"subcell" << '\t' <<  index2 << '\t' << ha::cells_in_cells_start_index[cell_l] << '\t' <<ha::cells_in_cells_end_index[cell_l] << '\t' << cell_sl <<  "\t" << ha::cell_positions[cell_sl*3 + 0] << "\t" << ha::cell_positions[cell_sl*3 + 1] << "\t" << ha::cell_positions[cell_sl*3 + 2] <<  "\t" <<ha::cell_dimensions[cell_sl*3 + 0] << "\t" <<ha::cell_dimensions[cell_sl*3 + 1] <<"\t" <<ha::cell_dimensions[cell_sl*3 + 2] << std::endl;
                         //std::cout << cell_l << "\t" << ha::cell_positions[cell_l*3 + 0] << "\t" << ha::cell_positions[cell_l*3 + 1] << "\t" << ha::cell_positions[cell_l*3 + 2] <<  "\t" <<ha::cell_dimensions[cell_l*3 + 0] << "\t" <<ha::cell_dimensions[cell_l*3 + 1] <<"\t" <<ha::cell_dimensions[cell_l*3 + 2] << std::endl;
               }

                      }

                     }
                  }
                }


               timer.stop();
               std::cout << "\t4 [ " << timer.elapsed_time() << " s ]" << std::endl;
               zlog << zTs() <<  "\tDIPOLE UPDATE4. Time taken: " << timer.elapsed_time() << " s"<< std::endl;


              timer.start();


                ha::num_zero_level_cells = ha::cells_level_end_index[0];
                //std::cout << ha::num_zero_level_cells << std::endl;
                ha::num_atoms_in_cell.resize(ha::total_num_cells,0);
                ha::interaction_list_start_index.resize(ha::num_zero_level_cells,false);
                ha::interaction_list_end_index.resize(ha::num_zero_level_cells,false);
                ha::mag_array_x.resize(ha::total_num_cells,0);
                ha::mag_array_y.resize(ha::total_num_cells,0);
                ha::mag_array_z.resize(ha::total_num_cells,0);

               if (ha::num_zero_level_cells != cells_num_cells) std::cout << "ERROR ZERO CELLS != NUM CELLS" << std::endl;
                int n_cells = 0;
            //    std::cout << ha::num_zero_level_cells << '\t' << cells_num_cells << '\t' << ha::cell_positions_mom.size()/4. << '\t' << cells_pos_and_mom_array.size()/4. <<"\t" <<  ha::num_atoms_in_cell.size() << '\t' << cells_num_atoms_in_cell.size() <<std::endl;
                for (int cell = 0; cell < cells_num_cells; cell ++){

                   ha::cell_positions_mom[4*cell+0] = cells_pos_and_mom_array[4*cell+0];
                   ha::cell_positions_mom[4*cell+1] = cells_pos_and_mom_array[4*cell+1];
                   ha::cell_positions_mom[4*cell+2] = cells_pos_and_mom_array[4*cell+2];
                   ha::cell_positions_mom[4*cell+3] = cells_pos_and_mom_array[4*cell+3];
               //    std::cout << ha::cell_positions_mom[4*cell+2] << '\t' << ha::cell_positions[3*cell+2] << std::endl;
                   ha::num_atoms_in_cell[cell] = cells_num_atoms_in_cell[cell];
                   //std::cout << ha::num_atoms_in_cell[cell] <<std::endl;
                   n_cells = n_cells+ ha::num_atoms_in_cell[cell];
                }

               timer.stop();
               std::cout << "\t5 [ " << timer.elapsed_time() << " s ]" << std::endl;
               zlog << zTs() <<  "\tDIPOLE UPDATE5. Time taken: " << timer.elapsed_time() << " s"<< std::endl;


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

                      n_cells = n_cells + ha::num_atoms_in_cell[cell];
                   }
                 }

                timer.stop();
                std::cout << "\t6 [ " << timer.elapsed_time() << " s ]" << std::endl;
                zlog << zTs() <<  "\tDIPOLE UPDATE6. Time taken: " << timer.elapsed_time() << " s"<< std::endl;


              timer.start();


                 std::vector < std::vector < double> >corners;
                 corners.resize(8);


                 for (int i = 0; i < 8; ++i) corners[i].resize(3,0.0);

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

               timer.stop();
               std::cout << "\t7 [ " << timer.elapsed_time() << " s ]" << std::endl;
               zlog << zTs() <<  "\tDIPOLE UPDATE7. Time taken: " << timer.elapsed_time() << " s"<< std::endl;


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


               timer.stop();
               std::cout << "\t8 [ " << timer.elapsed_time() << " s ]" << std::endl;
               zlog << zTs() <<  "\tDIPOLE UPDATE8. Time taken: " << timer.elapsed_time() << " s"<< std::endl;


               timer.start();

               //
                for(int lc=0; lc<cells_num_local_cells; lc++){
                //   start timer
                   int cell_i = cells::cell_id_array[lc];

                   if (cells_num_atoms_in_cell[cell_i] != 0){

                      const int start = ha::interaction_list_start_index[lc];
                      const int end = ha::interaction_list_end_index[lc];

                  //    std::cerr << cell_i << '\t' << "interaction" << "\t" << end - start << "\t" << cells_num_cells << "\t" << cells_num_atoms_in_cell[cell_i] << std::endl;

                      for(int j = start;j<end;j++){

                         int cell_j = ha::interaction_list[j];
                        // if (cell_i == 0) std::cout << cell_j <<  "\t" << ha::cell_positions[cell_j*3 + 0] << "\t" << ha::cell_positions[cell_j*3 + 1] << "\t" << ha::cell_positions[cell_j*3 + 2] <<  "\t" <<ha::cell_dimensions[cell_j*3 + 0] << "\t" <<ha::cell_dimensions[cell_j*3 + 1] <<"\t" <<ha::cell_dimensions[cell_j*3 + 2] << std::endl;
                         if ( cells_num_atoms_in_cell[cell_j] != 0 && cell_j < ha::num_zero_level_cells){
                        //   if ( cells_num_atoms_in_cell[cell_j] == 0 && cell_i == 0) std::cout << "HHHHM" << cell_j <<  "\t" << cells_num_atoms_in_cell[cell_j] << std::endl;
                            if(cell_i!=cell_j) {
                               ha::calc_inter(cell_i,cell_j,j, cells_atom_in_cell_coords_array_x, cells_atom_in_cell_coords_array_y, cells_atom_in_cell_coords_array_z);
                            }
                            else if (cell_i == cell_j){
                              ha::calc_intra(cell_i,cell_j,j, cells_atom_in_cell_coords_array_x, cells_atom_in_cell_coords_array_y, cells_atom_in_cell_coords_array_z);
                           }
                            if (ha::rij_tensor_xx[j]*ha::rij_tensor_xx[j] < 1e-15) ha::rij_tensor_xx[j] =0;
                            if (ha::rij_tensor_xy[j]*ha::rij_tensor_xy[j] < 1e-15) ha::rij_tensor_xy[j] =0;
                            if (ha::rij_tensor_xz[j]*ha::rij_tensor_xz[j] < 1e-15) ha::rij_tensor_xz[j] =0;
                            if (ha::rij_tensor_yy[j]*ha::rij_tensor_yy[j] < 1e-15) ha::rij_tensor_yy[j] =0;
                            if (ha::rij_tensor_yz[j]*ha::rij_tensor_yz[j] < 1e-15) ha::rij_tensor_yz[j] =0;
                            if (ha::rij_tensor_zz[j]*ha::rij_tensor_zz[j] < 1e-15) ha::rij_tensor_zz[j] =0;
                           //if (cell_i == 0 ) std::cout << "A" << '\t' << cell_i << '\t' << cell_j << "\t" << ha::rij_tensor_xx[j] << '\t' << ha::rij_tensor_xy[j] << '\t' << ha::rij_tensor_xz[j] << '\t' << ha::rij_tensor_yy[j] << '\t' << ha::rij_tensor_yz[j] <<'\t' << ha::rij_tensor_zz[j] << std::endl;

               //             //std::cout << cells_num_atoms_in_cell[cell_i] << '\t' << cells_num_atoms_in_cell[cell_j] << '\t' <<  cell_i << '\t' << cell_j << "\t" <<  ha::rij_tensor_xx[j] << '\t' << ha::rij_tensor_xy[j] << '\t' << ha::rij_tensor_xz[j] << '\t' << ha::rij_tensor_yy[j] << '\t' << ha::rij_tensor_yz[j] <<'\t' << ha::rij_tensor_zz[j] << std::endl;
                       }
                      }
                   }
                }
               timer.stop();
               std::cout << "\t9 [ " << timer.elapsed_time() << " s ]" << std::endl;
               zlog << zTs() <<  "\tDIPOLE UPDATE9. Time taken: " << timer.elapsed_time() << " s"<< std::endl;

         //   std::cout << "A" << std::endl;
      return;

   }

} // end of hierarchical namespace
