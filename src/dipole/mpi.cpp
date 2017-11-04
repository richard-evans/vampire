//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo 2016. All rights reserved.
//
//   Email: am1808@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

// Vampire headers
#include "cells.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

#include "atoms.hpp"

// cells module headers
#include "internal.hpp"

namespace dipole{
   //-----------------------------------------------------------------------------
   // Function to parallelise dipolar field calculation
   //-----------------------------------------------------------------------------

namespace internal{

   #ifdef MPICF
      /*------------------------------------------------*/
      /*Function to send and receive data between cells */
      /*------------------------------------------------*/
      int send_recv_cells_data(std::vector<int>& proc_cell_index_array1D,
                              std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_x,
                              std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_y,
                              std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_z,
                              std::vector< std::vector <int> >& cells_index_atoms_array,
                              std::vector<double>& cells_pos_and_mom_array,
                              std::vector<int>& cells_num_atoms_in_cell,
                              std::vector<int>& cells_cell_id_array,
                              std::vector<int>& cells_local_cell_array,
                              int cells_num_local_cells,
                              int cells_num_cells
                              ){
         // temporary variables to send and receive data
         int num_send_cells;
         std::vector<int> mpi_send_cells_num_atoms_in_cell;
         std::vector<int> mpi_send_cells_id;
         std::vector<double> mpi_send_cells_pos_mom;

         int num_recv_cells;
         std::vector<int> mpi_recv_cells_num_atoms_in_cell;
         std::vector<int> mpi_recv_cells_id;
         std::vector<double> mpi_recv_cells_pos_mom;

         // my_rank allocate memory to arrays for sending data
         num_send_cells = cells_num_local_cells;
         mpi_send_cells_id.resize(num_send_cells);
         mpi_send_cells_pos_mom.resize(4*num_send_cells);
         mpi_send_cells_num_atoms_in_cell.resize(num_send_cells);
         // populate temporary arrays with data
         for(int lc=0; lc<num_send_cells; lc++){
            int cell = cells_cell_id_array[lc];
            mpi_send_cells_id[lc]            = cells_cell_id_array[lc];
            mpi_send_cells_pos_mom[4*lc+0]   = cells_pos_and_mom_array[4*cell+0];
            mpi_send_cells_pos_mom[4*lc+1]   = cells_pos_and_mom_array[4*cell+1];
            mpi_send_cells_pos_mom[4*lc+2]   = cells_pos_and_mom_array[4*cell+2];
            mpi_send_cells_pos_mom[4*lc+3]   = cells_pos_and_mom_array[4*cell+3];
            mpi_send_cells_num_atoms_in_cell[lc] = cells_num_atoms_in_cell[cell];
         }

         // resize array for storing cell <-> cpu
         proc_cell_index_array1D.resize(cells_num_cells);
         // populate array with my_rank id
         for(int lc=0; lc<cells_num_cells; lc++){
            proc_cell_index_array1D[lc]=vmpi::my_rank;
         }

         for(int proc=0; proc<vmpi::num_processors; proc++){
            int root = proc;
            if(vmpi::my_rank == root){
               // my_rank send data to other cpus
               for(int cpu=0; cpu<vmpi::num_processors; cpu++){
                  if(cpu != vmpi::my_rank ){
                     MPI_Send(&num_send_cells, 1, MPI_INT, cpu, 100, MPI_COMM_WORLD);
                     MPI_Send(&mpi_send_cells_id[0], num_send_cells, MPI_INT, cpu, 101, MPI_COMM_WORLD);
                     MPI_Send(&mpi_send_cells_pos_mom[0], 4*num_send_cells, MPI_DOUBLE, cpu, 102, MPI_COMM_WORLD);
                     MPI_Send(&mpi_send_cells_num_atoms_in_cell[0], num_send_cells, MPI_INT, cpu, 112, MPI_COMM_WORLD);
                  }
               }
            }
            else{
               MPI_Recv(&num_recv_cells, 1, MPI_INT, root, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               // resize tmp recv arrays
               mpi_recv_cells_id.resize(num_recv_cells);
               mpi_recv_cells_pos_mom.resize(4*num_recv_cells);
               mpi_recv_cells_num_atoms_in_cell.resize(num_recv_cells);
               // receive data for arrays
               MPI_Recv(&mpi_recv_cells_id[0], num_recv_cells, MPI_INT, root, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_cells_pos_mom[0], 4*num_recv_cells, MPI_DOUBLE, root, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_cells_num_atoms_in_cell[0], num_recv_cells, MPI_INT, root, 112, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

               // resize arrays for storing data
               int size   = ceil(cells_pos_and_mom_array.size()/4.0);
               int size1D = cells_cell_id_array.size();
               cells_cell_id_array.resize(size1D + num_recv_cells);
               cells_pos_and_mom_array.resize(4*(size + num_recv_cells));
               proc_cell_index_array1D.resize(size + num_recv_cells);
               cells_num_atoms_in_cell.resize(size + num_recv_cells);

               // Save data into arrays
               int counter=0;
               for(int lc=size; lc<(size+num_recv_cells); lc++){
                  cells_cell_id_array[size1D+counter] = mpi_recv_cells_id[counter];
                  cells_pos_and_mom_array[4*lc+0]  = mpi_recv_cells_pos_mom[4*(counter)+0];
                  cells_pos_and_mom_array[4*lc+1]  = mpi_recv_cells_pos_mom[4*(counter)+1];
                  cells_pos_and_mom_array[4*lc+2]  = mpi_recv_cells_pos_mom[4*(counter)+2];
                  cells_pos_and_mom_array[4*lc+3]  = mpi_recv_cells_pos_mom[4*(counter)+3];
                  proc_cell_index_array1D[lc]      = root;

                  counter++;
               }
            }
         }
         int new_size = ceil(cells_pos_and_mom_array.size()/4.);
         cells_atom_in_cell_coords_array_x.resize(new_size);
         cells_atom_in_cell_coords_array_y.resize(new_size);
         cells_atom_in_cell_coords_array_z.resize(new_size);
         cells_index_atoms_array.resize(new_size);

         //Deallocate memory used for mpi arrays
         mpi_send_cells_id.clear();
         mpi_send_cells_pos_mom.clear();
         mpi_send_cells_num_atoms_in_cell.clear();
         mpi_recv_cells_id.clear();
         mpi_recv_cells_pos_mom.clear();
         mpi_recv_cells_num_atoms_in_cell.clear();

         return EXIT_SUCCESS;
      }

      /*------------------------------------------------*/
      /*Function to send and receive data atoms between cpus */
      /*------------------------------------------------*/
      int send_recv_atoms_data(std::vector<int>& proc_cell_index_array1D,
                              std::vector<int>& cells_cell_id_array,
                              std::vector<int>& cells_local_cell_array,
                              std::vector<double>& atom_pos_x,
                              std::vector<double>& atom_pos_y,
                              std::vector<double>& atom_pos_z,
                              std::vector<int>& atom_type_array, // atomic moments (from dipole;:internal::atom_type_array)
                              std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_x,
                              std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_y,
                              std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_z,
                              std::vector< std::vector <int> >& cells_index_atoms_array,
                              std::vector<double>& cells_pos_and_mom_array,
                              std::vector<int>& cells_num_atoms_in_cell,
                              int cells_num_local_cells,
                              int cells_num_cells,
                              double cells_macro_cell_size
                              ){

         // temporary variables to send and receive data
         std::vector<int> list_cpu_to_send_to;
         std::vector<int> list_cells_to_send;
         std::vector<int> list_cells_to_recv;
         int num_send_atoms;
         int num_send_cells;
         std::vector<int> mpi_send_atoms_cell;
         std::vector<int> mpi_send_num_atoms_in_cell;
         std::vector<int> mpi_send_atoms_id;
         std::vector<double> mpi_send_atoms_pos_x;
         std::vector<double> mpi_send_atoms_pos_y;
         std::vector<double> mpi_send_atoms_pos_z;
         std::vector<double> mpi_send_atoms_mom;
         std::vector<double> mpi_send_cells_pos_mom;

         int num_recv_atoms;
         int num_recv_cells;
         std::vector<int> mpi_recv_atoms_cell;
         std::vector<int> mpi_recv_num_atoms_in_cell;
         std::vector<int> mpi_recv_atoms_id;
         std::vector<double> mpi_recv_atoms_pos_x;
         std::vector<double> mpi_recv_atoms_pos_y;
         std::vector<double> mpi_recv_atoms_pos_z;
         std::vector<double> mpi_recv_atoms_mom;
         std::vector<double> mpi_recv_cells_pos_mom;


         for(int cpu=0; cpu<vmpi::num_processors; cpu++){
            num_send_atoms = 0;
            int size = ceil(cells_pos_and_mom_array.size()/4.0);
            int counter = 0;
            std::vector<bool> bool_array(cells_num_cells,1); /// bool arrays to check whether a cell has been already considered
            for(int lc=cells_num_cells; lc<size; lc++){
               int cpu_recv = proc_cell_index_array1D[lc];
               if( cpu == cpu_recv){
                  for(int i=0; i<cells_num_cells; i++){
                     if(cells_num_atoms_in_cell[i]>0){
                        //reciprocal of distance between cells
                        double rij_1 = 1.0/sqrt((cells_pos_and_mom_array[4*lc+0]-cells_pos_and_mom_array[4*i+0])*(cells_pos_and_mom_array[4*lc+0]-cells_pos_and_mom_array[4*i+0]) +
                                               	(cells_pos_and_mom_array[4*lc+1]-cells_pos_and_mom_array[4*i+1])*(cells_pos_and_mom_array[4*lc+1]-cells_pos_and_mom_array[4*i+1]) +
                                               	(cells_pos_and_mom_array[4*lc+2]-cells_pos_and_mom_array[4*i+2])*(cells_pos_and_mom_array[4*lc+2]-cells_pos_and_mom_array[4*i+2]));
                        //distance between cells
                        double rij = 1.0/rij_1;
                        if((rij/cells_macro_cell_size <= dipole::cutoff) || ((cells_pos_and_mom_array[4*lc+0]==cells_pos_and_mom_array[4*i+0]) && (cells_pos_and_mom_array[4*lc+1]==cells_pos_and_mom_array[4*i+1]) && (cells_pos_and_mom_array[4*lc+2]==cells_pos_and_mom_array[4*i+2]) ) ){
                           if(bool_array[i]!=0){
                              list_cpu_to_send_to.push_back(cpu_recv);
                              list_cells_to_send.push_back(i);
                              list_cells_to_recv.push_back(lc);
                              counter++;
                              bool_array[i]=0;
                           }
                        } // end if statement for comparison of distance
                     } // end of if statement fo num_atoms_in_cell
                  } // end of loop over num_cells
               }
            }
         }

         for(int proc=0; proc<vmpi::num_processors; proc++){
            if(vmpi::my_rank == proc){

               for(int cpu=0; cpu<vmpi::num_processors; cpu++){
                  int counter = 0;
                  int counter_cells = 0;
                  if(cpu != vmpi::my_rank){
                     num_send_atoms = 0;
                     num_send_cells = 0;
                     for(unsigned int i=0; i<list_cpu_to_send_to.size(); i++){
                        int cpu_recv = list_cpu_to_send_to[i];
                        int cell_send = list_cells_to_send[i];
                        if(cpu_recv==cpu){
                           num_send_cells++;
                           num_send_atoms += cells_num_atoms_in_cell[cell_send];
                           // resize arrays
                           mpi_send_atoms_id.resize((num_send_atoms));
                           mpi_send_atoms_cell.resize((num_send_atoms));
                           mpi_send_atoms_pos_x.resize((num_send_atoms));
                           mpi_send_atoms_pos_y.resize((num_send_atoms));
                           mpi_send_atoms_pos_z.resize((num_send_atoms));
                           mpi_send_atoms_mom.resize((num_send_atoms));
                           mpi_send_num_atoms_in_cell.resize((num_send_cells));
                           mpi_send_cells_pos_mom.resize((4*num_send_cells));

                           // store data about recv cell
                           mpi_send_num_atoms_in_cell[counter_cells] = cells_num_atoms_in_cell[cell_send];
                           mpi_send_cells_pos_mom[4*counter_cells+0] = cells_pos_and_mom_array[4*cell_send+0];
                           mpi_send_cells_pos_mom[4*counter_cells+1] = cells_pos_and_mom_array[4*cell_send+1];
                           mpi_send_cells_pos_mom[4*counter_cells+2] = cells_pos_and_mom_array[4*cell_send+2];
                           // save tmp atoms data to be sent
                           for(int j=0; j<cells_num_atoms_in_cell[cell_send]; j++){
                              mpi_send_atoms_cell[counter]  = cell_send;
                              mpi_send_atoms_id[counter]    = cells_index_atoms_array[cell_send][j];
                              int atom_id                   = mpi_send_atoms_id[counter];
                              mpi_send_atoms_pos_x[counter] = atom_pos_x[atom_id];
                              mpi_send_atoms_pos_y[counter] = atom_pos_y[atom_id];
                              mpi_send_atoms_pos_z[counter] = atom_pos_z[atom_id];
                              int type                      = atom_type_array[atom_id];
                              const double mus              = mp::material[type].mu_s_SI/9.27400915e-24;
                              mpi_send_atoms_mom[counter]   = mus;
                              counter++;
                           }
                           counter_cells++;
                        }
                     }
                     // Send num_of_atoms to be received and allocate memory in arrays
                     // my_rank send data to other cpus
                     int cpu_send = vmpi::my_rank;
                     int cpu_recv = cpu;
                     MPI_Send(&num_send_cells,           1,                MPI_INT,    cpu_recv, 111, MPI_COMM_WORLD);
                     MPI_Send(&num_send_atoms,           1,                MPI_INT,    cpu_recv, 103, MPI_COMM_WORLD);
                     MPI_Send(&mpi_send_atoms_id[0],     num_send_atoms,   MPI_INT,    cpu_recv, 104, MPI_COMM_WORLD);
                     MPI_Send(&mpi_send_atoms_pos_x[0],  num_send_atoms,   MPI_DOUBLE, cpu_recv, 105, MPI_COMM_WORLD);
                     MPI_Send(&mpi_send_atoms_pos_y[0],  num_send_atoms,   MPI_DOUBLE, cpu_recv, 106, MPI_COMM_WORLD);
                     MPI_Send(&mpi_send_atoms_pos_z[0],  num_send_atoms,   MPI_DOUBLE, cpu_recv, 107, MPI_COMM_WORLD);
                     MPI_Send(&mpi_send_atoms_mom[0],    num_send_atoms,   MPI_DOUBLE, cpu_recv, 108, MPI_COMM_WORLD);
                     MPI_Send(&mpi_send_atoms_cell[0],   num_send_atoms,   MPI_INT,    cpu_recv, 109, MPI_COMM_WORLD);
                     MPI_Send(&mpi_send_num_atoms_in_cell[0],  num_send_cells,   MPI_INT,    cpu_recv, 110, MPI_COMM_WORLD);
                     MPI_Send(&mpi_send_cells_pos_mom[0],  4*num_send_cells,   MPI_DOUBLE,   cpu_recv, 113, MPI_COMM_WORLD);
                  } //end if(cpu!=vmpi::my_rank)
               }


            }
            // if I am not =proc, then I recv data
            else{
               //fprintf(stderr," >>> rank %d is receiving data from proc %d <<<<< \n",vmpi::my_rank,proc);
               int cpu_send = proc;
               int cpu_recv = vmpi::my_rank;
               // my_rank receives number of objects that hav been sent from other cpus
               MPI_Recv(&num_recv_cells, 1, MPI_INT, cpu_send, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&num_recv_atoms, 1, MPI_INT, cpu_send, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

               // resize tmp recv arrays
               mpi_recv_atoms_id.resize(num_recv_atoms);
               mpi_recv_atoms_pos_x.resize(num_recv_atoms);
               mpi_recv_atoms_pos_y.resize(num_recv_atoms);
               mpi_recv_atoms_pos_z.resize(num_recv_atoms);
               mpi_recv_atoms_mom.resize(num_recv_atoms);
               mpi_recv_atoms_cell.resize(num_recv_atoms);
               mpi_recv_num_atoms_in_cell.resize(num_recv_cells);
               mpi_recv_cells_pos_mom.resize(4*num_recv_cells);

               // receive data for arrays
               MPI_Recv(&mpi_recv_atoms_id[0],     num_recv_atoms, MPI_INT,      cpu_send, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_atoms_pos_x[0],  num_recv_atoms, MPI_DOUBLE,   cpu_send, 105, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_atoms_pos_y[0],  num_recv_atoms, MPI_DOUBLE,   cpu_send, 106, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_atoms_pos_z[0],  num_recv_atoms, MPI_DOUBLE,   cpu_send, 107, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_atoms_mom[0],    num_recv_atoms, MPI_DOUBLE,   cpu_send, 108, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_atoms_cell[0],   num_recv_atoms, MPI_INT,      cpu_send, 109, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_num_atoms_in_cell[0],   num_recv_cells, MPI_INT,      cpu_send, 110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_cells_pos_mom[0],   4*num_recv_cells,   MPI_DOUBLE,   cpu_send, 113, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


               std::vector<int> recv_cell_id(num_recv_cells);
               std::vector<int> old_size_array(num_recv_cells,0);
               std::vector<bool> bool_array(num_recv_cells,1); /// bool arrays to check whether a cell has been already considered
               int cell = 0;
               for(int lc=0; lc<num_recv_cells; lc++){
                  int cell_recv_counter=0;
                  // resize arrays
                  for(unsigned int i=cells_num_cells; i<proc_cell_index_array1D.size(); i++){
                     cell = i;
                     int size = cells_index_atoms_array[cell].size();
                     old_size_array[lc]=size;
                     if((mpi_recv_cells_pos_mom[4*lc+0]==cells_pos_and_mom_array[4*cell+0]) &&
                        (mpi_recv_cells_pos_mom[4*lc+1]==cells_pos_and_mom_array[4*cell+1]) &&
                        (mpi_recv_cells_pos_mom[4*lc+2]==cells_pos_and_mom_array[4*cell+2]) && bool_array[lc]!=0){

                        recv_cell_id[lc] = cell;

                        cell_recv_counter++;
                        bool_array[lc]=0;
                     }
                     else if((mpi_recv_cells_pos_mom[4*lc+0]==cells_pos_and_mom_array[4*cell+0]) &&
                             (mpi_recv_cells_pos_mom[4*lc+1]==cells_pos_and_mom_array[4*cell+1]) &&
                             (mpi_recv_cells_pos_mom[4*lc+2]==cells_pos_and_mom_array[4*cell+2]) && bool_array[lc]==0){
                        cells_num_atoms_in_cell[cell]=0;
                     }
                  }
               }

               int counter_atoms=0;
               for(int lc=0; lc<num_recv_cells; lc++){
                  // if there is only one cell with that coord
                  int cell = recv_cell_id[lc];
                  int old_size = cells_index_atoms_array[cell].size();
                  cells_num_atoms_in_cell[cell] += mpi_recv_num_atoms_in_cell[lc];

                  cells_atom_in_cell_coords_array_x[cell].resize(mpi_recv_num_atoms_in_cell[lc]+old_size);
                  cells_atom_in_cell_coords_array_y[cell].resize(mpi_recv_num_atoms_in_cell[lc]+old_size);
                  cells_atom_in_cell_coords_array_z[cell].resize(mpi_recv_num_atoms_in_cell[lc]+old_size);
                  cells_index_atoms_array[cell].resize(mpi_recv_num_atoms_in_cell[lc]+old_size);

                  for(int atom=0; atom<mpi_recv_num_atoms_in_cell[lc]; atom++){
                     int id = atom+old_size;
                     cells_index_atoms_array[cell][id]           = mpi_recv_atoms_id[counter_atoms];
                     cells_atom_in_cell_coords_array_x[cell][id] = mpi_recv_atoms_pos_x[counter_atoms];
                     cells_atom_in_cell_coords_array_y[cell][id] = mpi_recv_atoms_pos_y[counter_atoms];
                     cells_atom_in_cell_coords_array_z[cell][id] = mpi_recv_atoms_pos_z[counter_atoms];

                     counter_atoms++;
                  }
               }
            } //end else statement
         }

         // Free memory
         list_cpu_to_send_to.clear();
         list_cells_to_send.clear();
         list_cells_to_recv.clear();
         mpi_send_atoms_cell.clear();
         mpi_send_num_atoms_in_cell.clear();
         mpi_send_atoms_id.clear();
         mpi_send_atoms_pos_x.clear();
         mpi_send_atoms_pos_y.clear();
         mpi_send_atoms_pos_z.clear();
         mpi_send_atoms_mom.clear();
         mpi_send_cells_pos_mom.clear();
         mpi_recv_atoms_cell.clear();
         mpi_recv_num_atoms_in_cell.clear();
         mpi_recv_atoms_id.clear();
         mpi_recv_atoms_pos_x.clear();
         mpi_recv_atoms_pos_y.clear();
         mpi_recv_atoms_pos_z.clear();
         mpi_recv_atoms_mom.clear();
         mpi_recv_cells_pos_mom.clear();

         return EXIT_SUCCESS;
      }

      /*------------------------------------------------*/
      /*Function to sort cells/atoms data after sharing */
      /*------------------------------------------------*/
      int sort_data(std::vector<int>& proc_cell_index_array1D,
                  std::vector<int>& cells_cell_id_array,
                  std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_x,
                  std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_y,
                  std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_z,
                  std::vector< std::vector <int> >& cells_index_atoms_array,
                  std::vector<double>& cells_pos_and_mom_array,
                  std::vector<int>& cells_num_atoms_in_cell,
                  int cells_num_local_cells,
                  int cells_num_cells
                  ){

         std::vector<int> same_cells;
         int counter = 0;
         int size = ceil(cells_pos_and_mom_array.size()/4.0);
         for(int j=cells_num_cells; j<size; j++){
            for(int i=0; i<cells_num_cells; i++){
               int lc = i;
               // add atoms coords to local cell with same coords
               if(((cells_pos_and_mom_array[4*lc+0]==cells_pos_and_mom_array[4*j+0]) &&
                   (cells_pos_and_mom_array[4*lc+1]==cells_pos_and_mom_array[4*j+1]) &&
                   (cells_pos_and_mom_array[4*lc+2]==cells_pos_and_mom_array[4*j+2]))&&
                   (cells_num_atoms_in_cell[j]>0)){
                  same_cells.push_back(j);
                  counter++;

                  int size_old = cells_num_atoms_in_cell[lc];
                  cells_atom_in_cell_coords_array_x[lc].resize(size_old+cells_num_atoms_in_cell[j]);
                  cells_atom_in_cell_coords_array_y[lc].resize(size_old+cells_num_atoms_in_cell[j]);
                  cells_atom_in_cell_coords_array_z[lc].resize(size_old+cells_num_atoms_in_cell[j]);
                  cells_index_atoms_array[lc].resize(size_old+cells_num_atoms_in_cell[j]);
                  for(int k=0; k<cells_num_atoms_in_cell[j]; k++){
                     double x = cells_atom_in_cell_coords_array_x[j][k];
                     double y = cells_atom_in_cell_coords_array_y[j][k];
                     double z = cells_atom_in_cell_coords_array_z[j][k];
                     int id   = cells_index_atoms_array[j][k];
                     cells_atom_in_cell_coords_array_x[lc][k+size_old] = x;
                     cells_atom_in_cell_coords_array_y[lc][k+size_old] = y;
                     cells_atom_in_cell_coords_array_z[lc][k+size_old] = z;
                     cells_index_atoms_array[lc][k+size_old]           = id;
                     cells_num_atoms_in_cell[lc]+=1;
                  }
               }
            }
         }

         // remove cells that have same coords
         for(int i=0; i<counter; i++){
            int lc = same_cells[i];

            cells_pos_and_mom_array[4*lc+0] = cells_pos_and_mom_array[4*(lc+1)+0];
            cells_pos_and_mom_array[4*lc+1] = cells_pos_and_mom_array[4*(lc+1)+1];
            cells_pos_and_mom_array[4*lc+2] = cells_pos_and_mom_array[4*(lc+1)+2];
            cells_pos_and_mom_array[4*lc+3] = cells_pos_and_mom_array[4*(lc+1)+3];
            proc_cell_index_array1D[lc]     = proc_cell_index_array1D[lc+1];
            cells_num_atoms_in_cell[lc]     = cells_num_atoms_in_cell[lc+1];
         }

         for(int i=counter-1;i>=0;i--){
            // reize cells_atom_in_cell_coords_array_x,y,z and cells_index_atoms_array
            int lc = same_cells[i];

            cells_atom_in_cell_coords_array_x.erase(cells_atom_in_cell_coords_array_x.begin()+lc);
            cells_atom_in_cell_coords_array_y.erase(cells_atom_in_cell_coords_array_y.begin()+lc);
            cells_atom_in_cell_coords_array_z.erase(cells_atom_in_cell_coords_array_z.begin()+lc);
            cells_index_atoms_array.erase(cells_index_atoms_array.begin()+lc);
         }

         //resize cells_pos_and_mom_array and proc_cell_index_array1D
         size       = cells_pos_and_mom_array.size();
         int size_new   = size - 4*counter;
         int size_new_1D= ceil(size_new/4.0);
         cells_pos_and_mom_array.resize(size_new);
         proc_cell_index_array1D.resize(size_new_1D);
         cells_num_atoms_in_cell.resize(size_new_1D);
         cells_atom_in_cell_coords_array_x.resize(size_new_1D);
         cells_atom_in_cell_coords_array_y.resize(size_new_1D);
         cells_atom_in_cell_coords_array_z.resize(size_new_1D);
         cells_index_atoms_array.resize(size_new_1D);

         // update cells_num_cells value
         cells_num_cells = ceil(cells_pos_and_mom_array.size()/4.0);

         return EXIT_SUCCESS;
      }

      /*--------------------------------------------------------*/
      /*Function to send cells demag factor to be computed      */
      /*--------------------------------------------------------*/
      int send_cells_demag_factor(std::vector<int>& cells_cell_id_array,
                                 std::vector<double>& N_tensor_array,
                                 int cells_num_local_cells
                                 ){

         // allocate memory to send data
         int num_send_cells = cells_num_local_cells;                          // number of objects to send
         std::vector<int> mpi_send_cells_id(num_send_cells,0);                // cells id to be sent
         std::vector<double> mpi_send_cells_demag_factor(6*num_send_cells,0.0);      // demag-factor and self term array to be sent

         // loop over local cells to save data to send
         for(int i=0; i<num_send_cells; i++){
            mpi_send_cells_id[i] = cells_cell_id_array[i];
            //store demag factor
            mpi_send_cells_demag_factor[6*i+0] = N_tensor_array[6*mpi_send_cells_id[i]+0];
            mpi_send_cells_demag_factor[6*i+1] = N_tensor_array[6*mpi_send_cells_id[i]+1];
            mpi_send_cells_demag_factor[6*i+2] = N_tensor_array[6*mpi_send_cells_id[i]+2];
            mpi_send_cells_demag_factor[6*i+3] = N_tensor_array[6*mpi_send_cells_id[i]+3];
            mpi_send_cells_demag_factor[6*i+4] = N_tensor_array[6*mpi_send_cells_id[i]+4];
            mpi_send_cells_demag_factor[6*i+5] = N_tensor_array[6*mpi_send_cells_id[i]+5];
         }

         // // Uncomment in case you want to check exchange of data between cores
         // std::cerr << "/* Data allocated on rank */" << vmpi::my_rank << '\t';

         // send cells id, demag factors and self term
         MPI_Send(&num_send_cells, 1, MPI_INT, 0, 120, MPI_COMM_WORLD);
         MPI_Send(&mpi_send_cells_id[0], num_send_cells, MPI_INT, 0, 121, MPI_COMM_WORLD);
         MPI_Send(&mpi_send_cells_demag_factor[0], 6*num_send_cells, MPI_DOUBLE, 0, 122, MPI_COMM_WORLD);

         // loop over CPUs
         for(int cpu=0; cpu<vmpi::num_processors; cpu++){
            // if I am root proc, then I receive cells demag factors
            if(vmpi::my_rank == 0){
               // allocate int to receive num of cells to be received
               int num_recv_cells;
               // Receive num_recv_cells
               MPI_Recv(&num_recv_cells, 1, MPI_INT, cpu, 120, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               // Allocate arrays for demag factor
               std::vector<int> mpi_recv_cells_id(num_recv_cells,0);
               std::vector<double> mpi_recv_cells_demag_factor(6*num_recv_cells,0.0);
               // Receive arrays
               MPI_Recv(&mpi_recv_cells_id[0], num_recv_cells, MPI_INT, cpu, 121, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_cells_demag_factor[0], 6*num_recv_cells, MPI_DOUBLE, cpu, 122, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

               // Save received data
               for(int i=0; i<num_recv_cells; i++){
                  int lc = mpi_recv_cells_id[i];
                  // save demag factor
                  N_tensor_array[6*lc+0] = mpi_recv_cells_demag_factor[6*i+0];
                  N_tensor_array[6*lc+1] = mpi_recv_cells_demag_factor[6*i+1];
                  N_tensor_array[6*lc+2] = mpi_recv_cells_demag_factor[6*i+2];
                  N_tensor_array[6*lc+3] = mpi_recv_cells_demag_factor[6*i+3];
                  N_tensor_array[6*lc+4] = mpi_recv_cells_demag_factor[6*i+4];
                  N_tensor_array[6*lc+5] = mpi_recv_cells_demag_factor[6*i+5];
               } // end saving data

               // // Uncomment in case you want to check exchange of data between cores
               // std::cerr << "/* \tData received on rank */" << vmpi::my_rank << '\t';

            } // end if I am root proc
         } // end loop over cpus

         // // Uncomment in case you want to check exchange of data between cores
         // int counter_total_cells_non_zero = 0;
         // if(vmpi::my_rank==0){
         //    std::cout << "I am rank 0 and I will print the tensor for each cell\n" << std::endl;
         //    // Print tensor
         //    for (int i = 0; i < cells::num_cells; i++)
         //    {
         //       if (dipole::internal::cells_num_atoms_in_cell[i] > 0)
         //       {
         //          std::cout << "\t*----------------------------------*" << std::endl;
         //          std::cout << "\tCell = " << i << "\tNat_cell_i = " << dipole::internal::cells_num_atoms_in_cell[i]  << "\n";
         //          std::cout << "\t" << N_tensor_array[6*i+0] << "\t" << N_tensor_array[6*i+1] << "\t" << N_tensor_array[6*i+2] << "\n";
         //          std::cout << "\t" << N_tensor_array[6*i+1] << "\t" << N_tensor_array[6*i+3] << "\t" << N_tensor_array[6*i+4] << "\n";
         //          std::cout << "\t" << N_tensor_array[6*i+2] << "\t" << N_tensor_array[6*i+4] << "\t" << N_tensor_array[6*i+5] << "\n";
         //          std::cout << "\t*----------------------------------*" << std::endl;
         //          std::cout << std::endl;
         //
         //          // increment counter for number of cells
         //          counter_total_cells_non_zero ++;
         //       }
         //    }
         // }

         // // Uncomment in case you want to check exchange of data between cores
         // std::cout << "Number of cells non zero on rank 0 = " << counter_total_cells_non_zero << std::endl;

         // Clear arrays used only for sending data
         mpi_send_cells_id.clear();
         mpi_send_cells_demag_factor.clear();

         return EXIT_SUCCESS;
      }

#endif

   } // end of namespace internal

   #ifdef MPICF
      /*--------------------------------------------------------*/
      /*Function to send cells field to be output in cfg file   */
      /*--------------------------------------------------------*/
      int send_cells_field(std::vector<int>& cells_cell_id_array,
                           std::vector<double>& dipole_cells_field_array_x, //B-field
                           std::vector<double>& dipole_cells_field_array_y,
                           std::vector<double>& dipole_cells_field_array_z,
                           int cells_num_local_cells
                  ){

         // allocate memory to send data
         int num_send_cells = cells_num_local_cells;                          // number of objects to send
         std::vector<int> mpi_send_cells_id(num_send_cells,0);                // cells id to be sent
         std::vector<double> mpi_send_cells_field(3*num_send_cells,0.0);      // B-field array to be sent

         // loop over local cells to save data to send
         for(int i=0; i<num_send_cells; i++){
            mpi_send_cells_id[i] = cells_cell_id_array[i];
            //store B-field
            mpi_send_cells_field[3*i+0] = dipole_cells_field_array_x[mpi_send_cells_id[i]];
            mpi_send_cells_field[3*i+1] = dipole_cells_field_array_y[mpi_send_cells_id[i]];
            mpi_send_cells_field[3*i+2] = dipole_cells_field_array_z[mpi_send_cells_id[i]];
         }

         // send cells id, B-field, Hd-field
         MPI_Send(&num_send_cells, 1, MPI_INT, 0, 114, MPI_COMM_WORLD);
         MPI_Send(&mpi_send_cells_id[0], num_send_cells, MPI_INT, 0, 115, MPI_COMM_WORLD);
         MPI_Send(&mpi_send_cells_field[0], 3*num_send_cells, MPI_DOUBLE, 0, 116, MPI_COMM_WORLD);

         // loop over CPUs
         for(int cpu=0; cpu<vmpi::num_processors; cpu++){
            // if I am root proc, then I receive cells dipolar field
            if(vmpi::my_rank == 0){
               // allocate int to receive num of cells to be received
               int num_recv_cells;
               // Receive num_recv_cells
               MPI_Recv(&num_recv_cells, 1, MPI_INT, cpu, 114, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               // Allocate arrays for field
               std::vector<int> mpi_recv_cells_id(num_recv_cells,0);
               std::vector<double> mpi_recv_cells_field(3*num_recv_cells,0.0);
               // Receive arrays
               MPI_Recv(&mpi_recv_cells_id[0], num_recv_cells, MPI_INT, cpu, 115, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_cells_field[0], 3*num_recv_cells, MPI_DOUBLE, cpu, 116, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

               // Save received data
               for(int i=0; i<num_recv_cells; i++){
                  int lc = mpi_recv_cells_id[i];
                  // save B-field
                  dipole_cells_field_array_x[lc] = mpi_recv_cells_field[3*i+0];
                  dipole_cells_field_array_y[lc] = mpi_recv_cells_field[3*i+1];
                  dipole_cells_field_array_z[lc] = mpi_recv_cells_field[3*i+2];
               } // end saving data
            } // end if I am root proc
         } // end loop over cpus

         // Clear arrays used only for sending data
         mpi_send_cells_id.clear();
         mpi_send_cells_field.clear();

         return EXIT_SUCCESS;
      }


   #endif

} // end namespace dipole
