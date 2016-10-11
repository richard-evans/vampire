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

         //---- DEFINED IN DIPOLE MODULE ------//
         //// Array to store list of processors and contained cells 1D
         //std::vector<int> proc_cell_index_array1D;

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
            ////fprintf(stderr,"cell_id = %d x = %f y = %f z = %f,  my_rank = %d\n",mpi_send_cells_id[lc],mpi_send_cells_pos_mom[4*lc+0],mpi_send_cells_pos_mom[4*lc+1],mpi_send_cells_pos_mom[4*lc+2],mpi_send_cells_pos_mom[4*lc+3],vmpi::my_rank);
         }

         // resize array for storing cell <-> cpu
         proc_cell_index_array1D.resize(cells_num_cells);
         // populate array with my_rank id
         for(int lc=0; lc<cells_num_cells; lc++){
            proc_cell_index_array1D[lc]=vmpi::my_rank;
         }
         ////fprintf(stderr,"proc_cell_index_array1D[cells_num_cells - 1] = %d on rank = %d\n",proc_cell_index_array1D[cells_num_cells - 1],vmpi::my_rank);

         //MPI::COMM_WORLD.Barrier();

         ////fprintf(stderr,"\n\n");

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

                     ////fprintf(stderr,"cpu %d sending num_send_cells = %d to cpu %d\n",root,num_send_cells, cpu);
                  }
               }
            }
            else{
       //// my_rank receives data form other cpus and saves them
               MPI_Recv(&num_recv_cells, 1, MPI_INT, root, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               // resize tmp recv arrays
               mpi_recv_cells_id.resize(num_recv_cells);
               mpi_recv_cells_pos_mom.resize(4*num_recv_cells);
               mpi_recv_cells_num_atoms_in_cell.resize(num_recv_cells);
               // receive data for arrays
               MPI_Recv(&mpi_recv_cells_id[0], num_recv_cells, MPI_INT, root, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_cells_pos_mom[0], 4*num_recv_cells, MPI_DOUBLE, root, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&mpi_recv_cells_num_atoms_in_cell[0], num_recv_cells, MPI_INT, root, 112, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

               ////fprintf(stderr,"cpu %d num_recv_cells = %d from cpu %d\n",vmpi::my_rank,num_recv_cells,root);

               // resize arrays for storing data
               int size   = ceil(cells_pos_and_mom_array.size()/4.0);
               int size1D = cells_cell_id_array.size();
               //////cells_local_cell_array.resize(size1D + num_recv_cells);
               //cells_pos_and_mom_array.resize(4*(size1D + num_recv_cells));
               //proc_cell_index_array1D.resize(size1D + num_recv_cells);
               cells_cell_id_array.resize(size1D + num_recv_cells);
               cells_pos_and_mom_array.resize(4*(size + num_recv_cells));
               proc_cell_index_array1D.resize(size + num_recv_cells);
               cells_num_atoms_in_cell.resize(size + num_recv_cells);
               //int size2D = cells_atom_in_cell_coords_array_x.size();
               ////fprintf(stderr,"size = %d and cells_atom_in_cell_coords_array_x.size() = %d on my_rank = %d\n",size,cells_atom_in_cell_coords_array_x.size(),vmpi::my_rank);

               // Save data into arrays
               int counter=0;
               for(int lc=size; lc<(size+num_recv_cells); lc++){
                  ////fprintf(stderr,"\tlc = %d, lc-size = %d, size+num_recv_cells = %d, cpu = %d\n",lc,lc-size,size+num_recv_cells,vmpi::my_rank);
                  cells_cell_id_array[size1D+counter] = mpi_recv_cells_id[counter];
                  //cells_local_cell_array[lc]      = mpi_recv_cells_id[counter];
                  cells_pos_and_mom_array[4*lc+0]  = mpi_recv_cells_pos_mom[4*(counter)+0];
                  cells_pos_and_mom_array[4*lc+1]  = mpi_recv_cells_pos_mom[4*(counter)+1];
                  cells_pos_and_mom_array[4*lc+2]  = mpi_recv_cells_pos_mom[4*(counter)+2];
                  cells_pos_and_mom_array[4*lc+3]  = mpi_recv_cells_pos_mom[4*(counter)+3];
                  proc_cell_index_array1D[lc]      = root;
                  cells_num_atoms_in_cell[lc]      = mpi_recv_cells_num_atoms_in_cell[counter];
                  ////fprintf(stderr,"\t>>> size = %d, cells_cell_id_array[%d] = %d, x = %f y = %f z = %f atoms_in_cell=%d on cpu = %d from root = %d proc_cell_index_array[%d] = %d\n",size+num_recv_cells,lc,cells_cell_id_array[size1D+counter],cells_pos_and_mom_array[4*lc+0],cells_pos_and_mom_array[4*lc+1],cells_pos_and_mom_array[4*lc+2],cells_num_atoms_in_cell[lc],vmpi::my_rank,root,lc,proc_cell_index_array1D[lc]);
                  counter++;
               }
            }
         }
         ////fprintf(stderr,"\n\n");
         //MPI::COMM_WORLD.Barrier();
         int new_size = cells_num_atoms_in_cell.size();
         cells_atom_in_cell_coords_array_x.resize(new_size);
         cells_atom_in_cell_coords_array_y.resize(new_size);
         cells_atom_in_cell_coords_array_z.resize(new_size);
         cells_index_atoms_array.resize(new_size);
         for(int i=cells_num_cells; i<new_size; i++){
            int n_atoms = cells_num_atoms_in_cell[i];
         cells_atom_in_cell_coords_array_x[i].resize(n_atoms);
         cells_atom_in_cell_coords_array_y[i].resize(n_atoms);
         cells_atom_in_cell_coords_array_z[i].resize(n_atoms);
         cells_index_atoms_array[i].resize(n_atoms);
         }

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
         //std::vector<int> num_send_atoms_to_cpu(vmpi::num_processors-1,0);
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


	      std::vector<MPI::Request> requests(0);
	      std::vector<MPI::Status> stati(0);

         for(int cpu=0; cpu<vmpi::num_processors; cpu++){
            num_send_atoms = 0;
            int size = ceil(cells_pos_and_mom_array.size()/4.0);
            int counter = 0;
            std::vector<bool> bool_array(cells_num_cells,1); /// bool arrays to check whether a cell has been already considered
            ////fprintf(stderr,"\nsize=%d, cells_pos_and_mom_array.size()=%lu, proc_cell_index_array1D.size()=%lu on my_rank=%d\n",size,cells_pos_and_mom_array.size(),proc_cell_index_array1D.size(),vmpi::my_rank);
            for(int lc=cells_num_cells; lc<size; lc++){
               int cpu_recv = proc_cell_index_array1D[lc];
               ////fprintf(stderr," #####<<<<< cpu=%d lc=%d cpu_recv=%d <-- my_rank=%d\n",cpu,lc,cpu_recv,vmpi::my_rank);
               if( cpu == cpu_recv){
                  for(int i=0; i<cells_num_cells; i++){
                     if(cells_num_atoms_in_cell[i]>0){
                        //reciprocal of distance between cells
                        double rij_1 = 1.0/sqrt((cells_pos_and_mom_array[4*lc+0]-cells_pos_and_mom_array[4*i+0])*(cells_pos_and_mom_array[4*lc+0]-cells_pos_and_mom_array[4*i+0]) +
                                               	(cells_pos_and_mom_array[4*lc+1]-cells_pos_and_mom_array[4*i+1])*(cells_pos_and_mom_array[4*lc+1]-cells_pos_and_mom_array[4*i+1]) +
                                               	(cells_pos_and_mom_array[4*lc+2]-cells_pos_and_mom_array[4*i+2])*(cells_pos_and_mom_array[4*lc+2]-cells_pos_and_mom_array[4*i+2]));
                        //distance between cells
                        double rij = 1.0/rij_1;
                        if(rij/cells_macro_cell_size <= dipole::cutoff){
                           if(counter < cells_num_local_cells && bool_array[i]!=0){
                              list_cpu_to_send_to.push_back(cpu_recv);
                              list_cells_to_send.push_back(i);
                              list_cells_to_recv.push_back(lc);
                              ////fprintf(stderr," counter=%d,cell %d on cpu %d is inside cutoff radius of cell %d on my_rank %d\n",counter,lc,proc_cell_index_array1D[lc],i,vmpi::my_rank);
                              counter++;
                              bool_array[i]=0;
                           }
                        } // end if statement for comparison of distance
                     } // end of if statement fo num_atoms_in_cell
                  } // end of loop over num_cells
               }
            }
         }

         //MPI::COMM_WORLD.Barrier();
         ////fprintf(stderr,"\n\n");
         //for(unsigned int i=0; i<list_cpu_to_send_to.size(); i++){
         //   //fprintf(stderr,"cpu_recv %d cell_send %d cell_recv %d on my_rank = %d\n",list_cpu_to_send_to[i],list_cells_to_send[i],list_cells_to_recv[i],vmpi::my_rank);
         //}
         ////fprintf(stderr,"\n\n saving data to send \n\n");
         //MPI::COMM_WORLD.Barrier();

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

                           ////fprintf(stderr,"\t\tnum_send_atoms = %d and num_send_cells = %d to cpu %d from rank =  %d\n",num_send_atoms,num_send_cells,cpu_recv,vmpi::my_rank);

                           // store data about recv cell
                           mpi_send_num_atoms_in_cell[counter_cells] = cells_num_atoms_in_cell[cell_send];
                           mpi_send_cells_pos_mom[4*counter_cells+0] = cells_pos_and_mom_array[4*cell_send+0];
                           mpi_send_cells_pos_mom[4*counter_cells+1] = cells_pos_and_mom_array[4*cell_send+1];
                           mpi_send_cells_pos_mom[4*counter_cells+2] = cells_pos_and_mom_array[4*cell_send+2];
                           // save tmp atoms data to be sent
                           for(int j=0; j<cells_num_atoms_in_cell[cell_send]; j++){
                              //mpi_send_num_atoms_in_cell[counter]  = cells_num_atoms_in_cell[cell_send];
                              //mpi_send_atoms_id[counter]    = cells_index_atoms_array[cell_id_array[cell_send]][j];
                              mpi_send_atoms_cell[counter]  = cell_send;
                              mpi_send_atoms_id[counter]    = cells_index_atoms_array[cell_send][j];
                              int atom_id                   = mpi_send_atoms_id[counter];
                              mpi_send_atoms_pos_x[counter] = atom_pos_x[atom_id];
                              mpi_send_atoms_pos_y[counter] = atom_pos_y[atom_id];
                              mpi_send_atoms_pos_z[counter] = atom_pos_z[atom_id];
                              int type                      = atom_type_array[atom_id];
                              const double mus              = mp::material[type].mu_s_SI/9.27400915e-24;
                              mpi_send_atoms_mom[counter]   = mus;
                              ////fprintf(stderr,"   atom %d x = %f y = %f z = %f on my_rank = %d\n",atom_id,atom_pos_x[atom_id],atom_pos_y[atom_id],atom_pos_z[atom_id],vmpi::my_rank);
                              ////fprintf(stderr,"\tcounter = %d cell_id = %d atom_id = %d  x = %f y = %f z = %f mus = %e,  my_rank = %d\n",counter,mpi_send_atoms_cell[counter],mpi_send_atoms_id[counter],mpi_send_atoms_pos_x[counter],mpi_send_atoms_pos_y[counter],mpi_send_atoms_pos_z[counter],mpi_send_atoms_mom[counter],vmpi::my_rank);
                              counter++;
                           }
                           counter_cells++;
                        }
                     }
                     ////fprintf(stderr,"counter = %d num_send_atoms = %d to cpu %d from rank =  %d\n\t mpi_send_atoms_id.size() = %d\n",counter,num_send_atoms,cpu,vmpi::my_rank,mpi_send_atoms_id.size());
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

                     //fprintf(stderr,"*****cpu %d sending num_send_atoms = %d to cpu %d******\n",cpu_send,num_send_atoms, cpu_recv);
                  } //end if(cpu!=vmpi::my_rank)
               }


            }
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


               // resize arrays for storing data
               int size = cells_atom_in_cell_coords_array_x.size();
               //////fprintf(stderr,"size = %d and cells_atom_in_cell_coords_array_x.size() = %d and cells_num_cells = %d on my_rank = %d\n",size,cells_atom_in_cell_coords_array_x.size(),cells_num_cells,vmpi::my_rank);
               int size_new = size+num_recv_cells;

               std::vector<int> recv_cell_id(num_recv_cells);
               int cell = 0;
               for(int lc=0; lc<num_recv_cells; lc++){
                  // resize arrays
                  for(unsigned int i=cells_num_cells; i<proc_cell_index_array1D.size(); i++){
                     cell = i;
                     if((mpi_recv_cells_pos_mom[4*lc+0]==cells_pos_and_mom_array[4*cell+0]) &&
                        (mpi_recv_cells_pos_mom[4*lc+1]==cells_pos_and_mom_array[4*cell+1]) &&
                        (mpi_recv_cells_pos_mom[4*lc+2]==cells_pos_and_mom_array[4*cell+2])){
                        //////fprintf(stderr," >>>><<<<<<lc=%d\ti=%d\tmy_rank=%d\n\trecv_x=%f\trecv_y=%f\trecv_z=%f\n\tx=%f\ty=%f\tz=%f\n",lc,i,vmpi::my_rank,mpi_recv_cells_pos_mom[4*lc+0],mpi_recv_cells_pos_mom[4*lc+1],mpi_recv_cells_pos_mom[4*lc+2],cells_pos_and_mom_array[4*cell+0],cells_pos_and_mom_array[4*cell+1],cells_pos_and_mom_array[4*cell+2]);
                        cells_atom_in_cell_coords_array_x[cell].resize(mpi_recv_num_atoms_in_cell[lc]);
                        cells_atom_in_cell_coords_array_y[cell].resize(mpi_recv_num_atoms_in_cell[lc]);
                        cells_atom_in_cell_coords_array_z[cell].resize(mpi_recv_num_atoms_in_cell[lc]);
                        cells_index_atoms_array[cell].resize(mpi_recv_num_atoms_in_cell[lc]);
                        recv_cell_id[lc] = cell;

                        //////fprintf(stderr,"\n cell %d is going to recv %lu data on rank=%d\n",cell,cells_index_atoms_array[cell].size(),vmpi::my_rank);
                     }
                  }
               }
               //for(unsigned int i=0; i<recv_cell_id.size(); i++){
               //   //fprintf(stderr,"\t***### i=%u, num_recv_cells=%d, recv_cell_id[%d]=%d, rank=%d ###***\n",i,num_recv_cells,i,recv_cell_id[i],vmpi::my_rank);
               //}

               ////fprintf(stderr,"\n\n >>> Resized cells_atom_in_cell_coords_array_x,y,x arrays on rank %d<<<< \n\n",vmpi::my_rank);

               int counter_atoms=0;
               for(int lc=0; lc<num_recv_cells; lc++){
                  int cell = recv_cell_id[lc];
                  cells_num_atoms_in_cell[cell] = mpi_recv_num_atoms_in_cell[lc];
                  ////fprintf(stderr,"\t***### lc=%u, cell==%d rank=%d ###***\n",lc,cell,vmpi::my_rank);
                  for(int atom=0; atom<mpi_recv_num_atoms_in_cell[lc]; atom++){
                     int cell_id = cell;
                     cells_index_atoms_array[cell][atom]           = mpi_recv_atoms_id[counter_atoms];
                     cells_atom_in_cell_coords_array_x[cell][atom] = mpi_recv_atoms_pos_x[counter_atoms];
                     cells_atom_in_cell_coords_array_y[cell][atom] = mpi_recv_atoms_pos_y[counter_atoms];
                     cells_atom_in_cell_coords_array_z[cell][atom] = mpi_recv_atoms_pos_z[counter_atoms];
                     ////fprintf(stderr,"counter_atoms = %d atom = %d atom_id = %d lc = %d x[%d][%d] = %f y[%d][%d] = %f z[%d][%d] = %f cell_id = %d num_atoms_in_cell = %d on my_rank = %d\n",counter_atoms,atom,cells_index_atoms_array[lc+size][atom],lc,lc+size,atom,cells_atom_in_cell_coords_array_x[lc+size][atom],lc+size,atom,cells_atom_in_cell_coords_array_y[lc+size][atom],lc+size,atom,cells_atom_in_cell_coords_array_z[lc+size][atom],cell_id,cells_num_atoms_in_cell[lc+size],vmpi::my_rank);
                     counter_atoms++;
                  }
               }
            } //end else statement
         }

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

         ////fprintf(stderr,"\n\n >>>>>> Starting sorting data <<<<<<<< \n\n");

         // if cells from other cpus have same coordinates => update
         std::vector<int> same_cells;
         int counter = 0;
         ////for(int i=0; i<cells_num_local_cells; i++){
            //int lc = cells_cell_id_array[i];
         int size = ceil(cells_pos_and_mom_array.size()/4.0);
         ////fprintf(stderr,"\n  *** cells_num_cells=%d siz=%d cells_index_atoms_array.size()=%lu on rank %d ***\n\n",cells_num_cells,size,cells_index_atoms_array.size(),vmpi::my_rank);
         for(int j=cells_num_cells; j<size; j++){
            for(int i=0; i<cells_num_cells; i++){
               int lc = i;
               if(((cells_pos_and_mom_array[4*lc+0]==cells_pos_and_mom_array[4*j+0]) &&
                   (cells_pos_and_mom_array[4*lc+1]==cells_pos_and_mom_array[4*j+1]) &&
                   (cells_pos_and_mom_array[4*lc+2]==cells_pos_and_mom_array[4*j+2]))&&
                   (cells_num_atoms_in_cell[j]>0)){
                  // add atoms coords to local cell with same coords
                  same_cells.push_back(j);
                  counter++;

                  ////fprintf(stderr,"cell %d\t%f %f %f on cpu %d = cell %d\t%f %f %f on cpu %d on my_rank=%d\n",lc,cells_pos_and_mom_array[4*lc+0],cells_pos_and_mom_array[4*lc+1],cells_pos_and_mom_array[4*lc+2],proc_cell_index_array1D[lc],j,cells_pos_and_mom_array[4*j+0],cells_pos_and_mom_array[4*j+1],cells_pos_and_mom_array[4*j+2],proc_cell_index_array1D[j],vmpi::my_rank);
                  ////fprintf(stderr," num_atoms_in_cell[%d] = %d on cpu = %d num_atoms_in_cell[%d] = %d on cpu = %d\n",lc,cells_num_atoms_in_cell[lc],j,cells_num_atoms_in_cell[j],proc_cell_index_array1D[j],vmpi::my_rank);

                  int size_old = cells_num_atoms_in_cell[lc];
                  cells_atom_in_cell_coords_array_x[lc].resize(size_old+cells_num_atoms_in_cell[j]);
                  cells_atom_in_cell_coords_array_y[lc].resize(size_old+cells_num_atoms_in_cell[j]);
                  cells_atom_in_cell_coords_array_z[lc].resize(size_old+cells_num_atoms_in_cell[j]);
                  cells_index_atoms_array[lc].resize(size_old+cells_num_atoms_in_cell[j]);
                  for(int k=0; k<cells_num_atoms_in_cell[j]; k++){
                     ////fprintf(stderr," num_atoms_in_cell[%d] = %d on my_rank = %d\n",lc,cells_num_atoms_in_cell[lc],vmpi::my_rank);
                     double x = cells_atom_in_cell_coords_array_x[j][k];
                     double y = cells_atom_in_cell_coords_array_y[j][k];
                     double z = cells_atom_in_cell_coords_array_z[j][k];
                     int id   = cells_index_atoms_array[j][k];
                     cells_atom_in_cell_coords_array_x[lc][k+size_old] = x;
                     cells_atom_in_cell_coords_array_y[lc][k+size_old] = y;
                     cells_atom_in_cell_coords_array_z[lc][k+size_old] = z;
                     cells_index_atoms_array[lc][k+size_old]           = id;
                     cells_num_atoms_in_cell[lc]+=1;
                     ////fprintf(stderr,"  atom = %d atom_id = %d lc = %d x[%d][%d] = %f y[%d][%d] = %f z[%d][%d] = %f num_atoms_in_cell = %d from cpu = %d\n",k+size_old,cells_index_atoms_array[j][k+size_old],lc,j,k+size_old,cells_atom_in_cell_coords_array_x[j][k+size_old],j,k+size_old,cells_atom_in_cell_coords_array_y[j][k+size_old],j,k+size_old,cells_atom_in_cell_coords_array_z[j][k+size_old],cells_num_atoms_in_cell[lc],proc_cell_index_array1D[j]);
                     ////fprintf(stderr,"  num_atoms_in_cell[%d] = %d on my_rank = %d\n",lc,cells_num_atoms_in_cell[lc],vmpi::my_rank);
                  }
               }
            }
         }

         //MPI::COMM_WORLD.Barrier();

         //for(int i=0; i<cells_num_cells; i++){
         //   if(cells_num_atoms_in_cell[i]>0){
         //      //fprintf(stderr,"  --> cell %d x = %f y = %f z= %f num_atoms = %d on my_rank = %d\n",i,cells_pos_and_mom_array[4*i+0],cells_pos_and_mom_array[4*i+1],cells_pos_and_mom_array[4*i+2],cells_num_atoms_in_cell[i],vmpi::my_rank);
         //   }
         //}

         //MPI::COMM_WORLD.Barrier();

         // print same_cells array
         //for(unsigned int i=0; i<same_cells.size(); i++){
         //   //fprintf(stderr," same_cells[%lu] = %d on my_rank = %d\n",i,same_cells[i],vmpi::my_rank);
         //}

         //MPI::COMM_WORLD.Barrier();

         ////fprintf(stderr,"\n\n >>>>>> Beore removing cellss <<<<<<<< \n\n");

         // remove cells that have same coords
         for(int i=0; i<counter; i++){
         //for(int i=counter-1;i>=0;i--){
            int lc = same_cells[i];
            ////fprintf(stderr," removed cell %d in position %d belonging to cpu = %d on my_rank = %d\n",lc,4*lc,proc_cell_index_array1D[lc],vmpi::my_rank);
            cells_pos_and_mom_array[4*lc+0] = cells_pos_and_mom_array[4*(lc+1)+0];
            cells_pos_and_mom_array[4*lc+1] = cells_pos_and_mom_array[4*(lc+1)+1];
            cells_pos_and_mom_array[4*lc+2] = cells_pos_and_mom_array[4*(lc+1)+2];
            cells_pos_and_mom_array[4*lc+3] = cells_pos_and_mom_array[4*(lc+1)+3];
            proc_cell_index_array1D[lc]     = proc_cell_index_array1D[lc+1];
            cells_num_atoms_in_cell[lc]     = cells_num_atoms_in_cell[lc+1];
         }

         ////fprintf(stderr,"\n\n >>>>>> After removing cellss <<<<<<<< \n\n");

         //MPI::COMM_WORLD.Barrier();

         ////fprintf(stderr,"\n\n >>>>>> Beore removing atoms <<<<<<<< \n\n");

         //for(int i=0; i<counter; i++){
         for(int i=counter-1;i>=0;i--){
            // reize cells_atom_in_cell_coords_array_x,y,z and cells_index_atoms_array
            int lc = same_cells[i];
            ////fprintf(stderr," removing atoms in cell %d on my_rank = %d\n",lc,vmpi::my_rank);
            cells_atom_in_cell_coords_array_x.erase(cells_atom_in_cell_coords_array_x.begin()+lc);
            cells_atom_in_cell_coords_array_y.erase(cells_atom_in_cell_coords_array_y.begin()+lc);
            cells_atom_in_cell_coords_array_z.erase(cells_atom_in_cell_coords_array_z.begin()+lc);
            cells_index_atoms_array.erase(cells_index_atoms_array.begin()+lc);
         }

         ////fprintf(stderr,"\n\n >>>>>> After removing atoms <<<<<<<< \n\n");

         //MPI::COMM_WORLD.Barrier();

         ////fprintf(stderr,"\n BEFORE: cells_pos_and_mom_array.size() =  %lu and proc_cell_index_array1D.size() = %lu on rank = %d\n",cells_pos_and_mom_array.size(),proc_cell_index_array1D.size(),vmpi::my_rank);
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

         //MPI::COMM_WORLD.Barrier();
         ////fprintf(stderr,"\n AFTER: cells_pos_and_mom_array.size() =  %lu and proc_cell_index_array1D.size() = %lu on rank = %d\n",cells_pos_and_mom_array.size(),proc_cell_index_array1D.size(),vmpi::my_rank);

         //MPI::COMM_WORLD.Barrier();

         // update cells_num_cells value
         cells_num_cells = ceil(cells_pos_and_mom_array.size()/4.0);
         ////fprintf(stderr,"\n\n new cells_num_cells = %d size_new_1D=%d on my_rank = %d\n\n",cells_num_cells,size_new_1D,vmpi::my_rank);

         //MPI::COMM_WORLD.Barrier();

         return EXIT_SUCCESS;
      }

   #endif

} // end namespace dipole
