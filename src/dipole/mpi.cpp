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
#include "vutil.hpp"

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

        // initiate timer
        vutil::vtimer_t timer;

        // start timer
        timer.start();

        // temporary variables to send and receive data
        std::vector<int> list_cpu_to_send_to;
        std::vector<int> list_cells_to_send;
        std::vector<int> list_cells_to_recv;
        int num_send_atoms;
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
        std::vector < int > mpi_recv_atoms_cell;
        std::vector < int > mpi_recv_num_atoms_in_cell;
        std::vector < int > mpi_recv_atoms_id;
        std::vector < double > mpi_recv_atoms_pos_x;
        std::vector < double > mpi_recv_atoms_pos_y;
        std::vector < double > mpi_recv_atoms_pos_z;
        std::vector < double > mpi_recv_atoms_mom;
        std::vector < double > mpi_recv_cells_pos_mom;

        std::vector < std::vector < int > > mpi_2d_send_atoms_id;
        std::vector < std::vector < int > > mpi_2d_send_atoms_cell;
        std::vector < std::vector < int > > mpi_2d_send_num_atoms_in_cell;
        std::vector < std::vector < double > > mpi_2d_send_atoms_pos_x;
        std::vector < std::vector < double > > mpi_2d_send_atoms_pos_y;
        std::vector < std::vector < double > > mpi_2d_send_atoms_pos_z;
        std::vector < std::vector < double > > mpi_2d_send_atoms_mom;
        std::vector < std::vector < double > > mpi_2d_send_cells_pos_mom;

        mpi_2d_send_atoms_id.resize(vmpi::num_processors);
        mpi_2d_send_atoms_cell.resize(vmpi::num_processors);
        mpi_2d_send_num_atoms_in_cell.resize(vmpi::num_processors);
        mpi_2d_send_atoms_pos_x.resize(vmpi::num_processors);
        mpi_2d_send_atoms_pos_y.resize(vmpi::num_processors);
        mpi_2d_send_atoms_pos_z.resize(vmpi::num_processors);
        mpi_2d_send_atoms_mom.resize(vmpi::num_processors);
        mpi_2d_send_cells_pos_mom.resize(vmpi::num_processors);

        std::vector <int>  receive_counts(vmpi::num_processors,0);
        std::vector <int>  receive_counts_cell(vmpi::num_processors,0);
        std::vector <int>  receive_counts_2(vmpi::num_processors,0);
        std::vector <int>  receive_displacements(vmpi::num_processors,0);
        std::vector <int>  receive_displacements_cells(vmpi::num_processors,0);
        std::vector <int>  receive_displacements_four_cells(vmpi::num_processors,0);
        std::vector <int>  counter(vmpi::num_processors,0);
        std::vector <int>  counter_cells(vmpi::num_processors,0);
        std::vector <int>  counter_four_cells(vmpi::num_processors,0);
        std::vector <int>  one_count(vmpi::num_processors,vmpi::num_processors);
        std::vector <int>  one_displacements(vmpi::num_processors,0);
        std::vector <int>  recv_counter(vmpi::num_processors,0);
        std::vector <int>  recv_counter_cells(vmpi::num_processors,0);
        std::vector <int>  recv_counter_four_cells(vmpi::num_processors,0);
        std::vector <int>  final_recieve_counter(vmpi::num_processors,0);
        std::vector <int>  final_recieve_counter_cells(vmpi::num_processors,0);
        std::vector <int>  final_recieve_counter_four_cells(vmpi::num_processors,0);

        //loop to calcualte if two cells are within the dipole cut off range.
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


         // Calculate total storage
         vmpi::barrier();
         //can maybe be replaced by proc == vmpir::my_rank
         for(int proc=0; proc<vmpi::num_processors; proc++){
            if(vmpi::my_rank == proc){
               //loops over cells tou are sending form your cpu (cpu send).
               // for(unsigned int i=0; i<list_cpu_to_send_to.size(); i++){
               //
               //    int cpu_recv = list_cpu_to_send_to[i];
               //    int cell_send = list_cells_to_send[i];
               //
               //    receive_counts[cpu_recv] += cells_num_atoms_in_cell[cell_send];
               //    receive_counts_cell[cpu_recv] += 1;
               // }

               // for (int i = 0; i < vmpi::num_processors; ++i)   {
               //    // resize arrays
               //    mpi_2d_send_atoms_id[i].resize(receive_counts[cpu_recv]);
               //    mpi_2d_send_atoms_pos_x[i].resize(receive_counts[cpu_recv]);
               //    mpi_2d_send_atoms_pos_y[i].resize(receive_counts[cpu_recv]);
               //    mpi_2d_send_atoms_pos_z[i].resize(receive_counts[cpu_recv]);
               //    mpi_2d_send_atoms_mom[i].resize(receive_counts[cpu_recv]);
               //    mpi_2d_send_atoms_cell[i].resize(receive_counts[cpu_recv]);
               //    mpi_2d_send_num_atoms_in_cell[i].resize(receive_counts[cpu_recv]);
               //    mpi_2d_send_cells_pos_mom[i].resize(4*receive_counts_cell[cpu_recv]);
               // }

            //}

               //loops over cells tou are sending form your cpu (cpu send).
               for(unsigned int i=0; i<list_cpu_to_send_to.size(); i++){

                  //determine cpu to send cell too
                  int cpu_recv = list_cpu_to_send_to[i];
                  //which cell youre sending
                  int cell_send = list_cells_to_send[i];

                  //determine total number of cells/atoms sent to each processor so far,
                  receive_counts[cpu_recv] += cells_num_atoms_in_cell[cell_send];
                  receive_counts_cell[cpu_recv] += 1;
               }

               for (int i = 0; i < vmpi::num_processors; ++i)   {
                  // resize arrays
                  mpi_2d_send_atoms_id[i].resize(receive_counts[i]);
                  mpi_2d_send_atoms_pos_x[i].resize(receive_counts[i]);
                  mpi_2d_send_atoms_pos_y[i].resize(receive_counts[i]);
                  mpi_2d_send_atoms_pos_z[i].resize(receive_counts[i]);
                  mpi_2d_send_atoms_mom[i].resize(receive_counts[i]);
                  mpi_2d_send_atoms_cell[i].resize(receive_counts[i]);
                  mpi_2d_send_num_atoms_in_cell[i].resize(receive_counts[i]);
                  mpi_2d_send_cells_pos_mom[i].resize(4*receive_counts_cell[i]);
               }

               for(unsigned int i=0; i<list_cpu_to_send_to.size(); i++){

                  //determine cpu to send cell too
                  int cpu_recv = list_cpu_to_send_to[i];
                  //which cell youre sending
                  int cell_send = list_cells_to_send[i];

                  // store data about recv cell
                  //stored in a 2d array for each cpu that receives data - for the gather command later.
                  int N_cell                                        = counter_cells[cpu_recv];
                  mpi_2d_send_num_atoms_in_cell[cpu_recv][N_cell]   = cells_num_atoms_in_cell[cell_send];
                  mpi_2d_send_cells_pos_mom[cpu_recv][4*N_cell+0]   = cells_pos_and_mom_array[4*cell_send+0];
                  mpi_2d_send_cells_pos_mom[cpu_recv][4*N_cell+1]   = cells_pos_and_mom_array[4*cell_send+1];
                  mpi_2d_send_cells_pos_mom[cpu_recv][4*N_cell+2]   = cells_pos_and_mom_array[4*cell_send+2];
                  mpi_2d_send_cells_pos_mom[cpu_recv][4*N_cell+3]   = cells_pos_and_mom_array[4*cell_send+3];
                  //std::cout <<"initiual\t" <<  N_cell << '\t' <<  cpu_recv << "\t" << cell_send << '\t' << cells_num_atoms_in_cell[cell_send] << '\t' << cells_pos_and_mom_array[4*cell_send+0] << '\t'<<std::endl;


                  // save tmp atoms data to be sent
                  for(int j=0; j<cells_num_atoms_in_cell[cell_send]; j++){

                     int N                                  = counter[cpu_recv];
                     int atom_id                            = cells_index_atoms_array[cell_send][j];
                     int type                               = atom_type_array[atom_id];
                     const double mus                       = mp::material[type].mu_s_SI/9.27400915e-24;
                     mpi_2d_send_atoms_cell[cpu_recv][N]    = cell_send;
                     mpi_2d_send_atoms_id[cpu_recv][N]      = atom_id;
                     mpi_2d_send_atoms_pos_x[cpu_recv][N]   = atom_pos_x[atom_id];
                     mpi_2d_send_atoms_pos_y[cpu_recv][N]   = atom_pos_y[atom_id];
                     mpi_2d_send_atoms_pos_z[cpu_recv][N]   = atom_pos_z[atom_id];
                     mpi_2d_send_atoms_mom[cpu_recv][N]     = mus;
                     counter[cpu_recv]++;
                     //std::cout << j << '\t' <<  atom_id << "\t" << atom_pos_x[atom_id] << '\t' << atom_pos_y[atom_id] << '\t' << atom_pos_z[atom_id] << '\t'<<std::endl;

                  }
                  counter_cells[cpu_recv]++;
                  counter_four_cells[cpu_recv] +=4;
                  //counts the number of cells
               }
            }
         }


         //calcualte the displacesments for the gather command
         for (int proc_rec = 1; proc_rec < vmpi::num_processors; proc_rec ++){
             one_displacements[proc_rec] = one_displacements[proc_rec-1] + vmpi::num_processors;
          }

          //calculate the total number of things sent.
         int total = one_displacements[vmpi::num_processors-1] + one_count[vmpi::num_processors - 1];

         //resize the receive arrays
         recv_counter.resize(total);
         recv_counter_cells.resize(total);
         recv_counter_four_cells.resize(total);


        for (int proc_recv = 0; proc_recv < vmpi::num_processors; proc_recv ++){
           MPI_Gatherv(&receive_counts[0],      vmpi::num_processors, MPI_INT, &recv_counter[0],         &one_count[0], &one_displacements[0], MPI_INT, proc_recv, MPI_COMM_WORLD);
           MPI_Gatherv(&receive_counts_cell[0], vmpi::num_processors, MPI_INT, &recv_counter_cells[0],   &one_count[0], &one_displacements[0], MPI_INT, proc_recv, MPI_COMM_WORLD);
        }


         int recv = 0;
         int send = 0;



         for (int i= 0; i < vmpi::num_processors*vmpi::num_processors; i ++){

            if (recv == vmpi::my_rank){
               final_recieve_counter[send]            = recv_counter[i];
               final_recieve_counter_cells[send]      = recv_counter_cells[i];
               final_recieve_counter_four_cells[send] = final_recieve_counter_cells[send]*4;
            }

            recv++;

            if (recv > vmpi::num_processors -1){
             recv = 0;
             send = send +1;
            }
       }




       for (int proc_rec = 1; proc_rec < vmpi::num_processors; proc_rec ++){
          receive_displacements[proc_rec]             = receive_displacements[proc_rec-1]          + final_recieve_counter[proc_rec-1];
          receive_displacements_cells[proc_rec]       = receive_displacements_cells[proc_rec-1]    + final_recieve_counter_cells[proc_rec-1];
          receive_displacements_four_cells[proc_rec]  = receive_displacements_cells[proc_rec]*4;
        }



        //last index+ number in final displacement!
        num_recv_atoms = receive_displacements[vmpi::num_processors -1]       + final_recieve_counter[vmpi::num_processors -1];
        num_recv_cells = receive_displacements_cells[vmpi::num_processors -1] + final_recieve_counter_cells[vmpi::num_processors -1];

        mpi_recv_atoms_id.resize(num_recv_atoms);
        mpi_recv_atoms_pos_x.resize(num_recv_atoms);
        mpi_recv_atoms_pos_y.resize(num_recv_atoms);
        mpi_recv_atoms_pos_z.resize(num_recv_atoms);
        mpi_recv_atoms_mom.resize(num_recv_atoms);
        mpi_recv_atoms_cell.resize(num_recv_atoms);
        mpi_recv_num_atoms_in_cell.resize(num_recv_cells);
        mpi_recv_cells_pos_mom.resize(4*num_recv_cells);





        for (int proc_recv = 0; proc_recv < vmpi::num_processors; proc_recv ++){

           int N        = counter[proc_recv];
           int N_cells  = counter_cells[proc_recv];
           //If N ==0 skip

           mpi_send_atoms_id.resize(N);
           mpi_send_atoms_pos_x.resize(N);
           mpi_send_atoms_pos_y.resize(N);
           mpi_send_atoms_pos_z.resize(N);
           mpi_send_atoms_mom.resize(N);
           mpi_send_atoms_cell.resize(N);
           mpi_send_num_atoms_in_cell.resize(N_cells);
           mpi_send_cells_pos_mom.resize(4*N_cells);


           for (int i = 0; i < N;i++){
             mpi_send_atoms_id[i]               = mpi_2d_send_atoms_id[proc_recv][i];
             mpi_send_atoms_pos_x[i]            = mpi_2d_send_atoms_pos_x[proc_recv][i];
             mpi_send_atoms_pos_y[i]            = mpi_2d_send_atoms_pos_y[proc_recv][i];
             mpi_send_atoms_pos_z[i]            = mpi_2d_send_atoms_pos_z[proc_recv][i];
             mpi_send_atoms_mom[i]              = mpi_2d_send_atoms_mom[proc_recv][i];
             mpi_send_atoms_cell[i]             = mpi_2d_send_atoms_cell[proc_recv][i];
            //std::cout << proc_recv << '\t' << i << '\t' <<  mpi_send_atoms_id[i] << "\t" << mpi_send_atoms_pos_x[i] << '\t' << mpi_send_atoms_pos_y[i] << '\t' << mpi_send_atoms_pos_z[i] << '\t'<<std::endl;
            //std::cout << i << '\t' <<  proc_recv << '\t' << cells_num_atoms_in_cell[cell_send] << '\t' << cells_pos_and_mom_array[4*cell_send+0] << '\t'<<std::endl;

            }
            //std::cout <<"N\t" <<  N_cells << std::endl;
            for (int i = 0; i < N_cells;i++){
               mpi_send_cells_pos_mom[4*i + 0]  = mpi_2d_send_cells_pos_mom[proc_recv][4*i + 0];
               mpi_send_cells_pos_mom[4*i + 1]  = mpi_2d_send_cells_pos_mom[proc_recv][4*i + 1];
               mpi_send_cells_pos_mom[4*i + 2]  = mpi_2d_send_cells_pos_mom[proc_recv][4*i + 2];
               mpi_send_cells_pos_mom[4*i + 3]  = mpi_2d_send_cells_pos_mom[proc_recv][4*i + 3];
               mpi_send_num_atoms_in_cell[i]    = mpi_2d_send_num_atoms_in_cell[proc_recv][i];
            //   std::cout << i << '\t' << proc_recv << '\t' << mpi_send_num_atoms_in_cell[i] << '\t' <<  mpi_send_cells_pos_mom[4*i + 0] << "\t" << mpi_send_cells_pos_mom[4*i + 1] << '\t' << mpi_send_cells_pos_mom[4*i + 2] << '\t' << mpi_send_cells_pos_mom[4*i + 3] << '\t'<<std::endl;
             }

             MPI_Gatherv(&mpi_send_atoms_id[0],          counter[proc_recv],              MPI_INT,    &mpi_recv_atoms_id[0],           &final_recieve_counter[0],             &receive_displacements[0],             MPI_INT,    proc_recv, MPI_COMM_WORLD);
             MPI_Gatherv(&mpi_send_atoms_pos_x[0],       counter[proc_recv],              MPI_DOUBLE, &mpi_recv_atoms_pos_x[0],        &final_recieve_counter[0],             &receive_displacements[0],             MPI_DOUBLE, proc_recv, MPI_COMM_WORLD);
             MPI_Gatherv(&mpi_send_atoms_pos_y[0],       counter[proc_recv],              MPI_DOUBLE, &mpi_recv_atoms_pos_y[0],        &final_recieve_counter[0],             &receive_displacements[0],             MPI_DOUBLE, proc_recv, MPI_COMM_WORLD);
             MPI_Gatherv(&mpi_send_atoms_pos_z[0],       counter[proc_recv],              MPI_DOUBLE, &mpi_recv_atoms_pos_z[0],        &final_recieve_counter[0],             &receive_displacements[0],             MPI_DOUBLE, proc_recv, MPI_COMM_WORLD);
             MPI_Gatherv(&mpi_send_atoms_mom[0],         counter[proc_recv],              MPI_DOUBLE, &mpi_recv_atoms_mom[0],          &final_recieve_counter[0],             &receive_displacements[0],             MPI_DOUBLE, proc_recv, MPI_COMM_WORLD);
             MPI_Gatherv(&mpi_send_atoms_cell[0],        counter[proc_recv],              MPI_INT,    &mpi_recv_atoms_cell[0],         &final_recieve_counter[0],             &receive_displacements[0],             MPI_INT,    proc_recv, MPI_COMM_WORLD);
             MPI_Gatherv(&mpi_send_num_atoms_in_cell[0], counter_cells[proc_recv],        MPI_INT,    &mpi_recv_num_atoms_in_cell[0],  &final_recieve_counter_cells[0],       &receive_displacements_cells[0],       MPI_INT,    proc_recv, MPI_COMM_WORLD);
             MPI_Gatherv(&mpi_send_cells_pos_mom[0],     counter_four_cells[proc_recv],   MPI_DOUBLE, &mpi_recv_cells_pos_mom[0],      &final_recieve_counter_four_cells[0],  &receive_displacements_four_cells[0],  MPI_DOUBLE, proc_recv, MPI_COMM_WORLD);
          }



          // std:: cout<<  "number of received cells:\t" << num_recv_cells <<std::endl;
          // std:: cout<<  "proc cell index 1D:\t" << proc_cell_index_array1D.size() <<std::endl;
          // std:: cout<<  "cells num cells:\t" << cells_num_cells <<std::endl;
         std::vector<int>  recv_cell_id(num_recv_cells);
         std::vector<int>  old_size_array(num_recv_cells,0);
         std::vector<bool> bool_array(num_recv_cells,1); /// bool arrays to check whether a cell has been already considered

         int cell = 0;
         for(int lc=0; lc<num_recv_cells; lc++){
            int cell_recv_counter=0;
            // resize arrays
         //   std::cout <<lc << '\t' <<  mpi_recv_cells_pos_mom[4*lc+0] << '\t' << mpi_recv_cells_pos_mom[4*lc+1] << '\t' << mpi_recv_cells_pos_mom[4*lc+2] << '\t' << mpi_recv_cells_pos_mom[4*lc+3] <<std::endl;
            for(unsigned int i=cells_num_cells; i<proc_cell_index_array1D.size(); i++){
               cell = i;
               int size = cells_index_atoms_array[cell].size();
               old_size_array[lc]=size;
               //std::cout << lc << '\t' << cell << "\t" << cells_num_cells << '\t' << proc_cell_index_array1D.size() << "\t" << mpi_recv_cells_pos_mom[4*lc+0]<< "\t" << cells_pos_and_mom_array[4*cell+0]  << "\t" << mpi_recv_cells_pos_mom[4*lc+1]<< "\t" << cells_pos_and_mom_array[4*cell+1]  << "\t" << mpi_recv_cells_pos_mom[4*lc+2]<< "\t" << cells_pos_and_mom_array[4*cell+2] << "\t" << bool_array[lc] << std::endl;

               if((mpi_recv_cells_pos_mom[4*lc+0] == cells_pos_and_mom_array[4*cell+0]) &&
               (mpi_recv_cells_pos_mom[4*lc+1]    == cells_pos_and_mom_array[4*cell+1]) &&
               (mpi_recv_cells_pos_mom[4*lc+2]    == cells_pos_and_mom_array[4*cell+2]) && bool_array[lc]!=0){
               //   std::cout << "A" << std::endl;
                  recv_cell_id[lc] = cell;
                  //std::cout << lc << '\t' << cell << std::endl;// "\t" << cells_num_cells << '\t' << proc_cell_index_array1D.size() << "\t" << num_recv_cells<< std::endl;


                  cell_recv_counter++;
                  bool_array[lc]=0;
               }
               else if((mpi_recv_cells_pos_mom[4*lc+0] == cells_pos_and_mom_array[4*cell+0]) &&
                       (mpi_recv_cells_pos_mom[4*lc+1] == cells_pos_and_mom_array[4*cell+1]) &&
                       (mpi_recv_cells_pos_mom[4*lc+2] == cells_pos_and_mom_array[4*cell+2]) && bool_array[lc]==0){
                     //     std::cout << "B" << std::endl;

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
            //std::cout <<cell << '\t' <<  lc << "\t" << num_recv_cells << '\t' << mpi_recv_num_atoms_in_cell[lc] << '\t' << cells_num_atoms_in_cell[cell] <<std::endl;

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



         //    } //end else statement
         // }

         //const double is = sizeof(int);
         //const double ds = sizeof(double);

         /*double mem_tot =  double(list_cpu_to_send_to.size())*is +
                           double(list_cells_to_send.size())*is +
                           double(list_cells_to_recv.size())*is +
                           double(mpi_send_atoms_cell.size())*is +
                           double(mpi_send_num_atoms_in_cell.size())*is +
                           double(mpi_send_atoms_id.size())*is +
                           double(mpi_send_atoms_pos_x.size())*ds +
                           double(mpi_send_atoms_pos_y.size())*ds +
                           double(mpi_send_atoms_pos_z.size())*ds +
                           double(mpi_send_atoms_mom.size())*ds +
                           double(mpi_send_cells_pos_mom.size())*ds +
                           double(mpi_recv_atoms_cell.size())*is +
                           double(mpi_recv_num_atoms_in_cell.size())*is +
                           double(mpi_recv_atoms_id.size())*is +
                           double(mpi_recv_atoms_pos_x.size())*ds +
                           double(mpi_recv_atoms_pos_y.size())*ds +
                           double(mpi_recv_atoms_pos_z.size())*ds +
                           double(mpi_recv_atoms_mom.size())*ds +
                           double(mpi_recv_cells_pos_mom.size())*ds;*/

         //double global_tot = 0.0;
         //MPI_Reduce(&mem_tot, &global_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
         //std::cout << "Total memory for tensor construction (all CPUS): " << global_tot*1.0e-6 << " MB" << std::endl;
         //zlog << zTs() << "Total memory for tensor construction (all CPUS): " << global_tot*1.0e-6 << " MB"<< std::endl;

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

         // stop timer
         timer.stop();

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
      int send_cells_demag_factor(std::vector<int>& cells_cell_id_array, // 1D array of cell id's to be sent to root process
                                  std::vector<double>& N_tensor_array,   // 6N tensor of local tensor components
                                  int cells_num_local_cells              // number of local cells on processor
      ){

         // allocate memory to send data
         const int num_send_cells = cells_num_local_cells;                        // number of objects to send
         std::vector<int> mpi_send_cells_id(num_send_cells,0);                    // array storing cell ids to be sent
         std::vector<double> mpi_send_cells_demag_factor(6*num_send_cells,0.0);   // array strong demag factor and self term array to be sent

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

         // arrays of MPI requests
         std::vector<MPI_Request> requests(0);

         MPI_Request req = MPI_REQUEST_NULL; // temporary variable for push_back operations
         MPI_Status status; // temporary variable for stati

         //------------------------------------------------------------
         // send cells id, demag factors and self term from all CPUs
         //------------------------------------------------------------
         requests.push_back(req);
         MPI_Isend(&num_send_cells, 1, MPI_INT, 0, 120, MPI_COMM_WORLD, &requests.back());
         requests.push_back(req);
         MPI_Isend(&mpi_send_cells_id[0], num_send_cells, MPI_INT, 0, 121, MPI_COMM_WORLD, &requests.back());
         requests.push_back(req);
         MPI_Isend(&mpi_send_cells_demag_factor[0], 6*num_send_cells, MPI_DOUBLE, 0, 122, MPI_COMM_WORLD, &requests.back());

         // loop over CPUs
         for(int cpu=0; cpu<vmpi::num_processors; cpu++){
            // if I am root proc, then I receive cells demag factors
            if(vmpi::my_rank == 0){
               // allocate int to receive num of cells to be received
               int num_recv_cells;
               // Receive num_recv_cells
               requests.push_back(req);
               MPI_Irecv(&num_recv_cells, 1, MPI_INT, cpu, 120, MPI_COMM_WORLD, &requests.back());
               MPI_Wait(&requests.back(), &status); // wait for number of cells to receive
               // Allocate arrays for demag factor
               std::vector<int> mpi_recv_cells_id(num_recv_cells,0);
               std::vector<double> mpi_recv_cells_demag_factor(6*num_recv_cells,0.0);
               // Receive arrays
               requests.push_back(req);
               MPI_Irecv(&mpi_recv_cells_id[0], num_recv_cells, MPI_INT, cpu, 121, MPI_COMM_WORLD, &requests.back());
               MPI_Wait(&requests.back(), &status); // wait for data to be received
               requests.push_back(req);
               MPI_Irecv(&mpi_recv_cells_demag_factor[0], 6*num_recv_cells, MPI_DOUBLE, cpu, 122, MPI_COMM_WORLD, &requests.back());
               MPI_Wait(&requests.back(), &status); // wait for data to be received

               // Save received data (only once for each cell, duplicates are discarded by being overwritten)
               for(int i=0; i<num_recv_cells; i++){

                  // get cell id to save demag tensor into
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

         // wait for all processors (and root) at this point
         vmpi::barrier();

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

      // array of MPI requests
      std::vector<MPI_Request> requests(0);

      MPI_Request req = MPI_REQUEST_NULL; // temporary variable for push_back operations
      MPI_Status status; // temporary variable for stati

      // send cells id, B-field, Hd-field
      requests.push_back(req);
      MPI_Isend(&num_send_cells, 1, MPI_INT, 0, 114, MPI_COMM_WORLD, &requests.back());
      requests.push_back(req);
      MPI_Isend(&mpi_send_cells_id[0], num_send_cells, MPI_INT, 0, 115, MPI_COMM_WORLD, &requests.back());
      requests.push_back(req);
      MPI_Isend(&mpi_send_cells_field[0], 3*num_send_cells, MPI_DOUBLE, 0, 116, MPI_COMM_WORLD, &requests.back());

      // loop over CPUs
      for(int cpu=0; cpu<vmpi::num_processors; cpu++){
         // if I am root proc, then I receive cells dipolar field
         if(vmpi::my_rank == 0){
            // allocate int to receive num of cells to be received
            int num_recv_cells;
            // Receive num_recv_cells
            requests.push_back(req);
            MPI_Irecv(&num_recv_cells, 1, MPI_INT, cpu, 114, MPI_COMM_WORLD, &requests.back());
            MPI_Wait(&requests.back(), &status); // wait for number of data to be received
            // Allocate arrays for field
            std::vector<int> mpi_recv_cells_id(num_recv_cells,0);
            std::vector<double> mpi_recv_cells_field(3*num_recv_cells,0.0);
            // Receive arrays
            requests.push_back(req);
            MPI_Irecv(&mpi_recv_cells_id[0], num_recv_cells, MPI_INT, cpu, 115, MPI_COMM_WORLD, &requests.back());
            MPI_Wait(&requests.back(), &status); // wait for data to be received
            requests.push_back(req);
            MPI_Irecv(&mpi_recv_cells_field[0], 3*num_recv_cells, MPI_DOUBLE, cpu, 116, MPI_COMM_WORLD, &requests.back());
            MPI_Wait(&requests.back(), &status); // wait for data to be received
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

      // wait for all processes here before de-allocating memory
      vmpi::barrier();

      // Clear arrays used only for sending data
      mpi_send_cells_id.clear();
      mpi_send_cells_field.clear();

      return EXIT_SUCCESS;
   }


   #endif

} // end namespace dipole
