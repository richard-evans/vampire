//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Rory Pond and Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iomanip>
#include <sstream>


// Vampire headers
#include "atoms.hpp"
#include "cells.hpp"
#include "config.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"

// config module headers
#include "internal.hpp"

namespace config{

namespace internal{

// Forward function declaration
void collate_int_data(std::vector<int>& source, std::vector<int>& collated);

//------------------------------------------------------------------------------
// Function to output atomic positions to disk
//------------------------------------------------------------------------------
void atoms_coords()
{

      //------------------------------------------
      // Output Meta Data from root process
      //------------------------------------------
      // set number of files
      // const int files = config::internal::num_io_groups; //unused variable

      if(config::internal::mode != legacy && vmpi::my_rank == 0){
         config::internal::write_coordinate_meta();
      }

      //------------------------------------------------
      // Create temporary buffers for atom information
      //------------------------------------------------
      uint64_t local_atoms = local_output_atom_list.size();

      std::vector<int> atom_type_buffer(local_atoms);
      for(unsigned int atom = 0; atom < local_atoms; atom++) atom_type_buffer[atom] = atoms::type_array[local_output_atom_list[atom]];

      std::vector<int> atom_category_buffer(local_atoms);
      for(unsigned int atom = 0; atom < local_atoms; atom++) atom_category_buffer[atom] = atoms::category_array[local_output_atom_list[atom]];

      //------------------------------------------
      // Output coordinate data
      //------------------------------------------

      // copy coordinate data to local buffer
      copy_data_to_buffer(atoms::x_coord_array, atoms::y_coord_array, atoms::z_coord_array, local_output_atom_list, config::internal::local_buffer);

      // Determine output filename
      std::stringstream file_sstr;

      // set simple file name for single file output
      if(config::internal::num_io_groups == 1) file_sstr << "atoms-coords.data";
      // otherwise set indexed files
      else file_sstr << "atoms-coords-" << std::setfill('0') << std::setw(6) << config::internal::io_group_id << ".data";

      // convert string stream to string
      std::string filename = file_sstr.str();

      // Calculate number of bytes to be written to disk
      const double spin_data_size = double(config::internal::total_output_atoms) * 1.0e-9 * (3.0*double(sizeof(double)) );
      const double coord_data_size = double(config::internal::total_output_atoms) * 1.0e-9 * (3.0*double(sizeof(double) + 2.0*double(sizeof(int)) ) );

      // Output informative messages of actual data size to be outputed to disk (in binary mode)
      zlog << zTs() << "Total coordinate data filesize: " << 1000.0 * coord_data_size << " MB" << std::endl;
      zlog << zTs() << "Total spin data filesize (per snapshot): " << 1000.0 * spin_data_size << " MB" << std::endl;

      // Output informative message to log file on root process
      zlog << zTs() << "Outputting atomic coordinate file to disk ";

      // Variable for calculating output bandwidth
      double io_time = 1.0e-12;

      //-----------------------------------------------------
      // Parallel mode output
      //-----------------------------------------------------
      #ifdef MPICF

      // Determine io mode and call appropriate function for data
      switch(config::internal::mode){

         // legacy
         case config::internal::legacy:
            io_time = config::internal::legacy_atoms_coords();
            break;

         case config::internal::mpi_io:{
            vutil::vtimer_t timer; // instantiate timer
            MPI_File fh; // MPI file handle
            MPI_Status status; // MPI io status
            // convert filename to character string for output
            char *cfilename = (char*)filename.c_str();
            // Open file on all processors
            MPI_File_open(MPI_COMM_WORLD, cfilename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
            // write number of atoms on root process
            if(vmpi::my_rank == 0) MPI_File_write(fh, &total_output_atoms, 1, MPI_UINT64_T, &status);

            // Calculate local byte offsets since MPI-IO is simple and doesn't update the file handle pointer after I/O
            MPI_Offset type_offset = config::internal::linear_offset + sizeof(uint64_t);
            MPI_Offset category_offset = config::internal::linear_offset + total_output_atoms * sizeof(int) + sizeof(uint64_t);
            MPI_Offset data_offset = config::internal::buffer_offset + 2 * total_output_atoms * sizeof(int) + sizeof(uint64_t);

            timer.start(); // start timer

            // Write data to disk
            MPI_File_write_at_all(fh, type_offset, &atom_type_buffer[0], atom_type_buffer.size(), MPI_INT, &status);
            MPI_File_write_at_all(fh, category_offset, &atom_category_buffer[0], atom_category_buffer.size(), MPI_INT, &status);
            MPI_File_write_at_all(fh, data_offset, &config::internal::local_buffer[0], config::internal::local_buffer.size(), MPI_DOUBLE, &status);

            timer.stop(); // Stop timer

            // Calculate elapsed time
            io_time = timer.elapsed_time();

            // Close file
            MPI_File_close(&fh);
            break;
         }

         case config::internal::fpprocess:
            io_time = write_coord_data(filename, config::internal::local_buffer, atom_type_buffer, atom_category_buffer);
            break;

         case config::internal::fpnode:{
            // Gather data from all processors in io group
            std::vector<int> collated_atom_type_buffer(0);
            collate_int_data(atom_type_buffer, collated_atom_type_buffer);
            std::vector<int> collated_atom_category_buffer(0);
            collate_int_data(atom_category_buffer, collated_atom_category_buffer);
            // Gather coord data from all processors in io group
            MPI_Gatherv(&local_buffer[0], local_buffer.size(), MPI_DOUBLE, &collated_buffer[0], &io_group_recv_counts[0], &io_group_displacements[0], MPI_DOUBLE, io_group_master_id, io_comm);

            // output data on master io processes
            if(config::internal::io_group_master) io_time = write_coord_data(filename, config::internal::collated_buffer, collated_atom_type_buffer, collated_atom_category_buffer);
            // find longest time in all io nodes
            double max_io_time = 0.0;
            // calculate actual bandwidth on root process
            MPI_Reduce(&io_time, &max_io_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            io_time = max_io_time;
            break;
         }

      }

      #else
         //-----------------------------------------------------
         // Serial mode output (ignores most io directives)
         //-----------------------------------------------------
         // check for legacy output
         if(config::internal::mode == config::internal::legacy) io_time = config::internal::legacy_atoms_coords();
         // otherwise use new one by default
         else io_time = write_coord_data(filename, config::internal::local_buffer, atom_type_buffer, atom_category_buffer);
      #endif

      // Output bandwidth to log file
      zlog << coord_data_size/io_time << " GB/s in " << io_time << " s" << std::endl;

      return;

   }

   //---------------------------------------------------------------------------
   // Function to collect data from local buffers in io communicator
   //---------------------------------------------------------------------------
   void collate_int_data(std::vector<int>& source, std::vector<int>& collated){

      #ifdef MPICF

      std::vector<int> recv_counts;
      std::vector<int> displacements;

      // Reduce number of local atoms to be outputted on master io process
      uint64_t source_size = source.size();
      uint64_t collated_size = 0;

      MPI_Reduce(&source_size, &collated_size, 1, MPI_UINT64_T, MPI_SUM, config::internal::io_group_master_id, config::internal::io_comm);

      // resize output array, recv_counts and displacements on io master only
      if(io_group_master){
         collated.resize(collated_size);
         recv_counts.resize(io_group_size,0);
         displacements.resize(io_group_size,0);
      }

      // get local data size (int required for MPI Gatherv call)
      int num_data = source.size();

      // gather number of data to be received from each processor in io group
      MPI_Gather(&num_data, 1, MPI_INT, &recv_counts[0], 1, MPI_INT, io_group_master_id, io_comm);

      // calculate displacements for gatherv and total number of points
      if(io_group_master){
         int disp = 0; // initial displacement
         for(int p=0; p<io_group_size; p++){
            displacements[p] = disp;
            disp += recv_counts[p]; // add number of counts to be recieved
         }
      }

      // Now collect data on master io process
      MPI_Gatherv(&source[0], source.size(), MPI_INT, &collated[0], &recv_counts[0], &displacements[0], MPI_INT, io_group_master_id, io_comm);

      #endif

      return;

   }

} // end of internal namespace
} // end of config namespace
