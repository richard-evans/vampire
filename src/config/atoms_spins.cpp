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

//------------------------------------------------------------------------------
// Function to output atomic spin positions to disk
//------------------------------------------------------------------------------
void atoms(){

   // instantiate timer for total output time including overheads
   vutil::vtimer_t total_timer;

   // start timer
   total_timer.start();

   //------------------------------------------
   // Output Meta Data
   //------------------------------------------
   // get system magnetization data
   const std::vector<double> magnetisation = stats::system_magnetization.get_magnetization();
   // calculate real time
   const double real_time = double(sim::time) * mp::dt_SI;

   if(config::internal::mode != legacy && vmpi::my_rank == 0){
      write_meta(real_time, sim::temperature, sim::H_vec[0], sim::H_vec[1], sim::H_vec[2], sim::H_applied, magnetisation[0], magnetisation[1], magnetisation[2]);
   }

   //------------------------------------------
   // Output spin data
   //------------------------------------------

   // copy data to local buffer
   copy_data_to_buffer(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, local_output_atom_list, config::internal::local_buffer);

   // Determine output filename
   std::stringstream file_sstr;

   // set simple file name for single file output
   if(config::internal::num_io_groups == 1) file_sstr << "spins-" << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter << ".data";
   // otherwise set indexed files
   else{
      file_sstr << "spins-" << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter << "-" <<
                               std::setfill('0') << std::setw(6) << config::internal::io_group_id << ".data";
   }

   // convert stringstream to string
   std::string filename = file_sstr.str();

   // Output informative message to log file on root process
   zlog << zTs() << "Outputting configuration file " << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter << " to disk " << std::flush;

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
         io_time = legacy_atoms();
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

         // Calculate local byte offset since MPI-IO is simple and doesn't update the file handle pointer after I/O
         MPI_Offset data_offset = config::internal::buffer_offset + sizeof(uint64_t);

         timer.start(); // start timer

         // Write data to disk
         MPI_File_write_at_all(fh, data_offset, &config::internal::local_buffer[0], config::internal::local_buffer.size(), MPI_DOUBLE, &status);
         //MPI_File_write_ordered(fh, &config::internal::local_buffer[0], config::internal::local_buffer.size(), MPI_DOUBLE, &status);

         timer.stop(); // Stop timer

         // Close file
         MPI_File_close(&fh);

         // Calculate elapsed time
         io_time = timer.elapsed_time();
         break;
      }

      case config::internal::fpprocess:
         io_time = write_data(filename, config::internal::local_buffer);
         break;

      case config::internal::fpnode:
         // Gather data from all processors in io group
         MPI_Gatherv(&local_buffer[0], local_buffer.size(), MPI_DOUBLE, &collated_buffer[0], &io_group_recv_counts[0], &io_group_displacements[0], MPI_DOUBLE, io_group_master_id, io_comm);
         // output data on master io processes
         if(config::internal::io_group_master) io_time = write_data(filename, config::internal::collated_buffer);
         double max_io_time = 0.0;
         // calculate actual bandwidth on root process
         MPI_Reduce(&io_time, &max_io_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
         io_time = max_io_time;
         break;

   }

   #else
      //-----------------------------------------------------
      // Serial mode output (ignores most io directives)
      //-----------------------------------------------------
      // check for legacy output
      if(config::internal::mode == config::internal::legacy) io_time = config::internal::legacy_atoms();
      // otherwise use new one by default
      else io_time = write_data(filename, config::internal::local_buffer);
   #endif

   // stop total timer
   total_timer.stop();

   // Output bandwidth to log file
   zlog << config::internal::io_data_size/io_time << " GB/s in " << io_time << " s [ " << total_timer.elapsed_time() << " s]" << std::endl;

   // increment file counter
   sim::output_atoms_file_counter++;

   return;

}

} // end of internal namespace
} // end of config namespace
