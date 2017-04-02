//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Rory Pond and Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

//#include <stdio.h>

// C++ standard library headers
//#include <cmath>
//#include <cstdlib>
//#include <ctime>
#include <iomanip>
//#include <iostream>
//#include <fstream>
#include <sstream>
//#include <string>

// Vampire headers
#include "atoms.hpp"
#include "cells.hpp"
#include "config.hpp"
//#include "dipole.hpp"
//#include "errors.hpp"
//#include "LLG.hpp"
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

   //------------------------------------------
   // Output Meta Data
   //------------------------------------------
   // get system magnetization data
   const std::vector<double> magnetisation = stats::system_magnetization.get_magnetization();
   // calculate real time
   const double real_time = double(sim::time) * mp::dt_SI;
   // set number of files
   const int files = config::internal::num_io_groups;

   if(config::internal::mode != legacy && vmpi::my_rank == 0){
      write_meta(real_time, sim::temperature, sim::H_vec[0], sim::H_vec[1], sim::H_vec[2], magnetisation[0], magnetisation[1], magnetisation[2], files);
   }

   //------------------------------------------
   // Output spin data
   //------------------------------------------

   // copy data to local buffer
   copy_data_to_buffer(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, local_output_atom_list, config::internal::local_buffer);

   // Determine output filename
   std::stringstream file_sstr;
   file_sstr << "atoms-spins-" << std::setfill('0') << std::setw(5) << config::internal::io_group_id;
   file_sstr << "-" << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter;
   file_sstr << ".data";
   std::string filename = file_sstr.str();

   // Output informative message to log file on root process
   zlog << zTs() << "Outputting configuration file " << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter << " to disk ";

   double bandwidth = 0.0;

   //-----------------------------------------------------
   // Parallel mode output
   //-----------------------------------------------------
   #ifdef MPICF

   // Determine io mode and call appropriate function for data
   switch(config::internal::mode){

      // legacy
      case config::internal::legacy:
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
         // Wait for Everyone
         vmpi::barrier();
         timer.start(); // start timer
         // Write data to disk
         MPI_File_write_ordered(fh, &config::internal::local_buffer[0], config::internal::local_buffer.size(), MPI_DOUBLE, &status);
         // Close file
         MPI_File_close(&fh);
         timer.stop(); // Stop timer
         bandwidth = config::internal::io_data_size / timer.elapsed_time();
         break;
      }

      case config::internal::fpprocess:
         bandwidth = write_data(filename, config::internal::local_buffer);
         break;

      case config::internal::fpnode:
         // Gather data from all processors in io group
         MPI_Gatherv(&local_buffer[0], local_buffer.size(), MPI_DOUBLE, &collated_buffer[0], &io_group_recv_counts[0], &io_group_displacements[0], MPI_DOUBLE, io_group_master_id, io_comm);
         // output data on master io processes
         if(config::internal::io_group_master) bandwidth = write_data(filename, config::internal::collated_buffer);
         // correct calculation for bandwidth (default calculation assumes all atoms are output to disk on each process)
         double total_bandwidth = 0.0;
         if(!io_group_master) bandwidth = 1e6; // set unachievable bandwidth for all processes not outputting data
         // calculate actual bandwidth on root process
         MPI_Reduce(&bandwidth, &total_bandwidth, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
         bandwidth = total_bandwidth / double(num_io_groups);
         break;

   }

   #else
      //-----------------------------------------------------
      // Serial mode output (ignores most io directives)
      //-----------------------------------------------------
      // check for legacy output
      if(config::internal::mode == config::internal::legacy) config::internal::legacy_atoms();
      // otherwise use new one by default
      else bandwidth = write_data(filename, config::internal::local_buffer);
   #endif

   // Output bandwidth to log file
   zlog << bandwidth << " GB/s" << std::endl;

   // increment file counter
   sim::output_atoms_file_counter++;

   return;

}


} // end of internal namespace
} // end of config namespace
