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
#include "hamr.hpp"
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

#ifdef MPICF
//------------------------------------------------------------------------------
// Function to write a window snapshot of an io group (one file per node): the
// ranges and spins of every process in the group are gathered on the io master,
// which offsets each process's positions by the position of its first atom in
// the group's file and writes them in process order
//------------------------------------------------------------------------------
double write_window_data_group(const std::string& filename, const std::vector<uint64_t>& ranges, const std::vector<double>& buffer){

   // number of range values and spin values of each process
   const int sizes[2] = {int(ranges.size()), int(buffer.size())};
   std::vector<int> all_sizes(io_group_master ? 2 * io_group_size : 0);
   MPI_Gather(sizes, 2, MPI_INT, all_sizes.data(), 2, MPI_INT, io_group_master_id, io_comm);

   std::vector<int> range_counts, range_displs, data_counts, data_displs;
   std::vector<uint64_t> all_ranges;
   std::vector<double> all_data;
   if(io_group_master){
      int range_total = 0, data_total = 0;
      for(int p = 0; p < io_group_size; p++){
         range_counts.push_back(all_sizes[2*p]);
         range_displs.push_back(range_total);
         range_total += all_sizes[2*p];
         data_counts.push_back(all_sizes[2*p+1]);
         data_displs.push_back(data_total);
         data_total += all_sizes[2*p+1];
      }
      all_ranges.resize(range_total);
      all_data.resize(data_total);
   }
   MPI_Gatherv(ranges.data(), sizes[0], MPI_UINT64_T, all_ranges.data(), range_counts.data(), range_displs.data(), MPI_UINT64_T, io_group_master_id, io_comm);
   MPI_Gatherv(buffer.data(), sizes[1], MPI_DOUBLE, all_data.data(), data_counts.data(), data_displs.data(), MPI_DOUBLE, io_group_master_id, io_comm);

   if(!io_group_master) return 0.0;

   // the atoms of each process follow those of the processes before it in the file (displacements count x, y, z)
   for(int p = 0; p < io_group_size; p++){
      for(int i = 0; i < range_counts[p]; i += 2) all_ranges[range_displs[p] + i] += io_group_displacements[p] / 3;
   }
   return write_window_data(filename, collated_buffer.size() / 3, all_ranges, all_data);

}
#endif

//------------------------------------------------------------------------------
// Function to output atomic spin positions to disk
//------------------------------------------------------------------------------
void atoms(const bool allow_window){

   // instantiate timer for total output time including overheads
   vutil::vtimer_t total_timer;

   // start timer
   total_timer.start();

   //------------------------------------------
   // With the HAMR active window only the atoms of grains integrated since the
   // previous snapshot can have changed: write just those (window snapshot),
   // in binary format, one file per node or per process (the same decision on
   // every process)
   //------------------------------------------
   bool window = false;
   if(allow_window && hamr::window_snapshot_possible()){
      if(config::internal::mode != legacy && config::internal::mode != mpi_io && config::internal::format == binary) window = true;
      else{
         static bool logged = false;
         if(!logged) zlog << zTs() << "Spin snapshots are written in full with hamr:active-window: window snapshots need the binary format "
                          << "and one file per node or per process" << std::endl;
         logged = true;
      }
   }

   //------------------------------------------
   // Output Meta Data
   //------------------------------------------
   // get system magnetization data
   const std::vector<double> magnetisation = stats::system_magnetization.get_magnetization();
   // calculate real time
   const double real_time = double(sim::time) * mp::dt_SI;

   if(config::internal::mode != legacy && vmpi::my_rank == 0){
      write_meta(real_time, sim::temperature, sim::H_vec[0], sim::H_vec[1], sim::H_vec[2], sim::H_applied, magnetisation[0], magnetisation[1], magnetisation[2], window);
   }

   //------------------------------------------
   // Output spin data
   //------------------------------------------

   // copy data to local buffer: all output atoms, or for a window snapshot the changed ones as
   // (first position, number) ranges of positions in the output order
   std::vector<uint64_t> window_ranges;
   std::vector<double> window_buffer;
   if(window){
      for(uint64_t p = 0; p < local_output_atom_list.size(); p++){
         const uint64_t atom = local_output_atom_list[p];
         if(!hamr::spin_changed_since_snapshot(atom)) continue;
         if(!window_ranges.empty() && window_ranges[window_ranges.size() - 2] + window_ranges.back() == p) window_ranges.back()++;
         else{
            window_ranges.push_back(p);
            window_ranges.push_back(1);
         }
         window_buffer.push_back(atoms::x_spin_array[atom]);
         window_buffer.push_back(atoms::y_spin_array[atom]);
         window_buffer.push_back(atoms::z_spin_array[atom]);
      }
   }
   else copy_data_to_buffer(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, local_output_atom_list, config::internal::local_buffer);

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
         if(window) io_time = write_window_data(filename, local_output_atom_list.size(), window_ranges, window_buffer);
         else io_time = write_data(filename, config::internal::local_buffer);
         break;

      case config::internal::fpnode:
         if(window) io_time = write_window_data_group(filename, window_ranges, window_buffer);
         else{
            // Gather data from all processors in io group
            MPI_Gatherv(&local_buffer[0], local_buffer.size(), MPI_DOUBLE, &collated_buffer[0], &io_group_recv_counts[0], &io_group_displacements[0], MPI_DOUBLE, io_group_master_id, io_comm);
            // output data on master io processes
            if(config::internal::io_group_master) io_time = write_data(filename, config::internal::collated_buffer);
         }
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
      // window snapshot
      if(window) io_time = write_window_data(filename, local_output_atom_list.size(), window_ranges, window_buffer);
      // check for legacy output
      else if(config::internal::mode == config::internal::legacy) io_time = config::internal::legacy_atoms();
      // otherwise use new one by default
      else io_time = write_data(filename, config::internal::local_buffer);
   #endif

   // stop total timer
   total_timer.stop();

   // Output bandwidth to log file
   const double data_size = window ? 1.0e-9 * double(sizeof(double)) * double(window_buffer.size()) : config::internal::io_data_size;
   zlog << data_size/io_time << " GB/s in " << io_time << " s [ " << total_timer.elapsed_time() << " s]";
   if(window) zlog << " (window snapshot, " << window_buffer.size() / 3 << " atoms in " << window_ranges.size() / 2 << " ranges)";
   zlog << std::endl;

   // the spins of the next snapshot are compared with this one (HAMR active window)
   hamr::snapshot_written();

   // increment file counter
   sim::output_atoms_file_counter++;

   return;

}

} // end of internal namespace
} // end of config namespace
