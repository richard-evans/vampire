//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) Rory Pond and Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
#include <stdio.h>
// C++ standard library headers
#include <iomanip>
#include <sstream>

// Vampire headers
#include "errors.hpp"
#include "config.hpp"
#include "vio.hpp"
#include "vutil.hpp"
#include "sim.hpp"

// config headers
#include "internal.hpp"

namespace config
{
namespace internal
{

// Forward function declarations
double write_data_text(std::string filename, const std::vector<double> &buffer);
double write_data_text_sld(std::string filename, const std::vector<double> &buffer_sld);
double write_data_binary(std::string filename, const std::vector<double> &buffer);

//--------------------------------------------------------------------------------------------------------
//  Function to copy and cast masked 3-vector data array to output buffer (serial and parallel versions)
//
//  The output buffer stores the data in
//
//                       | x y z | x y z | x y z | ... | x y z |
//
//  format which is then written to disk sequentially in binary or text mode.
//
//  Data which are to be output are predetermined in the mask for improved peformance. Data are also
//  cast to float to reduce storage requirements from 24 to 12 bytes per datum for improved write
//  performance and file size.
//
//  Parallel (MPI) mode
//  ------------------------
//  In parallel mode a temporary buffer is required to store the data on each process which is then
//  merged into the output buffer on the master io process, determined from the MPI_COMM_IO
//  communicator, doubling the memory requirement for data i/o.
//
//  output_buffer io_master | x y z | x y z | x y z | x y z | ... | x y z | x y z | x y z | x y z |
//
//                              ^       ^       ^       ^             ^       ^       ^       ^
//                              :       :       :       :             :       :       :       :
//                                                                    :       :       :       :
//  mpi_buffer process_1    | x y z | x y z | x y z | x y z |         :       :       :       :
//
//  mpi_buffer process_n                                          | x y z | x y z | x y z | x y z |
//
//  Best performance is likely achieved for 1 or 2 output processes/node.
//
//--------------------------------------------------------------------------------------------------------
//
// Data is imported as 3 x 1D vectors for x,y and z respectively. The mask identifies which data should
// be outputted as a sparse list (each mask id lists an array index of data to be output. The float buffer
// stores the final complete data to be output to disk and must be 3*mask.size().
//


//----------------------------------------------------------------------------------------------------
// Simple wrapper function to call output function for correct format
//----------------------------------------------------------------------------------------------------
//
double write_data(std::string filename, const std::vector<double> &buffer){

   double io_time = 0.0;

   switch (config::internal::format){

      case config::internal::binary:
         io_time = write_data_binary(filename, buffer);
         break;

      case config::internal::text:
         io_time = write_data_text(filename, buffer);

      //case config::internal::text_sld:
      // io_time = write_data_text_sld(filename, buffer);

         break;

   }

   return io_time;

}

//----------------------------------------------------------------------------------------------------
// Function to output spin data formatted as text
//----------------------------------------------------------------------------------------------------
//
double write_data_text(std::string filename, const std::vector<double> &buffer){

   // Declare and open output file
   std::ofstream ofile;
   ofile.open(filename.c_str());

   // determine number of data to output
   const uint64_t data_size = buffer.size() / 3;

   // output number of data
   ofile << data_size << "\n";

   // instantiate timer
   vutil::vtimer_t timer;

   // start timer
   timer.start();

   // output buffer to disk
   for(unsigned int index = 0; index < data_size; ++index){

      // add to anable high recision output: std::setprecision(17)
      ofile << buffer[3 * index + 0] << "\t"
            << buffer[3 * index + 1] << "\t"
            << buffer[3 * index + 2] << "\n";

   }

   // end timer
   timer.stop();

   // close output file
   ofile.close();

   // return bandwidth
   return timer.elapsed_time();

}

//----------------------------------------------------------------------------------------------------
// Function to output spin data in binary format
//----------------------------------------------------------------------------------------------------
//
double write_data_binary(std::string filename, const std::vector<double> &buffer)
{
   // Declare and open output file
   std::ofstream ofile;
   ofile.open(filename.c_str(), std::ios::binary);

   // determine number of data to output
   const uint64_t data_size = buffer.size() / 3;

   // instantiate timer
   vutil::vtimer_t timer;

   // output number of data
   ofile.write(reinterpret_cast<const char *>(&data_size), sizeof(uint64_t));

   // start timer
   timer.start();

   // output buffer to disk
   ofile.write(reinterpret_cast<const char *>(&buffer[0]), sizeof(double) * buffer.size());

   // end timer
   timer.stop();

   // close output file
   ofile.close();

   // return bandwidth
   return timer.elapsed_time();

}

//----------------------------------------------------------------------------------------------------
// Function to write a window snapshot (active HAMR window) in binary format: includes the spins of the atoms
// that changed since the previous snapshot, as ranges of positions in the order of the coordinate file
//
//    uint64                      number of atoms in the coordinate file
//    uint64                      number of ranges r
//    r x (uint64, uint64)        first position and number of atoms of each range
//    3 x (sum of numbers) double x, y, z of the atoms of the ranges, in order
//
//----------------------------------------------------------------------------------------------------
//
double write_window_data(std::string filename, const uint64_t num_atoms, const std::vector<uint64_t>& ranges, const std::vector<double>& buffer){

   // Declare and open output file
   std::ofstream ofile;
   ofile.open(filename.c_str(), std::ios::binary);

   // instantiate timer
   vutil::vtimer_t timer;

   // start timer
   timer.start();

   const uint64_t num_ranges = ranges.size() / 2;
   ofile.write(reinterpret_cast<const char *>(&num_atoms), sizeof(uint64_t));
   ofile.write(reinterpret_cast<const char *>(&num_ranges), sizeof(uint64_t));
   ofile.write(reinterpret_cast<const char *>(ranges.data()), sizeof(uint64_t) * ranges.size());
   ofile.write(reinterpret_cast<const char *>(buffer.data()), sizeof(double) * buffer.size());

   // end timer
   timer.stop();

   // close output file
   ofile.close();

   // return bandwidth
   return timer.elapsed_time();

}

} // end of namespace internal
} // end of namespace config
