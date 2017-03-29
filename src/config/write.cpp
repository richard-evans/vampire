//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------
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
// Rory these functions need fixing
namespace config
{
namespace internal
{
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
void copy_data_to_buffer(const std::vector<double> &x, // vector data
                         const std::vector<double> &y,
                         const std::vector<double> &z,
                         const std::vector<int> &mask,
                         std::vector<float> &buffer)
{

   // copy total number of output data to const for compiler
   const int data_size = mask.size();

   // loop over all atoms to be output
   for (int id = 0; id < data_size; ++id)
   {
      // determine next datum to be output
      const int index = mask[id];
      // copy and cast data to be output to main output buffer
      buffer[3 * id + 0] = float(x[index]);
      buffer[3 * id + 1] = float(y[index]);
      buffer[3 * id + 2] = float(z[index]);
   }

   return;
}

//----------------------------------------------------------------------------------------------------
// Simple wrapper function to call output function for correct format
//----------------------------------------------------------------------------------------------------
//
void write_data(const std::vector<float> &buffer, bool coord)
{
      printf("Writing data");
   std::string filename = data_filename(coord);
   // Output informative message to log file
   zlog << zTs() << "Outputting configuration file " << filename << " to disk ";

   switch (config::internal::output_data_format)
   {

   case config::internal::binary:
      printf("Writing binary");
      write_data_binary(filename, buffer);
      break;
   case config::internal::text:
      printf("Writing text");
      write_data_text(filename, buffer);
      break;
   }

   // increment file counter
   sim::output_atoms_file_counter++;

   return;
}

//----------------------------------------------------------------------------------------------------
// Function to output spin data formatted as text
//----------------------------------------------------------------------------------------------------
//
void write_data_text(std::string filename, const std::vector<float> &buffer)
{

   // Output informative message to log file
   //zlog << zTs() << "Outputting configuration file " << filename.str() << " to disk ";

   // Declare and open output file
   std::ofstream ofile;
   ofile.open(filename.c_str());

   // determine number of data to output
   const int buffer_size = buffer.size() / 3;

   // output number of data
   ofile << buffer_size << "\n";

   // output buffer to disk
   for (int index = 0; index < buffer_size; ++index)
   {
      ofile << buffer[3 * index + 0] << "\t"
            << buffer[3 * index + 1] << "\t"
            << buffer[3 * index + 2] << "\n";
   }

   // close output file
   ofile.close();

   // open file at end
   std::ifstream in(filename.c_str(), std::ios::binary | std::ios::ate);

   // get file size (bytes)
   double local_data_size = double(in.tellg());

   // close file
   in.close();

   return;
}

//----------------------------------------------------------------------------------------------------
// Function to output spin data in binary format
//----------------------------------------------------------------------------------------------------
//
void write_data_binary(std::string filename, const std::vector<float> &buffer)
{
   // Declare and open output file
   std::ofstream ofile;
   ofile.open(filename.c_str(), std::ios::binary);

   // determine number of data to output
   const int buffer_size = buffer.size() / 3;

   // output number of data
   ofile.write(reinterpret_cast<const char *>(&buffer_size), sizeof(unsigned int));

   // output buffer to disk
   ofile.write(reinterpret_cast<const char *>(&buffer[0]), sizeof(float) * buffer.size());

   // close output file
   ofile.close();

   // get file size (bytes)
   double local_data_size = double(sizeof(float) * buffer.size());

   return;
}

std::string data_filename(bool coords){
   // Set local output filename
   std::stringstream file_sstr;

   if (coords)
      file_sstr << "atom-coords-";
   else
      file_sstr << "atom-spins-";

   switch (config::internal::output_data_format)
   {

   case config::internal::binary:
      file_sstr << "binary-";
      break;
   case config::internal::text:
      file_sstr << "text-";
      break;
   }
   #ifdef MPICF
      if (!mpi_io)
         file_sstr << std::setfill('0') << std::setw(5) << vmpi::my_io_group;
   #endif
   if (!coords)
      file_sstr << "-" << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter;
   file_sstr << ".data";
   std::string filename = file_sstr.str();

   return filename;
}

} // end of namespace internal
} // end of namespace config
