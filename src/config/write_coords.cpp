//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
//   (c) Rory Pond and Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iomanip>
#include <sstream>

// Vampire headers
#include "config.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vutil.hpp"

// config headers
#include "internal.hpp"

namespace config
{
namespace internal
{

// Forward function declarations
double write_coord_data_text(std::string filename, const std::vector<double> &buffer, const std::vector<int> &type_buffer, const std::vector<int> &category_buffer);
double write_coord_data_binary(std::string filename, const std::vector<double> &buffer, const std::vector<int> &type_buffer, const std::vector<int> &category_buffer);

//----------------------------------------------------------------------------------------------------
// Simple wrapper function to call output function for correct format
//----------------------------------------------------------------------------------------------------
//
double write_coord_data(std::string filename, const std::vector<double>& buffer, const std::vector<int>& type_buffer, const std::vector<int>& category_buffer){

   double io_time;

   switch (config::internal::format){

      case config::internal::binary:
         io_time = write_coord_data_binary(filename, buffer, type_buffer, category_buffer);
         break;

      case config::internal::text:
         io_time = write_coord_data_text(filename, buffer, type_buffer, category_buffer);
         break;

   }

   return io_time;

}

//----------------------------------------------------------------------------------------------------
// Function to output coord data formatted as text
//----------------------------------------------------------------------------------------------------
//
double write_coord_data_text(std::string filename, const std::vector<double> &buffer, const std::vector<int> &type_buffer, const std::vector<int> &category_buffer){

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

      ofile << type_buffer[index] << "\t" <<
               category_buffer[index] << "\t" <<
               buffer[3 * index + 0] << "\t" <<
               buffer[3 * index + 1] << "\t" <<
               buffer[3 * index + 2] << "\n";

   }

   // end timer
   timer.stop();

   // close output file
   ofile.close();

   // return elapsed time for io
   return timer.elapsed_time();

}

//----------------------------------------------------------------------------------------------------
// Function to output coord data in binary format
//----------------------------------------------------------------------------------------------------
//
double write_coord_data_binary(std::string filename, const std::vector<double> &buffer, const std::vector<int> &type_buffer, const std::vector<int> &category_buffer){
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

   // output buffers to disk
   ofile.write(reinterpret_cast<const char *>(&type_buffer[0]), sizeof(int) * type_buffer.size());
   ofile.write(reinterpret_cast<const char *>(&category_buffer[0]), sizeof(int) * category_buffer.size());
   ofile.write(reinterpret_cast<const char *>(&buffer[0]), sizeof(double) * buffer.size());

   // end timer
   timer.stop();

   // close output file
   ofile.close();

   // return elapsed time for io
   return timer.elapsed_time();

}

} // end of namespace internal
} // end of namespace config
