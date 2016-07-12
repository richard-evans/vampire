//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iomanip>
#include <fstream>
#include <sstream>

// Vampire headers
#include "config.hpp"
#include "vio.hpp"
#include "vutil.hpp"

// config headers
#include "internal.hpp"

//----------------------------------------------
// Packed struct data type for binary buffer
//----------------------------------------------
struct coord_t{
   // standard types to ensure cross-compiler portability
   uint32_t material;
   uint32_t category;
   // position vectors
   float x;
   float y;
   float z;
};

namespace config{
   namespace internal{

// Rory these functions need fixing

      // forward function declarations
      void write_coordinate_data_text(const std::vector<double>& spins_cx, // spin coordinates (Angstroms)
                                      const std::vector<double>& spins_cy,
                                      const std::vector<double>& spins_cz,
                                      const std::vector<int>& material, // material id
                                      const std::vector<int>& category, // category id
                                      const std::string filename);

      void write_coordinate_data_binary(const std::vector<double>& spins_cx, // spin coordinates (Angstroms)
                                        const std::vector<double>& spins_cy,
                                        const std::vector<double>& spins_cz,
                                        const std::vector<int>& material, // material id
                                        const std::vector<int>& category, // category id
                                        const std::string filename);

      //----------------------------------------------------------------------------------------------------
      // Simple wrapper function to call output function for correct format
      //----------------------------------------------------------------------------------------------------
      //
      void write_coordinate_data(const std::vector<double>& spins_cx, // spin coordinates (Angstroms)
                                 const std::vector<double>& spins_cy,
                                 const std::vector<double>& spins_cz,
                                 const std::vector<int>& material, // material id
                                 const std::vector<int>& category){ // category id

         // determine filename
         std::stringstream filename;
         filename << "atoms-coords.cfg";

         // Output informative message to log file
         zlog << zTs() << "Outputting configuration file " << filename.str() << " to disk ";

         // copy total number of output data to const for compiler
         const unsigned int num_data = config::internal::local_output_atom_list.size();

         // initialize temporary buffers
         std::vector<double> x_buffer(num_data);
         std::vector<double> y_buffer(num_data);
         std::vector<double> z_buffer(num_data);
         std::vector<double> mat_buffer(num_data);
         std::vector<double> cat_buffer(num_data);

         // loop over all atoms to be output
         for(int id=0; id < num_data; ++id){

            // determine next datum to be output
            const int index = config::internal::local_output_atom_list[id];

            // copy and cast data to be output to main output buffer
            //buffer[id].material = material[index];
            //buffer[id].category = category[index];
            //buffer[id].x = float(spins_cx[index]);
            //buffer[id].y = float(spins_cy[index]);
            //buffer[id].z = float(spins_cz[index]);

         }

         // pack cordinate_buffer(); + MPI

         switch(config::internal::output_data_format){

            case config::internal::binary:
               //write_coordinate_data_binary(buffer, filename.str());
               //write_coordinate_data_binary(buffer, filename.str());
               break;

            case config::internal::text:
               //write_coordinate_data_text(buffer, filename.str());
               break;

         }

         return;

      }

      //----------------------------------------------------------------------------------------------------
      // Function to output coordinate data formatted as text
      //----------------------------------------------------------------------------------------------------
      //
      void write_coordinate_data_text(const std::vector<double>& spins_cx, // spin coordinates (Angstroms)
                                      const std::vector<double>& spins_cy,
                                      const std::vector<double>& spins_cz,
                                      const std::vector<int>& material, // material id
                                      const std::vector<int>& category, // category id
                                      const std::string filename){ // file name

         #ifdef MPICF


            // Set CPUID on non-root process
         //   filename << std::setfill('0') << std::setw(5) << config vmpi::my_rank << "-";

         #else

            // mpi::syncronize_coords(cx, cy, cz, material, category);

            // instantiate timer
            vutil::vtimer_t timer;

            // start timer
            timer.start();

            // Declare and open output file
            std::ofstream ofile;
            ofile.open (filename.c_str());

            // copy total number of output data to const for compiler
            const int num_data = config::internal::local_output_atom_list.size();

            // output number of data
            ofile << num_data << "\n";

            // loop over all atoms to be output
            for(int id=0; id < num_data; ++id){

               // determine next datum to be output
               const int index = config::internal::local_output_atom_list[id];

               // output data
               ofile << material[index] << "\t" << category[index] << "\t" << spins_cx[index] << "\t" << spins_cy[index] << "\t" << spins_cz[index] << "\n";

            }

            // close output file
            ofile.close();

            // stop the timer
            double total_time = timer.elapsed_time(); // seconds

            // open file at end
            std::ifstream in(filename.c_str(), std::ios::binary | std::ios::ate);

            // get file size (bytes)
            double data_size = double(in.tellg());

            // close file
            in.close();

            // calculate data rate and output to log
            zlog << 1.0e-6*data_size/total_time << " MB/s" << std::endl;

         #endif

         return;

      }

      //----------------------------------------------------------------------------------------------------
      // Function to output coordinate data in binary format
      //----------------------------------------------------------------------------------------------------
      //
      void write_coordinate_data_binary(const std::vector<double>& spins_cx, // spin coordinates (Angstroms)
                                        const std::vector<double>& spins_cy,
                                        const std::vector<double>& spins_cz,
                                        const std::vector<int>& material, // material id
                                        const std::vector<int>& category, // category id
                                        const std::string filename){ // file name

         #ifdef MPICF


            // Set CPUID on non-root process
         //   filename << std::setfill('0') << std::setw(5) << config vmpi::my_rank << "-";

         #else

            //--------------------------
            // fill temporary buffer
            //--------------------------





            // instantiate timer
            vutil::vtimer_t timer;

            // start timer
            timer.start();

            // Declare and open output file
            std::ofstream ofile;
            ofile.open (filename.c_str(),std::ios::binary);

            // determine number of data to output
            //const uint64_t buffer_size = buffer.size();

            // output number of data
            //ofile.write(reinterpret_cast<const char*>(&buffer_size),sizeof(uint64_t));

            // output buffer to disk
            //ofile.write(reinterpret_cast<const char*>(&buffer[0]),sizeof(coord_t)*buffer.size());

            // close output file
            ofile.close();

            // stop the timer
            //double total_time = timer.elapsed_time(); // seconds

            // get file size (bytes)
            //double data_size = double(sizeof(coord_t)*buffer.size());

            // calculate data rate and output to log
            //zlog << 1.0e-6*data_size/total_time << " MB/s" << std::endl;

         #endif

         return;

      }

   } // end of namespace internal
} // end of namespace config
