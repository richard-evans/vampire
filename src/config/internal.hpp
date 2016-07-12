#ifndef CONFIG_INTERNAL_H_
#define CONFIG_INTERNAL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// configuration output implementation. These functions should 
// not be accessed outside of the configuration output code.
//---------------------------------------------------------------------
namespace config{
   namespace internal{
      //-----------------------------------------------------------------------------
      // Shared data structures used for configuration output
      //-----------------------------------------------------------------------------
      enum data_format_t {text, binary};
      
      //-----------------------------------------------------------------------------
      // Shared variables used for configuration output
      //-----------------------------------------------------------------------------
      extern bool output_atoms; // flag to control atomic output
      extern bool output_cells; // flag to control macrocell output
      extern bool output_meta; // flag to enable output of coordinate/meta data
      extern bool output_coords; // flag to control atomic coordinate output
      extern int output_rate; // relative rate of file output
      extern int output_file_counter; // configuration file number
      extern int output_rate_counter; // configuration output counter

      extern data_format_t output_data_format;
      
      extern double atoms_output_min[3]; // Spatial range fr atomic output
      extern double atoms_output_max[3];
      extern std::vector<int> local_output_atom_list; // list of atoms to be output according to spatial range      
      extern int total_output_atoms; // Total number of atoms to be output on local node

      // Buffer variables to store copies of data in float format for reduced file size
      extern std::vector<float> output_spin_buffer; // buffer to store float cast spin array for output to disk

      extern int num_io_nodes; // total number of i/o processes
      extern bool io_master; // flag set to true if I output data
      extern bool set_num_io_nodes_to_ppn; // flag to post initialise num_io_nodes

      //-----------------------------------------------------------------------------
      // Shared functions used for configuration output
      //-----------------------------------------------------------------------------

      // Function to write meta data for each configuration
      void write_meta(const double simulation_time, // time (seconds)
                      const double temperature, // system temperature (Kelvin)
                      const double applied_field_x, // applied field components (Tesla)
                      const double applied_field_y,
                      const double applied_field_z,
                      const double magnetization_x, // magnetization components (normalized)
                      const double magnetization_y,
                      const double magnetization_z);
      
      void copy_data_to_buffer(const std::vector<double>& x, // vector data
                               const std::vector<double>& y,
                               const std::vector<double>& z,
                               const std::vector<int>& mask,
                               std::vector<float>& buffer);

      void write_data(std::string filename, const std::vector<float>& buffer);

      void write_coordinate_meta();

      void write_coordinate_data(const std::vector<double>& spins_cx, // spin coordinates (Angstroms)
                                 const std::vector<double>& spins_cy,
                                 const std::vector<double>& spins_cz,
                                 const std::vector<int>& material, // material id
                                 const std::vector<int>& category);

      namespace mpi{
         void initialize();
      }

   } // end of internal namespace
} // end of config namespace

#endif //CONFIG_INTERNAL_H_
