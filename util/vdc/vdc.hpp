//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans, Daniel Meilak 2017-2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef VDC_H_
#define VDC_H_

// C++ standard library headers
#include <iostream>
#include <vector>
#include <unordered_map>
#include <functional>
#include <string>

namespace vdc{

   // simple struct for xy coordinates
   struct xy_t{
      double x;
      double y;
   };

   // input filename
   extern std::string input_file;

   // program option flags
   extern bool verbose;
   extern bool xyz;
   extern bool grains; // flag to enable grain calculations
   extern bool povray;
   extern bool povcells;
   extern bool povsticks;
   extern bool cells;
   extern bool cellsf;
   extern bool vtk;
   extern bool ssc; // flag to specify spin-spin correlation
   extern bool txt;
   extern bool x_vector;
   extern bool z_vector;

   // keyword variables
   extern std::string colour_keyword;
   extern std::string custom_colourmap_file;
   extern std::vector<std::vector<double>> colourmap;
   extern std::vector<std::string> colourmaps;
   extern bool x_axis_colour;
   extern bool default_camera_pos;

   // enumerated integers for option selection
   enum format_t{ binary = 0, text = 1};
   enum slice_type{ box, box_void, sphere, cylinder};
   extern format_t format;

   // list of input file parameters set in command line (to check for double usage)
   extern std::vector<std::string> cmdl_parameters;

   // simple struct to store material parameters
   struct material_t{
      int id = 0;
      double moment = 1.0;
      std::string name = "material";
      std::string element = "H";
   };

   // struct to hold input parameters, value and line in input file for error check
   struct input_t{
      std::string key;
      std::vector<std::string> value;
      int line_number;
   };

   // struct to hold slice information
   struct slice_t{
      enum slice_type type;
      std::vector<double> param;
      std::vector<double> bound;
   };

   // vector of slices
   extern std::vector<slice_t> slices;

   // unordered map of input keys to function wrappers
   extern const std::unordered_map<std::string,std::function<void(const input_t&)>> key_list;

   extern uint64_t num_atoms;

   // two sets of ids required
   extern unsigned int vdc_start_file_id;
   extern unsigned int vdc_final_file_id;
   extern unsigned int start_file_id;
   extern unsigned int final_file_id;

   extern double system_size[3];
   extern double system_centre[3];

   // slice parameters for cutting the original system
   extern std::vector<double> slice_parameters;
   extern std::vector<int> remove_materials;
   extern std::vector<int> afm_materials;
   extern std::vector<int> atoms_list;
   extern std::vector<int> nm_atoms_list;
   extern std::vector<int> sliced_atoms_list;
   extern std::vector<int> sliced_nm_atoms_list;

   extern std::vector<material_t> materials;

   extern std::vector<int> category;
   extern std::vector<int> type;
   extern std::vector<int> grain;

   extern std::vector<double> coordinates;
   extern std::vector<double> spins;

   // axis vectors for povray colouring
   extern std::vector<double> vector_z;
   extern std::vector<double> vector_y;
   extern std::vector<double> vector_x;

   // povray camera settings
   extern std::vector<double> camera_pos;
   extern std::vector<double> camera_look_at;
   extern double camera_zoom;
   extern std::string background_colour;

   // povray sticks settings
   extern double sticks_cutoff;

   // povray shape sizes
   extern std::vector<double> atom_sizes;
   extern std::vector<double> arrow_sizes;

   // non-magnetic atom data
   extern uint64_t num_nm_atoms;
   extern std::vector<int> nm_category;
   extern std::vector<int> nm_type;
   extern std::vector<int> nm_grain;
   extern std::vector<double> nm_coordinates;

   // cell data
   extern double cell_size[3]; // Angstroms
   extern unsigned int total_cells;
   extern unsigned int nx_cells;
   extern unsigned int ny_cells;
   extern unsigned int nz_cells;

   extern std::vector<int> atom_cell_id;
   extern std::vector<int> num_atoms_in_cell;
   extern std::vector<double> cell_coords;
   extern std::vector< std::vector< std::vector <double> > > cell_magnetization;

   // grain data
   extern std::vector < std::vector <xy_t> > grain_vertices_array;

   // array to store subsidiary data file names
   extern std::vector <std::string> coord_filenames;
   extern std::vector <std::string> spin_filenames;
   extern std::vector <std::string> nm_filenames;

   // arrays for storing time-averaged spin-spin correlations
   extern std::vector<double> ssc_counts; // number of counts
   extern std::vector<double> ssc_correl; // sum of correlations
   extern double ssc_magnetization; // sum snapshot magnetizations
   extern double ssc_snapshots; // number of snapshots
   extern double ssc_num_bins;  // number of bins for correlations
   extern double ssc_bin_width; // width of each bin (Agstroms)
   extern double ssc_inv_bin_width; // 1/bin width

   //==========================================
   // Forward function declarations
   //==========================================

   // main
   void command( int argc, char* argv[]);
   void read_and_set();
   void process_coordinates();
   void process_spins();

   // non-magnetic
   void read_nm_metadata();
   void read_nm_data();
   void slice_nm_system();

   // XYZ
   void output_xyz_file();

   // atoms
   void output_atoms_txt_file();

   // VTK
   void output_vtk_file(unsigned int spin_file_id);

   // TXT
   void output_txt_file(unsigned int spin_file_id);

   // Povray
   void initialise_povray();
   void output_inc_file(unsigned int spin_file_id);
   void output_povray_file();
   void output_cells_inc_file(unsigned int spin_file_id);
   void output_povray_cells_file();
   void output_sticks_file();

   // grains
   void load_grain_vertices();
   void determine_atom_grain_id();
   //void generate_povray_grains()

   // Colour
   void rgb( const double& sx, const double& sy, const double& sz, double &red, double &green, double &blue);
   void initialise_colourwheel();

   // SSC
   void initialise_ssc();
   void output_average_ssc_file();
   void output_ssc_file(unsigned int spin_file_id);

   // CELL
   void initialise_cells();
   void output_cell_file(unsigned int spin_file_id);

   // setting functions
   void set_frame_start(const input_t &input);
   void set_frame_final(const input_t &input);
   void set_remove_materials(const input_t &input);
   void set_afm(const input_t &input);
   void set_slice(const input_t &input);
   void set_slice_void(const input_t &input);
   void set_slice_sphere(const input_t &input);
   void set_slice_cylinder(const input_t &input);
   void set_vector_z(const input_t &input);
   void set_vector_x(const input_t &input);
   void set_cell_size(const input_t &input);
   void set_colourmap(const input_t &input);
   void set_custom_colourmap(const input_t &input);
   void set_3D(const input_t &input);
   void set_camera_position(const input_t &input);
   void set_camera_look_at(const input_t &input);
   void set_camera_zoom(const input_t &input);
   void set_sticks_cutoff(const input_t &input);
   void set_background_colour(const input_t &input);
   void set_atom_sizes(const input_t &input);
   void set_arrow_sizes(const input_t &input);
}

#endif //VDC_H_
