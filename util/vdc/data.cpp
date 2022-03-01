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

// C++ standard library headers

// program header
#include "vdc.hpp"

namespace vdc{

   // input filename
   std::string input_file = "vdc_input";

   // program option flags
   bool verbose  = false; // flag to specify verbosity of output to user
   bool xyz      = false; // flag to specify crystal.xyz file output
   bool povray   = false; // flag to specify povray file output
   bool povcells = false; // flag to specify povray cells output
   bool cells    = false; // flag to specify cells output
   bool cellsf   = false; // flag to output cell file
   bool vtk      = false; // flag to specify vtk output
   bool ssc      = false; // flag to specify spin-spin correlation
   bool txt      = false; // flag to specify plain text output
   bool x_vector = false; // flag to specify direction of povray colouring
   bool z_vector = false; // flag to specify plane for povray colouring

   // keyword variables
   std::string colour_keyword = "cbwr";
   std::string custom_colourmap_file;
   std::vector<std::vector<double>> colourmap(256, std::vector<double>(3));
   std::vector<std::string> colourmaps = {"c2", "bwr", "cbwr", "rainbow"};
   bool x_axis_colour = false;
   bool default_camera_pos = true;

   // list of input file parameters set in command line (to check for double usage)
   std::vector<std::string> cmdl_parameters;

   format_t format;

   uint64_t num_atoms = 0;

   unsigned int vdc_start_file_id = 0;
   unsigned int vdc_final_file_id = 99999999;
   unsigned int start_file_id = 0;
   unsigned int final_file_id = 99999999;

   // system size and centre
   double system_size[3] = {0.0, 0.0, 0.0};
   double system_centre[3] = {0.0, 0.0, 0.0};

   std::vector<material_t> materials(0);

   std::vector<int> category(0);
   std::vector<int> type(0);

   std::vector<double> coordinates(0);
   std::vector<double> spins(0);

   // slice parameters for cutting the original system
   std::vector<double> slice_parameters = {0.0,1.0,0.0,1.0,0.0,1.0};
   std::vector<int> remove_materials(0);
   std::vector<int> afm_materials(0);
   std::vector<int> atoms_list(0);
   std::vector<int> nm_atoms_list(0);
   std::vector<int> sliced_atoms_list(0);
   std::vector<int> sliced_nm_atoms_list(0);

   // user defined slices
   std::vector<slice_t> slices;

   // axis vectors for povray colouring
   std::vector<double> vector_z = {0.0,0.0,1.0};
   std::vector<double> vector_y = {0.0,1.0,0.0};
   std::vector<double> vector_x = {1.0,0.0,0.0};

   // povray camera settings
   std::vector<double> camera_pos    = {0.0,0.0,1.0};
   std::vector<double> camera_look_at = {0.0,0.0,0.0};
   double camera_zoom = 1.0;
   std::string background_colour = "Gray30";

   // povray shape sizes
   std::vector<double> atom_sizes  = {1.2};
   std::vector<double> arrow_sizes = {2.0};

   // non-magnetic atom data
   uint64_t num_nm_atoms = 0;
   std::vector<int> nm_category(0);
   std::vector<int> nm_type(0);
   std::vector<double> nm_coordinates(0);

   // cell data
   double cell_size[3] = { 10.0, 10.0, 10.0 }; // Angstroms
   unsigned int total_cells = 0;
   unsigned int nx_cells = 1;
   unsigned int ny_cells = 1;
   unsigned int nz_cells = 1;

   std::vector<int> atom_cell_id;
   std::vector<int> num_atoms_in_cell;
   std::vector<double> cell_coords;
   std::vector<std::vector<std::vector <double>>> cell_magnetization;

   // array to store subsidiary data file names
   std::vector<std::string> coord_filenames(0);
   std::vector<std::string> spin_filenames(0);
   std::vector<std::string> nm_filenames(0);

   // arrays for storing time-averaged spin-spin correlations
   std::vector<double> ssc_counts(0); // number of counts
   std::vector<double> ssc_correl(0); // sum of correlations
   double ssc_magnetization; // sum snapshot magnetizations
   double ssc_snapshots; // number of snapshots
   double ssc_num_bins;  // number of bins for correlations
   double ssc_bin_width; // width of each bin (Agstroms)
   double ssc_inv_bin_width; // 1/bin width

   // unordered map of input key string to function wrapper
   const std::unordered_map<std::string, std::function<void(const input_t&)>> key_list = {

      // frame start and end
      {"frame-start", set_frame_start},
      {"start-frame", set_frame_start},
      {"frame-final", set_frame_final},
      {"final-frame", set_frame_final},
      {"frame-end"  , set_frame_final},
      // remove materials
      {"remove-material" , set_remove_materials},
      {"remove-materials", set_remove_materials},
      // antiferromagnetic materials
      {"antiferromagnetic-materials", set_afm},
      {"afm", set_afm},
      // slices
      {"slice"         , set_slice},
      {"slice-void"    , set_slice_void},
      {"void-slice"    , set_slice_void},
      {"slice-sphere"  , set_slice_sphere},
      {"sphere-slice"  , set_slice_sphere},
      {"slice-cylinder", set_slice_cylinder},
      // cells
      {"cell-size", set_cell_size},
      // colourmap
      {"vector-z" , set_vector_z},
      {"z-vector" , set_vector_z},
      {"vector-x" , set_vector_x},
      {"x-vector" , set_vector_x},
      {"colourmap", set_colourmap},
      {"colormap" , set_colourmap},
      {"custom-colormap" , set_custom_colourmap},
      {"custom-colourmap", set_custom_colourmap},
      {"3d", set_3D},
      // povray camera settings
      {"camera-position", set_camera_position},
      {"camera-look-at" , set_camera_look_at},
      {"camera-lookat"  , set_camera_look_at},
      {"camera-zoom"    , set_camera_zoom},
      // povray background colour
      {"background-colour", set_background_colour},
      {"background-color" , set_background_colour},
      // povray shape sizes
      {"atom-sizes" , set_atom_sizes},
      {"atom-size"  , set_atom_sizes},
      {"arrow-sizes", set_arrow_sizes},
      {"arrow-size" , set_arrow_sizes}
   };

} // end of namespace vdc
