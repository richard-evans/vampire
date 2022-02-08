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

   // program option flags
   bool verbose  = false; // flag to specify verbosity of output to user
   bool xyz      = false; // flag to specify crystal.xyz file output
   bool povray   = false; // flag to specify povray file output
   bool cells    = false; // flag to specify cells output
   bool vtk      = false; // flag to specify vtk output
   bool ssc      = false; // flag to specify spin-spin correlation
   bool txt      = false; // flag to specify plain text output
   bool x_vector = false; // flag to specify direction of povray colouring
   bool z_vector = false; // flag to specify plane for povray colouring

   // keyword variables
   std::string colour_keyword = "CBWR";
   std::string custom_colourmap_file;
   bool x_axis_colour = false;
   std::string slice_type = "no-slice";

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

   // axis vectors for povray colouring
   std::vector<double> vector_z = {0.0,0.0,1.0};
   std::vector<double> vector_y = {0.0,1.0,0.0};
   std::vector<double> vector_x = {1.0,0.0,0.0};

   // non-magnetic atom data
   uint64_t num_nm_atoms = 0;
   std::vector<int> nm_category(0);
   std::vector<int> nm_type(0);
   std::vector<double> nm_coordinates(0);

   // cell data
   unsigned int total_cells = 0;
   unsigned int nx_cells = 1;
   unsigned int ny_cells = 1;
   unsigned int nz_cells = 1;

   std::vector<int> atom_cell_id;
   std::vector<double> cell_coords;
   std::vector< std::vector< std::vector <double> > > cell_magnetization;

   // array to store subsidiary data file names
   std::vector <std::string> coord_filenames(0);
   std::vector <std::string> spin_filenames(0);
   std::vector <std::string> nm_filenames(0);

   // arrays for storing time-averaged spin-spin correlations
   std::vector<double> ssc_counts(0); // number of counts
   std::vector<double> ssc_correl(0); // sum of correlations
   double ssc_magnetization; // sum snapshot magnetizations
   double ssc_snapshots; // number of snapshots
   double ssc_num_bins;  // number of bins for correlations
   double ssc_bin_width; // width of each bin (Agstroms)
   double ssc_inv_bin_width; // 1/bin width

} // end of namespace vdc
