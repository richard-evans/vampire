//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
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
   bool verbose = true; // flag to specify verbosity of output to user
   bool xyz = true; // flag to specify crystal.xyz file output
   bool povray = true; // flag to specify povray file output
   bool cells = false; // flag to specify cells output
   bool vtk = true;
   bool x_vector = false; // flag to specify direction of povray colouring
   bool z_vector = false; // flag to specify plane for povray colouring

   // keyword variables
   std::string colour_keyword = "C2";
   std::string custom_colourmap_file;
   bool z_axis_colour = true;

   format_t format;

   uint64_t num_atoms = 0;

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

} // end of namespace vdc
