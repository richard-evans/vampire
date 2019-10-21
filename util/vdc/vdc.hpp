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

#ifndef VDC_H_
#define VDC_H_

// C++ standard library headers
#include <iostream>
#include <vector>

namespace vdc{

   // program option flags
   extern bool verbose;
   extern bool xyz;
   extern bool povray;
   extern bool cells;
   extern bool vtk;
   extern bool txt;
   extern bool x_vector;
   extern bool z_vector;

   // keyword variables
   extern std::string colour_keyword;
   extern std::string custom_colourmap_file;
   extern bool x_axis_colour;

   // enumerated integers for option selection
   enum format_t{ binary = 0, text = 1};
   extern format_t format;

   // simple struct to store material parameters
   struct material_t{
      int id;
      double moment;
      std::string name;
      std::string element;
   };

   extern uint64_t num_atoms;

   extern unsigned int vdc_start_file_id;
   extern unsigned int vdc_final_file_id;
   extern unsigned int start_file_id;
   extern unsigned int final_file_id;

   extern double system_size[3];
   extern double system_centre[3];

   extern std::vector<material_t> materials;

   extern std::vector<int> category;
   extern std::vector<int> type;

   extern std::vector<double> coordinates;
   extern std::vector<double> spins;

   // axis vectors for povray colouring
   extern std::vector<double> vector_z;
   extern std::vector<double> vector_y;
   extern std::vector<double> vector_x;

   // non-magnetic atom data
   extern uint64_t num_nm_atoms;
   extern std::vector<int> nm_category;
   extern std::vector<int> nm_type;
   extern std::vector<double> nm_coordinates;

   extern unsigned int total_cells;
   extern unsigned int nx_cells;
   extern unsigned int ny_cells;
   extern unsigned int nz_cells;

   extern std::vector<int> atom_cell_id;
   extern std::vector<double> cell_coords;
   extern std::vector< std::vector< std::vector <double> > > cell_magnetization;

   // array to store subsidiary data file names
   extern std::vector <std::string> coord_filenames;
   extern std::vector <std::string> spin_filenames;
   extern std::vector <std::string> nm_filenames;

   // Functions
   int command( int argc, char* argv[]);
   void process_coordinates();
   void process_spins();

   // forward function declarations
   void read_nm_metadata();
   void read_nm_data();

   void output_xyz_file();
   void output_inc_file(unsigned int spin_file_id);
   void output_povray_file();
   void output_vtk_file(unsigned int spin_file_id);
   void output_txt_file(unsigned int file_id);

   void initialise_cells();
   void output_cell_file(unsigned int spin_file_id);

   void rgb( const double& sx, const double& sy, const double& sz, double &red, double &green, double &blue);

   int colourwheel ( std::vector<std::vector<double>>& colourmap );

}

#endif //VDC_H_
