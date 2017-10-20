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

   extern unsigned int start_file_id;
   extern unsigned int final_file_id;

   extern double system_size[3];
   extern double system_centre[3];

   extern std::vector<material_t> materials;

   extern std::vector<int> category;
   extern std::vector<int> type;

   extern std::vector<double> coordinates;
   extern std::vector<double> spins;

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
   void process_coordinates();
   void process_spins();

   // forward function declarations
   void read_nm_metadata();
   void read_nm_data();

   void output_xyz_file();
   void output_inc_file(unsigned int spin_file_id);
   void output_povray_file();

   void initialise_cells();
   void output_cell_file(unsigned int spin_file_id);

   void colourwheel( std::vector<std::vector<double>>& colourmap );
   void rgb( const double& sx, const double& sy, const double& sz, double &red, double &green, double &blue);
   double scale( double scale_min, double scale_max, double start_min, double start_max, double x );
   void rgb2hsl( double& red, double& green, double& blue, double& hue, double& light, double& saturation );
   void hsl2rgb( double& red, double& green, double& blue, double& hue, double& light, double& saturation );
   void interpolate(double xy_angle, double& red, double& green, double& blue );
   void rgb2hsi( double& red, double& green, double& blue, double& h2, double& intensity, double& saturation );
   void hsi2rgb( double& red, double& green, double& blue, double& hue, double& intensity, double& saturation );
}

#endif //VDC_H_
