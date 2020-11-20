#ifndef LTMP_INTERNAL_H_
#define LTMP_INTERNAL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// local temperature pulse implementation. These functions should
// not be accessed outside of the local temperature pulse code.
//---------------------------------------------------------------------
namespace ltmp{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Shared variables used for the localised temperature pulse calculation
      //-----------------------------------------------------------------------------
      extern bool enabled; // enable localised temperature pulse calculation
      extern bool initialised; /// flag set if initialised
      extern bool lateral_discretisation; /// enable lateral temperature profile
      extern bool vertical_discretisation; /// enable vertical temperature profile
      extern bool output_microcell_data; /// enable verbose output data for temperature cells
      extern bool temperature_rescaling; /// enable rescaled temperature calculation
      extern bool gradient; /// enable temperature gradient

      extern double micro_cell_size; /// lateral size of local temperature microcells (A)
      extern double laser_spot_size; /// laser spot size for lateral profile (A)
      extern double penetration_depth; /// vertical laser penetration depth
      extern double thermal_conductivity; //J/s/m/K

      extern double pump_power; // laser pump power
      extern double pump_time; // laser pump time (s)
      extern double TTG;  // electron-lattice coupling constant
      extern double TTCe; // electron heat capacity (T=0)
      extern double TTCl; // lattice heat capcity
      extern double dt; // time step

      extern double minimum_temperature; // Minimum temperature in temperature gradient
      extern double maximum_temperature; // Maximum temperature in temperature gradient

      extern int num_local_atoms; /// number of local atoms (ignores halo atoms in parallel simulation)
      extern int num_cells; /// number of temperature cells
      extern int my_first_cell; /// first cell on my CPU
      extern int my_last_cell; /// last cell on my CPU

      extern std::vector<int> atom_temperature_index; /// defines which temperature cell applies to atom (including Te or Tp)
      extern std::vector<double> atom_sigma; /// unrolled list of thermal prefactor sqrt(2kBalpha/gamma*mu_s*dt)
      extern std::vector<double> atom_rescaling_root_Tc; /// unrolled list of material Curie temperature for rescaling calculation
      extern std::vector<double> atom_rescaling_alpha; /// unrolled list of material rescaling exponent

      extern std::vector<int> cell_neighbour_list; // list of cell interactions for heat transfer
      extern std::vector<int> cell_neighbour_start_index; // start index of interactions for cell
      extern std::vector<int> cell_neighbour_end_index; // end index of interactions for cell

      extern std::vector<double> x_field_array; /// arrays to store atomic temperature field
      extern std::vector<double> y_field_array;
      extern std::vector<double> z_field_array;

      extern std::vector<double> root_temperature_array; /// stored as pairs sqrt(Te), sqrt(Tp) (2 x number of cells) MIRRORED on all CPUs
      extern std::vector<double> cell_position_array; /// position of cells in x,y,z (3*n) MIRRORED on all CPUs // dont need this
      extern std::vector<double> delta_temperature_array; /// stored as pairs dTe, dTp LOCAL CPU only
      extern std::vector<double> attenuation_array; /// factor reducng incident laser fluence for each cell LOCAL CPU only

      extern std::vector<double> material_kerr_sensitivity_depth; // unrolled list of kerr sensitivity depths for each material

      void write_microcell_data();
      void open_vertical_temperature_profile_file();
      void open_lateral_temperature_profile_file();
      void write_cell_temperature_data();
      void calculate_local_temperature_pulse(const double time_from_start);
      void calculate_local_temperature_gradient();

   } // end of iternal namespace
} // end of st namespace

#endif //LTMP_INTERNAL_H_
