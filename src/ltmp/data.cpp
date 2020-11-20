//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "ltmp.hpp"

// Local temperature pulse headers
#include "internal.hpp"

namespace ltmp{

   //-----------------------------------------------------------------------------------------------
   // Externally visible variables
   //-----------------------------------------------------------------------------------------------
   abs_t absorption_profile; // class variable containing tabulated absorption profile

   namespace internal{

      //-----------------------------------------------------------------------------
      // Shared variables used for the local temperature pulse calculation
      //-----------------------------------------------------------------------------
      bool enabled=false; /// enable localised local temperature pulse
      bool initialised=false; /// flag set if initialised
      bool lateral_discretisation=false; /// enable lateral temperature profile
      bool vertical_discretisation=true; /// enable vertical temperature profile
      bool output_microcell_data=false; /// enable verbose output data for temperature cells
      bool temperature_rescaling=false; /// enable rescaled temperature calculation
      bool gradient=false; /// enable temperature gradient

      double micro_cell_size = 10.0; /// lateral size of local temperature microcells (A)
      double laser_spot_size = 350.0; /// laser spot size for lateral profile (A)
      double penetration_depth = 200.0; /// vertical laser penetration depth
      double thermal_conductivity = 11.0; //J/s/m/K

      double pump_power; // laser pump power
      double pump_time; // laser pump time (s)
      double TTG;  // electron-lattice coupling constant
      double TTCe; // electron heat capacity (T=0)
      double TTCl; // lattice heat capcity
      double dt; // time step

      double minimum_temperature = 0.0; // Minimum temperature in temperature gradient
      double maximum_temperature = 0.0; // Maximum temperature in temperature gradient

      int num_local_atoms; /// number of local atoms (ignores halo atoms in parallel simulation)
      int num_cells; /// number of temperature cells
      int my_first_cell; /// first cell on my CPU
      int my_last_cell; /// last cell on my CPU

      std::vector<int> atom_temperature_index; /// defines which temperature cell applies to atom (including Te or Tp)
      std::vector<double> atom_sigma; /// unrolled list of thermal prefactor sqrt(2kBalpha/gamma*mu_s*dt)
      std::vector<double> atom_rescaling_root_Tc; /// unrolled list of material Curie temperature for rescaling calculation
      std::vector<double> atom_rescaling_alpha; /// unrolled list of material rescaling exponent

      std::vector<int> cell_neighbour_list; // list of cell interactions for heat transfer
      std::vector<int> cell_neighbour_start_index; // start index of interactions for cell
      std::vector<int> cell_neighbour_end_index; // end index of interactions for cell

      std::vector<double> x_field_array; /// arrays to store atomic temperature field
      std::vector<double> y_field_array;
      std::vector<double> z_field_array;

      std::vector<double> root_temperature_array; /// stored as pairs sqrt(Te), sqrt(Tp) (2 x number of cells) MIRRORED on all CPUs
      std::vector<double> cell_position_array; /// position of cells in x,y,z (3*n) MIRRORED on all CPUs // dont need this
      std::vector<double> delta_temperature_array; /// stored as pairs dTe, dTp LOCAL CPU only
      std::vector<double> attenuation_array; /// factor reducing incident laser fluence for each cell LOCAL CPU only

      //std::vector<double> material_kerr_sensitivity_depth; // unrolled list of kerr sensitivity depths for each material

   } // end of internal namespace
} // end of ltmp namespace

