//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and Rory Pond 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk and rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef VIO_INTERNAL_H_
#define VIO_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the vio module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <string>
#include <vector>
#include <iostream>

// Vampire headers

#ifdef WIN_COMPILE
   #include <direct.h>
#endif

namespace vin{

   //----------------------------------------------------------------------------------
   //Funciton protypes for functions inside: match.cpp
   //----------------------------------------------------------------------------------
   int match(std::string const, std::string const, std::string const, std::string const, int const);
   int match(std::string const, std::string const, std::string const, std::string const, int const, int const, int const);
   int match_create(std::string const, std::string const, std::string const, int const);
   int match_dimension(std::string const, std::string const, std::string const, int const);
   int match_sim(std::string const, std::string const, std::string const, int const);
   int match_vout_list(std::string const, std::string const, int const, std::vector<unsigned int> &);
   int match_material(std::string const, std::string const, std::string const, int const, int const, int const, std::string const, std::string const);
   int match_config(std::string const, std::string const, std::string const, int const);

   //----------------------------------------------------------------------------------
   //Function protypes for functions inside: read.cpp
   //----------------------------------------------------------------------------------
   int read(std::string const filename);
   int read_mat_file(std::string const, int const);

}

namespace vout{

   extern std::string zLogProgramName; /// Program Name
   extern std::string zLogHostName; /// Host Name
   extern bool        zLogInitialised; /// Initialised flag
   #ifdef WIN_COMPILE
   	extern int      zLogPid; /// Process ID
   #else
   	extern pid_t    zLogPid; /// Process ID
   #endif

   // namespaced io lists (to avoid collisions)
   namespace grain{
      // defined enumerated types
      enum output_t {
         time_steps,
         real_time,
         temperature,
         electron_temperature,
         phonon_temperature,
         applied_field,
         applied_field_unit_vector,
         constraint_phi,
         constraint_theta,
         magnetisation,
         material_magnetisation,
         material_height_magnetisation,
         mean_magnetisation_length,
         mean_specific_heat,
         mean_susceptibility,
         mean_torque
      };

      // internal variables
      extern int output_rate; // rate of output compared to calculation

      // grain output list
      extern std::vector<grain::output_t> output_list;

   }

   //-------------------------------------------------------------------------
   // New output functions
   //-------------------------------------------------------------------------
   void write_grain_file();

   // formatting wrapper functions
   std::string generic_output_double(const std::string str, const double d, const bool header);

   //-------------------------------------------------------------------------
   // New match functions
   //-------------------------------------------------------------------------
   int match_vout_grain_list(std::string const word, std::string const value, int const line, std::vector<grain::output_t> & output_list);


   //-------------------------------------------------------------------------
   // Function protypes for functions inside: outputfunctions.cpp
   //-------------------------------------------------------------------------
   void time(std::ostream& stream,bool header);
   void real_time(std::ostream& stream,bool header);
   void temperature(std::ostream& stream,bool header);
   void Happ(std::ostream& stream,bool header);
   void Hvec(std::ostream& stream,bool header);
   void mvec(std::ostream& stream,bool header);
   void magm(std::ostream& stream,bool header);
   void mean_magm(std::ostream& stream,bool header);
   void mat_mvec(std::ostream& stream,bool header);
   void mat_mean_magm(std::ostream& stream,bool header);
   void mdoth(std::ostream& stream,bool header);
   void systorque(std::ostream& stream,bool header);
   void mean_systorque(std::ostream& stream,bool header);
   void constraint_phi(std::ostream& stream,bool header);
   void constraint_theta(std::ostream& stream,bool header);
   void material_constraint_phi(std::ostream& stream,bool header);
   void material_constraint_theta(std::ostream& stream,bool header);
   void material_mean_systorque(std::ostream& stream,bool header);
   void material_torque(std::ostream& stream, bool header);
   void standard_deviation(std::ostream& stream,bool header);
   void mean_system_susceptibility(std::ostream& stream,bool header);
   void system_binder_cumulant(std::ostream& stream,bool header);
   void phonon_temperature(std::ostream& stream,bool header);
   void material_temperature(std::ostream& stream,bool header);
   void material_applied_field_strength(std::ostream& stream,bool header);
   void material_fmr_field_strength(std::ostream& stream,bool header);
   void mat_mdoth(std::ostream& stream,bool header);
   void total_energy(std::ostream& stream,bool header);
   void mean_total_energy(std::ostream& stream,bool header);
   void total_anisotropy_energy(std::ostream& stream,bool header);
   void mean_total_anisotropy_energy(std::ostream& stream,bool header);
   void lfa_ms(std::ostream& stream,bool header);
   void x_track_pos(std::ostream& stream,bool header);
   void z_track_pos(std::ostream& stream,bool header);
   void fractional_electric_field_strength(std::ostream& stream,bool header);
   //void total_cubic_anisotropy_energy(std::ostream& stream,bool header);
   //void mean_total_cubic_anisotropy_energy(std::ostream& stream,bool header);
   //void total_surface_anisotropy_energy(std::ostream& stream,bool header);
   //void mean_total_surface_anisotropy_energy(std::ostream& stream,bool header);
   void total_exchange_energy(std::ostream& stream,bool header);
   void mean_total_exchange_energy(std::ostream& stream,bool header);
   void total_applied_field_energy(std::ostream& stream,bool header);
   void mean_total_applied_field_energy(std::ostream& stream,bool header);
   void total_magnetostatic_energy(std::ostream& stream,bool header);
   void mean_total_magnetostatic_energy(std::ostream& stream,bool header);
   //void total_so_anisotropy_energy(std::ostream& stream,bool header);
   //void mean_total_so_anisotropy_energy(std::ostream& stream,bool header);
   void height_mvec(std::ostream& stream,bool header);
   void material_height_mvec(std::ostream& stream,bool header);
   void height_mvec_actual(std::ostream& stream,bool header);
   void material_height_mvec_actual(std::ostream& stream,bool header);
   void fmr_field_strength(std::ostream& stream,bool header);
	void mean_mvec(std::ostream& stream,bool header);
	void mat_mean_mvec(std::ostream& stream,bool header);
   void mean_material_susceptibility(std::ostream& stream,bool header);
   void material_binder_cumulant(std::ostream& stream,bool header);
   void mean_height_magnetisation_length(std::ostream& stream,bool header);
   void mean_height_magnetisation(std::ostream& stream,bool header);

   void MPITimings(std::ostream& stream,bool header);
   void mean_system_specific_heat(std::ostream& stream,bool header);
   void mean_material_specific_heat(std::ostream& stream,bool header);
   void material_total_energy(std::ostream& stream,bool header);
   void material_mean_total_energy(std::ostream& stream,bool header);
   void resistance(std::ostream& stream, bool header);
   void current(std::ostream& stream, bool header);
   void domain_wall_position(std::ostream& stream,bool header);
   void MRresistance(std::ostream& stream, bool header);

   //-------------------------------------------------------------------------
   // Funciton protypes for functions inside: datalog.cpp
   //-------------------------------------------------------------------------
   void data();
   void zLogTsInit(std::string tmp);



}

#endif //VIO_INTERNAL_H_
