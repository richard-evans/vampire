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
   int match_create(std::string const, std::string const, std::string const, int const);
   int match_dimension(std::string const, std::string const, std::string const, int const);
   int match_sim(std::string const, std::string const, std::string const, int const);
   int match_vout_list(std::string const, std::string const, int const, std::vector<unsigned int> &);
   int match_vout_grain_list(std::string const, std::string const, int const, std::vector<unsigned int> &);
   int match_material(std::string const, std::string const, std::string const, int const, int const, int const, std::string const, std::string const);
   int match_config(std::string const, std::string const, std::string const, int const);

   //----------------------------------------------------------------------------------
   //Function protypes for functions inside: read.cpp
   //----------------------------------------------------------------------------------
   int read(std::string const filename);
   int read_mat_file(std::string const, int const);

}

namespace vout{
   //-------------------------------------------------------------------------
   // Funciton protypes for functions inside: outputfunctions.cpp
   //-------------------------------------------------------------------------
   void time(std::ostream& stream);
   void real_time(std::ostream& stream);
   void temperature(std::ostream& stream);
   void Happ(std::ostream& stream);
   void Hvec(std::ostream& stream);
   void mvec(std::ostream& stream);
   void magm(std::ostream& stream);
   void mean_magm(std::ostream& stream);
   void mat_mvec(std::ostream& stream);
   void mat_mean_magm(std::ostream& stream);
   void grain_mvec(std::ostream& stream);
   void grain_magm(std::ostream& stream);
   void mdoth(std::ostream& stream);
   void grain_mat_mvec(std::ostream& stream);
   void systorque(std::ostream& stream);
   void mean_systorque(std::ostream& stream);
   void constraint_phi(std::ostream& stream);
   void constraint_theta(std::ostream& stream);
   void material_constraint_phi(std::ostream& stream);
   void material_constraint_theta(std::ostream& stream);
   void material_mean_systorque(std::ostream& stream);
   void mean_system_susceptibility(std::ostream& stream);
   void phonon_temperature(std::ostream& stream);
   void material_temperature(std::ostream& stream);
   void material_applied_field_strength(std::ostream& stream);
   void material_fmr_field_strength(std::ostream& stream);
   void mat_mdoth(std::ostream& stream);
   void total_energy(std::ostream& stream);
   void mean_total_energy(std::ostream& stream);
   void total_anisotropy_energy(std::ostream& stream);
   void mean_total_anisotropy_energy(std::ostream& stream);
   void total_cubic_anisotropy_energy(std::ostream& stream);
   void mean_total_cubic_anisotropy_energy(std::ostream& stream);
   void total_surface_anisotropy_energy(std::ostream& stream);
   void mean_total_surface_anisotropy_energy(std::ostream& stream);
   void total_exchange_energy(std::ostream& stream);
   void mean_total_exchange_energy(std::ostream& stream);
   void total_applied_field_energy(std::ostream& stream);
   void mean_total_applied_field_energy(std::ostream& stream);
   void total_magnetostatic_energy(std::ostream& stream);
   void mean_total_magnetostatic_energy(std::ostream& stream);
   void total_so_anisotropy_energy(std::ostream& stream);
   void mean_total_so_anisotropy_energy(std::ostream& stream);
   void height_mvec(std::ostream& stream);
   void material_height_mvec(std::ostream& stream);
   void height_mvec_actual(std::ostream& stream);
   void material_height_mvec_actual(std::ostream& stream);
   void fmr_field_strength(std::ostream& stream);
	void mean_mvec(std::ostream& stream);
	void mat_mean_mvec(std::ostream& stream);

   void MPITimings(std::ostream& stream);

   //-------------------------------------------------------------------------
   // Funciton protypes for functions inside: datalog.cpp
   //-------------------------------------------------------------------------
   void data();
   void zLogTsInit(std::string tmp);

}

#endif //VIO_INTERNAL_H_
