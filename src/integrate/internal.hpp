//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef INTEGRATE_INTERNAL_H_
#define INTEGRATE_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the integrate module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <cstdint>
#include <vector>

// Vampire headers
#include "integrate.hpp"
#include "vtypes.hpp"

//-----------------------------------------------------------------------------
// Persistent integration arrays (previously hdr/LLG.hpp, hdr/lsf.hpp and
// hdr/lsf_rk4.hpp). LSF arrays (c) David R Papp 2024.
//-----------------------------------------------------------------------------
namespace LLG_arrays{

//==========================================================
// Namespace to store persistant LLG integration arrays
//==========================================================

   extern std::vector <double> x_euler_array;
   extern std::vector <double> y_euler_array;
   extern std::vector <double> z_euler_array;

   extern std::vector <double> x_heun_array;
   extern std::vector <double> y_heun_array;
   extern std::vector <double> z_heun_array;

   extern std::vector <double> x_spin_storage_array;
   extern std::vector <double> y_spin_storage_array;
   extern std::vector <double> z_spin_storage_array;

   extern std::vector <double> x_initial_spin_array;
   extern std::vector <double> y_initial_spin_array;
   extern std::vector <double> z_initial_spin_array;

   extern bool LLG_set;

} // end of LLG_arrays namespace

namespace LSF_arrays{

//==========================================================
// Namespace to store persistant LSF integration arrays
//==========================================================

   extern std::vector<double> x_lsf_array;
   extern std::vector<double> y_lsf_array;
   extern std::vector<double> z_lsf_array;

   extern std::vector<double> x_euler_array;
   extern std::vector<double> y_euler_array;
   extern std::vector<double> z_euler_array;

   extern std::vector<double> x_heun_array;
   extern std::vector<double> y_heun_array;
   extern std::vector<double> z_heun_array;

   extern std::vector<double> x_spin_storage_array;
   extern std::vector<double> y_spin_storage_array;
   extern std::vector<double> z_spin_storage_array;

   extern std::vector<double> x_initial_spin_array;
   extern std::vector<double> y_initial_spin_array;
   extern std::vector<double> z_initial_spin_array;

   extern bool LSF_set;

   extern std::vector<double> mod_S;

   extern std::vector<double> tx;
   extern std::vector<double> ty;
   extern std::vector<double> tz;

} // end of LSF_arrays namespace

namespace LSF_RK4_arrays{

//==========================================================
// Namespace to store persistant LSF-RK4 integration arrays
//==========================================================

   extern std::vector<double> x_lsf_array;
   extern std::vector<double> y_lsf_array;
   extern std::vector<double> z_lsf_array;

   extern std::vector<double> x_initial_spin_array;
   extern std::vector<double> y_initial_spin_array;
   extern std::vector<double> z_initial_spin_array;

   extern std::vector<double> x_k1_array;
   extern std::vector<double> y_k1_array;
   extern std::vector<double> z_k1_array;

   extern std::vector<double> x_k2_array;
   extern std::vector<double> y_k2_array;
   extern std::vector<double> z_k2_array;

   extern std::vector<double> x_k3_array;
   extern std::vector<double> y_k3_array;
   extern std::vector<double> z_k3_array;

   extern std::vector<double> x_k4_array;
   extern std::vector<double> y_k4_array;
   extern std::vector<double> z_k4_array;

   extern bool LSF_RK4_set;

   extern std::vector<double> mod_S;

   extern std::vector<double> tx;
   extern std::vector<double> ty;
   extern std::vector<double> tz;

} // end of LSF_RK4_arrays namespace

namespace LLG_rotation_arrays{

//==========================================================
// Namespace to store persistant arrays for rotation based
// LLG integrators (RKMK2, RKMK4, Depondt-Mertens)
//==========================================================

   extern std::vector<double> x_initial_spin_array;
   extern std::vector<double> y_initial_spin_array;
   extern std::vector<double> z_initial_spin_array;

   extern std::vector<double> x_k1_array;
   extern std::vector<double> y_k1_array;
   extern std::vector<double> z_k1_array;

   extern std::vector<double> x_k2_array;
   extern std::vector<double> y_k2_array;
   extern std::vector<double> z_k2_array;

   extern std::vector<double> x_k3_array;
   extern std::vector<double> y_k3_array;
   extern std::vector<double> z_k3_array;

   extern std::vector<double> x_spin_storage_array;
   extern std::vector<double> y_spin_storage_array;
   extern std::vector<double> z_spin_storage_array;

   extern bool LLG_rotation_set;

} // end of LLG_rotation_arrays namespace

namespace integrate{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-----------------------------------------------------------------------------
      // internal materials class for storing material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

          private:

          public:

             //------------------------------
             // material parameter variables
             //------------------------------
             vtype::set_double_t lsf_second_order_coefficient; // Second order LSF coefficient
             vtype::set_double_t lsf_fourth_order_coefficient; // Fourth order LSF coefficient
             vtype::set_double_t lsf_sixth_order_coefficient;  // Sixth order LSF coefficient

       }; // end of internal::mp class

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      extern bool enabled; // bool to enable module

      extern std::vector<internal::mp_t> mp; // array of material properties

      extern std::vector<double> lsf_second_order_coefficient; // LSF Landau coefficients per material (J)
      extern std::vector<double> lsf_fourth_order_coefficient;
      extern std::vector<double> lsf_sixth_order_coefficient;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      // integrator dispatch
      void integrate_serial(uint64_t n_steps);
      int integrate_mpi(uint64_t n_steps);
      void multiscale_simulation_steps(int n_steps);
      void increment_time();

      // LLG Heun and midpoint integrators
      int LLGinit();
      int LLG_Heun();
      int LLG_Heun_mpi();
      int LLG_Heun_cuda();
      int LLG_Midpoint();
      int LLG_Midpoint_mpi();
      int LLG_Midpoint_cuda();

      // LLG integrator with quantum thermostat
      int LLGQinit();
      void llg_quantum_step();

      // rotation based LLG integrators (ported from JAMS)
      void initialize_rotation_integrators();
      int LLG_rotation_init();
      void llg_rotation_noise_half_step(const int start_index, const int end_index);
      void llg_rotation_copy_storage_to_spins(const int start_index, const int end_index);
      void llg_rotation_stage(void (*stage)(const int, const int), const bool calculate_external_fields);
      void llg_rotation_stage_mpi(void (*stage)(const int, const int), const bool calculate_external_fields);

      // Runge-Kutta-Munthe-Kaas integrators
      void rkmk_stage_1(const int start_index, const int end_index);
      void rkmk2_stage_2(const int start_index, const int end_index);
      void rkmk4_stage_2(const int start_index, const int end_index);
      void rkmk4_stage_3(const int start_index, const int end_index);
      void rkmk4_stage_4(const int start_index, const int end_index);
      void llg_rkmk2_step();
      void llg_rkmk4_step();
      void llg_rkmk2_mpi_step();
      void llg_rkmk4_mpi_step();

      // Depondt-Mertens integrator
      void dm_predictor_stage(const int start_index, const int end_index);
      void dm_corrector_stage(const int start_index, const int end_index);
      void llg_dm_step();
      void llg_dm_mpi_step();

      // longitudinal spin fluctuation integrators
      int LSFinit();
      void calculate_lsf_magnetic_field(const int start_index, const int end_index);
      void lsf_step();
      int LSF_mpi();
      int LSF_RK4_init();
      void calculate_lsf_rk4_magnetic_field(const int start_index, const int end_index);
      void lsf_rk4_step();
      int LSF_RK4_mpi();

   } // end of internal namespace

} // end of integrate namespace

#endif //INTEGRATE_INTERNAL_H_
