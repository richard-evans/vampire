#ifndef SIM_INTERNAL_H_
#define SIM_INTERNAL_H_
//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2025. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Vampire headers
#include "vtypes.hpp"

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// simulation methods implementation. These functions should
// not be accessed outside of the simulate module.
//---------------------------------------------------------------------

namespace sim{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Internal data types used for simulation module
      //-----------------------------------------------------------------------------
      struct mp_t{
         vtype::set_double_t stt_asm; // spin tranfer torque asymmetry
         vtype::set_double_t stt_rj;  // spin tranfer relaxation torque
         vtype::set_double_t stt_pj;  // spin transfer precession torque
         vtype::set_int_t stt_pm;  // spin transfer polarisation material
         vtype::vec_t stt_pv;  // spin transfer polarisation vector x,y,z
         vtype::set_double_t sot_asm; // spin orbit torque asymmetry
         vtype::set_double_t sot_rj;  // spin orbit relaxation torque
         vtype::set_double_t sot_pj;  // spin orbit precession torque
         vtype::set_int_t sot_pm;  // spin orbit polarisation material
         vtype::vec_t sot_pv;  // spin orbit polarisation vector x,y,z
         vtype::set_double_t vcmak;   // voltage controlled anisotropy coefficient
         vtype::set_double_t lsf_second_order_coefficient; // Second order LSF coefficient
         vtype::set_double_t lsf_fourth_order_coefficient; // Fourth order LSF coefficient
         vtype::set_double_t lsf_sixth_order_coefficient; // Sixth order LSF coefficient
      };

      //-----------------------------------------------------------------------------
      // Internal shared variables used for the simulation
      //-----------------------------------------------------------------------------
      extern bool enable_spin_torque_fields; // flag to enable spin torque fields
      extern bool enable_vcma_fields;        // flag to enable voltage-controlled anisotropy fields
      extern bool enable_local_stt_polarizers; // flag to enable localised stt polarization vectors
      extern bool enable_local_sot_polarizers; // flag to enable localised sot polarization vectors

      extern std::vector<sim::internal::mp_t> mp; // array of material properties

      extern std::vector<double> stt_asm; // array of spin transfer torque asymmetry
      extern std::vector<double> stt_rj; // array of adiabatic spin torques
      extern std::vector<double> stt_pj; // array of non-adiabatic spin torques
      extern std::vector<int> stt_pm; // array of polarisation material for spin torques
      extern vtype::vec_t stt_polarization_unit_vector; // stt spin polarization direction
      extern std::vector<vtype::vec_t>  stt_material_polarization_unit_vector; // array of stt spin polarization direction

      extern std::vector<double> sot_asm; // array of spin orbit torque asymmetry
      extern std::vector<double> sot_rj; // array of adiabatic spin torques
      extern std::vector<double> sot_pj; // array of non-adiabatic spin torques
      extern std::vector<int> sot_pm; // array of polarisation material for spin torques
      extern vtype::vec_t sot_polarization_unit_vector; // sot spin polarization direction
      extern std::vector<vtype::vec_t>  sot_material_polarization_unit_vector; // array of sot spin polarization direction

      extern std::vector<double> vcmak;   // voltage controlled anisotropy coefficient

      extern std::vector<double> lsf_second_order_coefficient;
      extern std::vector<double> lsf_fourth_order_coefficient; // LSF coefficients
      extern std::vector<double> lsf_sixth_order_coefficient;

      // shared Functions
      void llg_quantum_step();

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      extern void initialize_modules();
      extern void increment_time();
      extern void lsf_step();
      extern void lsf_rk4_step();

   } // end of internal namespace
} // end of sim namespace

#endif //SIM_INTERNAL_H_
