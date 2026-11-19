//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "sim.hpp"

// Internal sim header
#include "internal.hpp"

namespace sim{
   //----------------------------------------------------------------------------
   // Shared variables used with main vampire code
   //---------------------------------------------------------------------------


   int num_monte_carlo_preconditioning_steps(0);

   uint64_t time         = 0; // time step counter
   uint64_t total_time   = 10000; // total time steps (non-loop code)
   uint64_t loop_time    = 10000; // loop time steps (hysteresis/temperature loops)
   uint64_t partial_time = 1000; // same as time-step-increment
   uint64_t equilibration_time = 0; // equilibration time steps

   namespace internal{

      //----------------------------------------------------------------------------
      // Shared variables used within sim module
      //----------------------------------------------------------------------------
      bool enable_spin_torque_fields = false; // flag to enable spin torque fields
      bool enable_local_stt_polarizers = false; // flag to enable localised stt polarization vectors
      bool enable_local_sot_polarizers = false; // flag to enable localised sot polarization vectors

      bool enable_vcma_fields        = false; // flag to enable voltage-controlled anisotropy fields

      std::vector<sim::internal::mp_t> mp; // array of material properties

      std::vector<double> stt_asm; // array of spin transfer torque asymmetry
      std::vector<double> stt_rj;  // array of adiabatic spin torques
      std::vector<double> stt_pj;  // array of non-adiabatic spin torques
      std::vector<int> stt_pm;     // array of spin orbit torque polarizer materials
      vtype::vec_t stt_polarization_unit_vector; // stt spin polarization direction
      std::vector<vtype::vec_t> stt_material_polarization_unit_vector; // stt spin polarization direction for each material

      std::vector<double> sot_asm; // array of spin orbit torque asymmetry
      std::vector<double> sot_rj;  // array of adiabatic spin torques
      std::vector<double> sot_pj;  // array of non-adiabatic spin torques
      std::vector<int> sot_pm;     // array of spin orbit torque polarizer materials
      vtype::vec_t sot_polarization_unit_vector; // sot spin polarization direction
      std::vector<vtype::vec_t>  sot_material_polarization_unit_vector; // sot spin polarization direction for each material

      std::vector<double> vcmak;   // voltage controlled anisotropy coefficient


   } // end of internal namespace

   //------------------------------------------------------------------------
   // getter functions to give access to internal variables
   //------------------------------------------------------------------------
   std::vector<double> get_stt_polarization_unit_vector(){
      std::vector<double> sttpv(3);
      sttpv[0] = sim::internal::stt_polarization_unit_vector.x;
      sttpv[1] = sim::internal::stt_polarization_unit_vector.y;
      sttpv[2] = sim::internal::stt_polarization_unit_vector.z;
      return sttpv;
   }

   std::vector<double> get_stt_rj(){
      return sim::internal::stt_rj;
   }

   std::vector<double> get_stt_pj(){
      return sim::internal::stt_pj;
   }

   bool spin_torque_fields_enabled(){
      return sim::internal::enable_spin_torque_fields;
   }

   bool vcma_fields_enabled(){
      return sim::internal::enable_vcma_fields;
   }

} // end of sim namespace
