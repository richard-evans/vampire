//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2014. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Vampire headers
#include "anisotropy.hpp"
#include "atoms.hpp"
#include "gpu.hpp"
#include "integrate.hpp"
#include "cuda.hpp"
#include "dipole.hpp"
#include "environment.hpp"
#include "errors.hpp"
#include "exchange.hpp"
#include "material.hpp"
#include "micromagnetic.hpp"
#include "mgpu.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "spintorque.hpp"
#include "spintransport.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

namespace gpu{

   //-------------------------------------------------------------------------------
   // Function to stop with an error message if the simulation needs features the
   // GPU version does not calculate, which would otherwise be silently missing
   //-------------------------------------------------------------------------------
   void check_supported_features(){

      if(err::check==true) std::cout << "gpu::check_supported_features has been called" << std::endl;

      std::vector<std::string> unsupported;

      // integrators: llg-heun runs on the GPU, and monte-carlo on a single GPU. With any
      // other integrator the host spins would be overwritten by the (never updated)
      // device spins at every statistics update.
      if(integrate::integrator == integrate::integrator_t::monte_carlo){
         #ifdef CUDA_MULTIGPU
            unsupported.push_back("the monte-carlo integrator in the multi-GPU version");
         #endif
         // the Monte Carlo energy includes exchange, uniaxial anisotropy and the global applied field
         bool cubic = false;
         for(int m = 0; m < mp::num_materials; m++) cubic = cubic || anisotropy::get_kc4(m) != 0.0 || anisotropy::get_kc6(m) != 0.0;
         if(cubic)                               unsupported.push_back("cubic anisotropy with the monte-carlo integrator");
         if(anisotropy::is_neel_enabled())       unsupported.push_back("Neel anisotropy with the monte-carlo integrator");
         if(anisotropy::is_rotational_enabled()) unsupported.push_back("rotational anisotropy with the monte-carlo integrator");
      }
      else if(integrate::integrator != integrate::integrator_t::llg_heun) unsupported.push_back("integrators other than llg-heun");

      // exchange
      if(::exchange::biquadratic) unsupported.push_back("biquadratic exchange");
      if(::exchange::four_spin)   unsupported.push_back("four spin exchange");

      // anisotropy
      if(anisotropy::is_triaxial_enabled())                    unsupported.push_back("triaxial anisotropy");
      if(anisotropy::is_lattice_enabled())                     unsupported.push_back("lattice anisotropy");
      if(anisotropy::is_cubic_rotation_enabled())              unsupported.push_back("rotated cubic anisotropy (cubic-anisotropy-direction)");
      if(anisotropy::is_biaxial_fourth_order_simple_enabled()) unsupported.push_back("fourth order biaxial anisotropy");
      if(sim::vcma_fields_enabled())                           unsupported.push_back("voltage controlled anisotropy");

      // fields
      if(sim::local_applied_field)        unsupported.push_back("material specific applied fields");
      if(sim::ext_demag)                  unsupported.push_back("the thin film demagnetising field (sim:demagnetisation-factor)");
      if(sim::local_temperature)          unsupported.push_back("material specific temperatures");
      if(::environment::enabled)          unsupported.push_back("the environment module");
      if(::micromagnetic::discretisation_type > 0) unsupported.push_back("micromagnetic and multiscale simulations");
      if(::micromagnetic::bias_magnets_enabled()) unsupported.push_back("micromagnetic bias magnets");
      if(program::program == program::program_t::localised_temperature_pulse) unsupported.push_back("the localised-temperature-pulse program");
      if(::dipole::activated)             unsupported.push_back("dipole fields");
      if(sim::spin_torque_fields_enabled()) unsupported.push_back("spin transfer and spin orbit torque fields");
      if(::st::is_enabled())              unsupported.push_back("the spin torque module");
      if(::spin_transport::is_enabled())  unsupported.push_back("spin transport");
      if(program::lagrange_multiplier_enabled)unsupported.push_back("the Lagrange multiplier");

      // statistics
      if(::stats::calculate_system_torque || ::stats::calculate_grain_torque || ::stats::calculate_material_torque)
         unsupported.push_back("torque statistics");
      if(::stats::calculate_system_spin_temp || ::stats::calculate_grain_spin_temp || ::stats::calculate_material_spin_temp)
         unsupported.push_back("spin temperature statistics");
      if(::stats::calculate_system_lattice_temp || ::stats::calculate_grain_lattice_temp || ::stats::calculate_material_lattice_temp)
         unsupported.push_back("lattice temperature statistics");
      if(::stats::calculate_system_spin_length || ::stats::calculate_material_spin_length || ::stats::calculate_height_spin_length)
         unsupported.push_back("spin length statistics");
      if(::stats::calculate_system_binder_cumulant || ::stats::calculate_material_binder_cumulant)
         unsupported.push_back("binder cumulant statistics");
      if(::stats::calculate_system_sld_energy || ::stats::calculate_grain_sld_energy || ::stats::calculate_material_sld_energy)
         unsupported.push_back("spin-lattice energy statistics");

      // atoms that are not integrated (gpu:skip-non-magnetic-atoms) must not enter
      // the fields or the statistics of the other atoms
      if(gpu::skip_non_magnetic_atoms){
         if(::stats::calculate_system_energy || ::stats::calculate_grain_energy || ::stats::calculate_material_energy)
            unsupported.push_back("energy statistics with gpu:skip-non-magnetic-atoms (they include the atoms that are not integrated)");
         if(integrate::integrator == integrate::integrator_t::monte_carlo)
            unsupported.push_back("gpu:skip-non-magnetic-atoms with the monte-carlo integrator");
         if(::exchange::get_exchange_type() != 0){
            unsupported.push_back("gpu:skip-non-magnetic-atoms with vectorial or tensorial exchange");
         }
         else{
            bool coupled = false;
            for(int atom = 0; atom < atoms::num_atoms && !coupled; atom++){
               const bool keep_i = gpu::skipped_material(atoms::type_array[atom]);
               for(int nn = atoms::neighbour_list_start_index[atom]; nn <= atoms::neighbour_list_end_index[atom]; nn++){
                  const bool keep_j = gpu::skipped_material(atoms::type_array[atoms::neighbour_list_array[nn]]);
                  if((keep_i || keep_j) && std::fabs(atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij) > gpu::exchange_tolerance){
                     coupled = true;
                     break;
                  }
               }
            }
            if(coupled) unsupported.push_back("gpu:skip-non-magnetic-atoms with exchange interactions of non-magnetic (keep) atoms");
         }
      }

      if(unsupported.empty()) return;

      std::string message = "Error: GPU acceleration does not support the following features of this simulation:\n";
      for(size_t i = 0; i < unsupported.size(); i++) message += "\t" + unsupported[i] + "\n";
      message += "Run without GPU acceleration (gpu:acceleration = false) to use them.";
      err::fatal(message);

   }

   //-------------------------------------------------------------------------------
   // Function to return true for a material whose atoms are not integrated
   // (non-magnetic = keep with gpu:skip-non-magnetic-atoms)
   //-------------------------------------------------------------------------------
   bool skipped_material(const int material){
      return gpu::skip_non_magnetic_atoms && mp::material[material].non_magnetic == 2;
   }

   //-------------------------------------------------------------------------------
   // Function to call correct initialization function depending on the backend (cuda)
   //-------------------------------------------------------------------------------
   void initialize(){

      // stop if the simulation needs features the GPU version does not calculate
      check_supported_features();

      if(gpu::skip_non_magnetic_atoms){
         int skipped = 0;
         for(int atom = 0; atom < atoms::num_atoms; atom++) skipped += gpu::skipped_material(atoms::type_array[atom]);
         zlog << zTs() << "Not integrating " << skipped << " of " << atoms::num_atoms << " atoms of non-magnetic (keep) materials on rank " << vmpi::my_rank << " (gpu:skip-non-magnetic-atoms)" << std::endl;
      }

      // multi-GPU version: select the GPU of this rank first
      #ifdef CUDA_MULTIGPU
         mgpu::initialize();
      #endif

      #ifdef CUDA
         initialized = vcuda::initialize(gpu::cpu_stats);
      #endif

      // Check for no initialization
      if(initialized == false){
         std::cerr << "Error: no support for GPU acceleration, please recompile using -D CUDA compile" << std::endl;
         std::cerr << "time flag with cuda development libraries v6.0 or greater." << std::endl;
         std::cerr << "Aborting program." << std::endl;
         zlog << zTs() << "Error: no support for GPU acceleration, please recompile using -D CUDA compile time flag with cuda development libraries v6.0 or greater." << std::endl;
         err::vexit();
      }

      #ifdef CUDA_MULTIGPU
         mgpu::initialize_device();
      #endif

      return;
   }


   //-------------------------------------------------------------------------------
   // Function to call correct initialization function depending on the backend (cuda)
   //-------------------------------------------------------------------------------
   void initialize_dipole(){

      // check for gpu functions and data initialised
      if( !gpu::initialized ) return;

      #ifdef CUDA
         vcuda::initialize_dipole();
      #endif

      return;
   }


} // end of namespace gpu
