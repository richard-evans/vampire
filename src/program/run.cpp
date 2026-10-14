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

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "errors.hpp"
#include "program.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// program module headers
#include "internal.hpp"

namespace program{

   //---------------------------------------------------------------------------
   // Function to run the program selected with the sim:program keyword
   //---------------------------------------------------------------------------
   void run(){

      switch(program::program){
         case program::program_t::benchmark:
            if(vmpi::my_rank==0){
               std::cout << "Benchmark..." << std::endl;
               zlog << "Benchmark..." << std::endl;
            }
            program::bmark();
            break;

         case program::program_t::time_series:
            if(vmpi::my_rank==0){
               std::cout << "Time-Series..." << std::endl;
               zlog << "Time-Series..." << std::endl;
            }
            program::time_series();
            break;

         case program::program_t::hysteresis_loop:
            if(vmpi::my_rank==0){
               std::cout << "Hysteresis-Loop..." << std::endl;
               zlog << "Hysteresis-Loop..." << std::endl;
            }
            program::hysteresis();
            break;

         case program::program_t::static_hysteresis_loop:
            if(vmpi::my_rank==0){
               std::cout << "Static-Hysteresis-Loop..." << std::endl;
               zlog << "Static-Hysteresis-Loop..." << std::endl;
            }
            program::static_hysteresis();
            break;

         case program::program_t::curie_temperature:
            if(vmpi::my_rank==0){
               std::cout << "Curie-Temperature..." << std::endl;
               zlog << "Curie-Temperature..." << std::endl;
            }
            program::curie_temperature();
            break;

         case program::program_t::field_cool:
            if(vmpi::my_rank==0){
               std::cout << "Field-Cool..." << std::endl;
               zlog << "Field-Cool..." << std::endl;
            }
            program::field_cool();
            break;

         case program::program_t::laser_pulse:
            if(vmpi::my_rank==0){
               std::cout << "Temperature-Pulse..." << std::endl;
               zlog << "Temperature-Pulse..." << std::endl;
            }
            program::temperature_pulse();
            break;

         case program::program_t::hamr_simulation:
            if(vmpi::my_rank==0){
               std::cout << "HAMR-Simulation..." << std::endl;
               zlog << "HAMR-Simulation..." << std::endl;
            }
            program::hamr();
            break;

         case program::program_t::cmc_anisotropy:
            if(vmpi::my_rank==0){
               std::cout << "CMC-Anisotropy..." << std::endl;
               zlog << "CMC-Anisotropy..." << std::endl;
            }
            program::cmc_anisotropy();
            break;

         case program::program_t::hybrid_cmc:
            if(vmpi::my_rank==0){
               std::cout << "Hybrid-CMC..." << std::endl;
               zlog << "Hybrid-CMC..." << std::endl;
            }
            program::hybrid_cmc();
            break;

         case program::program_t::reverse_hybrid_cmc:
            if(vmpi::my_rank==0){
               std::cout << "Reverse-Hybrid-CMC..." << std::endl;
               zlog << "Reverse-Hybrid-CMC..." << std::endl;
            }
            program::reverse_hybrid_cmc();
            break;

         case program::program_t::lagrange_multiplier:
            if(vmpi::my_rank==0){
               std::cout << "LaGrange-Multiplier..." << std::endl;
               zlog << "LaGrange-Multiplier..." << std::endl;
            }
            program::lagrange_multiplier();
            break;

         case program::program_t::partial_hysteresis_loop:
            if(vmpi::my_rank==0){
               std::cout << "partial-hysteresis-loop..." << std::endl;
               zlog << "partial-hysteresis-loop..." << std::endl;
            }
            program::partial_hysteresis_loop();
            break;

         case program::program_t::localised_temperature_pulse:
            if(vmpi::my_rank==0){
               std::cout << "localised-temperature-pulse..." << std::endl;
               zlog << "localised-temperature-pulse..." << std::endl;
            }
            program::localised_temperature_pulse();
            break;

         case program::program_t::effective_damping:
            if(vmpi::my_rank==0){
               std::cout << "effective-damping..." << std::endl;
               zlog << "effective-damping..." << std::endl;
            }
            program::effective_damping();
            break;

         case program::program_t::fmr:
            if(vmpi::my_rank==0){
               std::cout << "fmr..." << std::endl;
               zlog << "fmr..." << std::endl;
            }
            program::fmr();
            break;

         case program::program_t::localised_field_cool:
            if(vmpi::my_rank==0){
               std::cout << "localised-field-cool..." << std::endl;
               zlog << "localised-field-cool..." << std::endl;
            }
            program::local_field_cool();
            break;

         case program::program_t::electrical_pulse:
            if(vmpi::my_rank==0){
               std::cout << "electrical-pulse..." << std::endl;
               zlog << "electrical-pulse..." << std::endl;
            }
            program::electrical_pulse();
            break;
         case program::program_t::field_pulse:
            if(vmpi::my_rank==0){
               std::cout << "field-pulse..." << std::endl;
               zlog << "field-pulse..." << std::endl;
            }
            program::field_pulse();
            break;
         case program::program_t::laser_electrical_pulse:
            if(vmpi::my_rank==0){
               std::cout << "laser-electrical-pulse..." << std::endl;
               zlog << "laser-electrical-pulse..." << std::endl;
            }
            program::laser_electrical_pulse();
            break;

         case program::program_t::diagnostic_boltzmann:
            if(vmpi::my_rank==0){
               std::cout << "Diagnostic-Boltzmann..." << std::endl;
               zlog << "Diagnostic-Boltzmann..." << std::endl;
            }
            program::boltzmann_dist();
            break;
         case program::program_t::setting:
            if(vmpi::my_rank==0){
               std::cout << "Setting..." << std::endl;
               zlog << "Setting..." << std::endl;
            }
            program::setting_process();
            break;
         //------------------------------------------------------------------------
         case program::program_t::domain_wall:
            if(vmpi::my_rank==0){
               std::cout << "Domain walls..." << std::endl;
               zlog << "Domain walls..." << std::endl;
            }
            program::domain_wall();
            break;
         //------------------------------------------------------------------------
         case program::program_t::exchange_stiffness:
            if(vmpi::my_rank==0){
               std::cout << "exchange stiffness..." << std::endl;
               zlog << "exchange stiffness..." << std::endl;
            }
            program::exchange_stiffness();
            break;
         //------------------------------------------------------------------------
         case program::program_t::mm_a_calculation:
            if(vmpi::my_rank==0){
               std::cout << "mm-A-calculation..." << std::endl;
               zlog << "mm-A-calculation..." << std::endl;
            }
            program::mm_A_calculation();
            break;
         //------------------------------------------------------------------------
         case program::program_t::field_sweep:
            if(vmpi::my_rank==0){
               std::cout << "field-sweep..." << std::endl;
               zlog << "field-sweep..." << std::endl;
            }
            program::field_sweep();
         break;
         //------------------------------------------------------------------------
         case program::program_t::disk_tracks:
            if(vmpi::my_rank==0){
               std::cout << "Tracks..." << std::endl;
               zlog << "Tracks..." << std::endl;
            }
            program::tracks();
            break;
         //------------------------------------------------------------------------
         case program::program_t::diagnostic_boltzmann_micromagnetic_llg:
            if(vmpi::my_rank==0){
               std::cout << "diagnostic-boltzmann-micromganetic-llg..." << std::endl;
               zlog << "diagnostic-boltzmann-micromganetic-llg..." << std::endl;
            }
            program::boltzmann_dist_micromagnetic_llg();
            break;
         // JRH call spin_waves program -------------------------------------------
         case program::program_t::spin_waves:
            if(vmpi::my_rank==0){
               std::cout << "Spin-waves..." << std::endl;
               zlog << "Spin-waves..." << std::endl;
            }
            program::spin_waves();
            break;
         //------------------------------------------------------------------------
         case program::program_t::hamr_cool:
            if(vmpi::my_rank==0){
               std::cout << "HAMR cool..." << std::endl;
               zlog << "HAMR cool..." << std::endl;
            }
            program::hamr_cool();
            break;

         //------------------------------------------------------------------------
         default:{
            std::cerr << "Unknown Internal Program ID "<< static_cast<int>(program::program) << " requested, exiting" << std::endl;
            zlog << "Unknown Internal Program ID "<< static_cast<int>(program::program) << " requested, exiting" << std::endl;
            err::vexit();
            }
      }

      return;

   }

} // end of program namespace
