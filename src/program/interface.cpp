//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "program.hpp"
#include "sim.hpp"
#include "errors.hpp"
#include "vio.hpp"

// program module headers
#include "internal.hpp"

namespace program{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for program module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false - note use of sim: module name for semantic reasons
      // program variables are generally related to simulations but a logically separated for ease of the programmer
      std::string prefix="sim";
      if(key!=prefix) return false;
      //-------------------------------------------------------------------
      std::string test="program";
      if(word==test){
         test="benchmark";
         if(value==test){
            program::program = program_t::benchmark;
            return true;
         }
         test="time-series";
         if(value==test){
            program::program = program_t::time_series;
            return true;
         }
         test="hysteresis-loop";
         if(value==test){
            program::program = program_t::hysteresis_loop;
            return true;
         }
         test="static-hysteresis-loop";
         if(value==test){
            program::program = program_t::static_hysteresis_loop;
            return true;
         }
         test="curie-temperature";
         if(value==test){
            program::program = program_t::curie_temperature;
            return true;
         }
         test="field-cool";
         if(value==test){
            program::program = program_t::field_cool;
            return true;
         }
         test="laser-pulse";
         if(value==test){
            program::program = program_t::laser_pulse;
            return true;
         }
         test="hamr-simulation";
         if(value==test){
            program::program = program_t::hamr_simulation;
            return true;
         }
         test="cmc-anisotropy";
         if(value==test){
            program::program = program_t::cmc_anisotropy;
            return true;
         }
         test="hybrid-cmc";
         if(value==test){
            program::program = program_t::hybrid_cmc;
            return true;
         }
         test="reverse-hybrid-cmc";
         if(value==test){
            program::program = program_t::reverse_hybrid_cmc;
            return true;
         }
         test="LaGrange-Multiplier";
         if(value==test){
            program::program = program_t::lagrange_multiplier;
            return true;
         }
         test="partial-hysteresis-loop";
         if(value==test){
            program::program = program_t::partial_hysteresis_loop;
            return true;
         }
         test="localised-temperature-pulse";
         if(value==test){
            program::program = program_t::localised_temperature_pulse;
            return true;
         }
         test="effective-damping";
         if(value==test){
            program::program = program_t::effective_damping;
            return true;
         }
         test="fmr";
         if(value==test){
            program::program = program_t::fmr;
            return true;
         }
         test="localised-field-cool";
         if(value==test){
            program::program = program_t::localised_field_cool;
            return true;
         }
         test="electrical-pulse";
         if(value==test){
            program::program = program_t::electrical_pulse;
            return true;
         }
         test="field-pulse";
         if(value==test){
            program::program = program_t::field_pulse;
            return true;
         }
         test="laser-electrical-pulse";
         if(value==test){
            program::program = program_t::laser_electrical_pulse;
            return true;
         }
         test="diagnostic-boltzmann";
         if(value==test){
            program::program = program_t::diagnostic_boltzmann;
            return true;
         }
         test="setting";
         if(value==test){
            program::program = program_t::setting;
            return true;
         }
         test="domain-wall";
         if(value==test){
            program::program = program_t::domain_wall;
            return true;
         }
         test="exchange-stiffness";
         if(value==test){
            program::program = program_t::exchange_stiffness;
            return true;
         }
         test="mm-A-calculation";
         if(value==test){
            program::program = program_t::mm_a_calculation;
            return true;
         }
         test="field-sweep";
         if(value==test){
            program::program = program_t::field_sweep;
            return true;
         }
         test="disk-tracks";
         if(value==test){
            program::program = program_t::disk_tracks;
            return true;
         }
         test="diagnostic-boltzmann-micromagnetic-llg";
         if(value==test){
            program::program = program_t::diagnostic_boltzmann_micromagnetic_llg;
            return true;
         }
         // JRH - Program option for spinwaves
         test="spin-waves";
         if(value==test){
            program::program = program_t::spin_waves;
            return true;
         }
         // HAMR cool
         test="hamr-cool";
         if(value==test){
            program::program = program_t::hamr_cool;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cout << word << '\t' << test << std::endl;
            std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"benchmark\"" << std::endl;
            std::cerr << "\t\"cmc-anisotropy\"" << std::endl;
            std::cerr << "\t\"curie-temperature\"" << std::endl;
            std::cerr << "\t\"domain-wall\"" << std::endl;
            std::cerr << "\t\"effective-damping\"" << std::endl;
            std::cerr << "\t\"electrical-pulse\"" << std::endl;
            std::cerr << "\t\"exchange-stiffness\"" << std::endl;
            std::cerr << "\t\"field-cool\"" << std::endl;
            std::cerr << "\t\"field-pulse\"" << std::endl;
            std::cerr << "\t\"laser-electrical-pulse\"" << std::endl;
            std::cerr << "\t\"laser-pulse\"" << std::endl;
            std::cerr << "\t\"localised-field-cool\"" << std::endl;
            std::cerr << "\t\"localised-temperature-pulse\"" << std::endl;
            std::cerr << "\t\"time-series\"" << std::endl;
            std::cerr << "\t\"hysteresis-loop\"" << std::endl;
            std::cerr << "\t\"partial-hysteresis-loop\"" << std::endl;
            std::cerr << "\t\"hybrid-cmc\"" << std::endl;
            std::cerr << "\t\"reverse-hybrid-cmc\"" << std::endl;
            std::cerr << "\t\"static-hysteresis-loop\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //-------------------------------------------------------------------
      test = "electrical-pulse-time";
      if(word == test){
         double pt = atof(value.c_str()); // convert string to uint64_t
         // Test for valid range
         vin::check_for_valid_positive_value(pt, word, line, prefix, unit, "time", 0.0, 1.0,"input","0.0 - 1s");
         // save sanitized value
         program::internal::electrical_pulse_time = pt;
         return true;
      }
      //-------------------------------------------------------------------
      test = "electrical-pulse-rise-time";
      if(word == test){
         double rt = atof(value.c_str()); // convert string to uint64_t
         // Test for valid range
         vin::check_for_valid_positive_value(rt, word, line, prefix, unit, "time", 0.0, 1.0,"input","0.0 - 1s");
         // save sanitized value
         program::internal::electrical_pulse_rise_time = rt;
         return true;
      }
      //-------------------------------------------------------------------
      test = "electrical-pulse-fall-time";
      if(word == test){
         double ft = atof(value.c_str()); // convert string to uint64_t
         // Test for valid range
         vin::check_for_valid_positive_value(ft, word, line, prefix, unit, "time", 0.0, 1.0,"input","0.0 - 1s");
         // save sanitized value
         program::internal::electrical_pulse_fall_time = ft;
         return true;
      }
      //--------------------------------------------------------------------
      test = "electrical-pulse-delay";
      if(word == test){
         double pd = atof(value.c_str());
         vin::check_for_valid_positive_value(pd, word, line, prefix, unit, "time", 0.0, 1.0,"input","0.0 - 1s");
         program::internal::electrical_pulse_delay = pd;
         return true;
      }
      //--------------------------------------------------------------------
      test="exchange-stiffness-maximum-angle";
      if(word==test){
         double ma = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(ma, word, line, prefix, unit, "", 0.0, 180.1,"input","0 - 180 degrees");
         program::internal::exchange_stiffness_max_constraint_angle = ma;
         return true;
      }
      //--------------------------------------------------------------------
      test="exchange-stiffness-angle-increment";
      if(word==test){
         double ai = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(ai, word, line, prefix, unit, "", 1.0, 90.0,"input","1 - 90 degrees");
         program::internal::exchange_stiffness_delta_constraint_angle = ai;
         return true;
      }
      //--------------------------------------------------------------------
      test="effective-damping-rotation-angle";
      if(word==test){
         double ra = atof(value.c_str());
         vin::check_for_valid_value(ra, word, line, prefix, unit, "", 0.0, 180.0,"input","0 - 180 degrees");
         program::internal::effective_damping_rotation_angle = ra;
         return true;
      }
      //-------------------------------------------------------------------
      test = "field-pulse-time";
      if(word == test){
         double pt = atof(value.c_str()); // convert string to uint64_t
         // Test for valid range
         vin::check_for_valid_positive_value(pt, word, line, prefix, unit, "time", 0.0, 1.0,"input","0.0 - 1s");
         // save sanitized value
         program::internal::field_pulse_time = pt;
         return true;
      }
      //--------------------------------------------------------------------
      test="track-Ms";
      if(word==test){
         double m=atof(value.c_str());
         vin::check_for_valid_value(m, word, line, prefix, unit, "magnetisation", 0.0, 360.0,"input","0.0 - 360.0");
         program::track_Ms = m;
         return true;
      }
      //--------------------------------------------------------------------
      test="track-file";
      if(word==test){
         program::track_ms_file=true; // Save checkpoint
         return true;
      }
      //--------------------------------------------------------------------
      test="track-bit-size-x";
      if(word==test){
         double m=atof(value.c_str());
         vin::check_for_valid_value(m, word, line, prefix, unit, "length", 0.0, 10000000.0,"input","0.0 - 10,000,000.0");
         program::track_bit_size = m;
         return true;
      }
      //--------------------------------------------------------------------
      test="track-bit-size-y";
      if(word==test){
         double m=atof(value.c_str());
         vin::check_for_valid_value(m, word, line, prefix, unit, "length", 0.0, 10000000.0,"input","0.0 - 10,000,000.0");
         program::track_bit_depth = m;
         return true;
      }
      //--------------------------------------------------------------------
      test="track-bit-size-z";
      if(word==test){
         double m=atof(value.c_str());
         vin::check_for_valid_value(m, word, line, prefix, unit, "length", 0.0, 10000000.0,"input","0.0 - 10,000,000.0");
         program::track_bit_width = m;
         return true;
      }
      //--------------------------------------------------------------------
      test="track-fly-height";
      if(word==test){
         double m=atof(value.c_str());
         vin::check_for_valid_value(m, word, line, prefix, unit, "length", 0.0, 36000.0,"input","0.0 - 36,000.0");
         program::track_fly_height = m;
         return true;
      }
      //--------------------------------------------------------------------
      test="track-gap-between-bits";
      if(word==test){
         double m=atof(value.c_str());
         vin::check_for_valid_value(m, word, line, prefix, unit, "length", 0.0, 36000.0,"input","0.0 - 36,000.0");
         program::track_bit_gap = m;
         return true;
      }
      //--------------------------------------------------------------------
      test="track-gap-between-tracks";
      if(word==test){
         double m=atof(value.c_str());
         vin::check_for_valid_value(m, word, line, prefix, unit, "length", 0.0, 36000.0,"input","0.0 - 36,000.0");
         program::track_track_gap = m;
         return true;
      }
      //--------------------------------------------------------------------
      test="num-tracks";
      if(word==test){
         int n=atoi(value.c_str());
         vin::check_for_valid_int(n, word, line, prefix, 0, 100,"input","0 - 100");
         program::track_num_tracks = n;
         return true;
      }
      //--------------------------------------------------------------------
      test="num-bits-per-track";
      if(word==test){
         int n=atoi(value.c_str());
         vin::check_for_valid_int(n, word, line, prefix, 0, 100,"input","0 - 100");
         program::track_num_bits_per_track = n;
         return true;
      }
      //--------------------------------------------------------------------
      test="x-velocity";
      if(word==test){
         double m=atof(value.c_str());
         vin::check_for_valid_value(m, word, line, prefix, unit, "length", 0.0, 360.0,"input","0.0 - 360.0");
         program::cross_track_velocity = m;
         return true;
      }
      //--------------------------------------------------------------------
      test="z-velocity";
      if(word==test){
         double m=atof(value.c_str());
         vin::check_for_valid_value(m, word, line, prefix, unit, "length", 0.0, 360.0,"input","0.0 - 360.0");
         program::down_track_velocity = m;
         return true;
      }
      //--------------------------------------------------------------------
      test="initial-x-position";
      if(word==test){
         double m=atof(value.c_str());
         vin::check_for_valid_value(m, word, line, prefix, unit, "length", -100000000.0, 10000000.0,"input","-100,000,000.0 - 10,000,000.0");
         program::initial_cross_track_position = m;
         return true;
      }
      //--------------------------------------------------------------------
      test="initial-z-position";
      if(word==test){
         double m=atof(value.c_str());
         vin::check_for_valid_value(m, word, line, prefix, unit, "length", -100000000.0, 10000000.0,"input","-100,000,000.0 - 10,000,000.0");
         program::initial_down_track_position = m;
         return true;
      }
      //--------------------------------------------------------------------
      test="LFA";
      if(word==test){
         program::LFA = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="LFA-track-field-step";
      if(word==test){
         double m=atof(value.c_str());
         vin::check_for_valid_value(m, word, line, prefix, unit, "none", 0.0, 36000.0,"input","0.0 - 36,000.0");
         program::LFA_scan_field_step = m;
         return true;
      }
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      // add prefix string
      std::string prefix="material:";

      // Check for material id > current array size and if so dynamically expand mp array
      if((unsigned int) super_index + 1 > internal::mp.size() && super_index + 1 < 101) internal::mp.resize(super_index + 1);

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of program namespace
