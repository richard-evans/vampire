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
#include <string>
#include <vector>

// vampire headers
#include "micromagnetic.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// vio module headers
#include "internal.hpp"

namespace vout{

   //---------------------------------------------------------------------------
   // Function to process output: and screen: list keywords
   //---------------------------------------------------------------------------
   bool match_vout_list(std::string const word, std::string const value, int const line, std::vector<output_t> & output_list){

        std::string prefix="output:";

        std::string test="time-steps";
        if(word==test){
            output_list.push_back(output_t::time_steps);
            return true;
        }
        //--------------------------------------------------------------------
        test="real-time";
        if(word==test){
            output_list.push_back(output_t::real_time);
            return true;
        }
        //--------------------------------------------------------------------
        test="temperature";
        if(word==test){
            output_list.push_back(output_t::temperature);
            return true;
        }
        //--------------------------------------------------------------------
        test="applied-field-strength";
        if(word==test){
            output_list.push_back(output_t::applied_field_strength);
            return true;
        }
        //--------------------------------------------------------------------
        test="applied-field-unit-vector";
        if(word==test){
            output_list.push_back(output_t::applied_field_unit_vector);
            return true;
        }
        //--------------------------------------------------------------------
        test="applied-field-alignment";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(output_t::applied_field_alignment);
            return true;
        }
        //--------------------------------------------------------------------
        test="magnetisation";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(output_t::magnetisation);
            return true;
        }
        //-------------------------------------------------------------------
        test="magnetisation-length";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(output_t::magnetisation_length);
            return true;
        }
        //--------------------------------------------------------------------
        test="mean-magnetisation-length";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(output_t::mean_magnetisation_length);
            return true;
        }
        //--------------------------------------------------------------------
        test="mean-magnetisation";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(output_t::mean_magnetisation);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-magnetisation";
        if(word==test){
            stats::calculate_material_magnetization=true;
            output_list.push_back(output_t::material_magnetisation);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-magnetisation-length";
        if(word==test){
            stats::calculate_material_magnetization=true;
            output_list.push_back(output_t::material_mean_magnetisation_length);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-magnetisation";
        if(word==test){
            stats::calculate_material_magnetization=true;
            output_list.push_back(output_t::material_mean_magnetisation);
            return true;
        }
        //--------------------------------------------------------------------
        test="total-torque";
        if(word==test){
            stats::calculate_system_torque = true;
            output_list.push_back(output_t::total_torque);
            return true;
        }
        //--------------------------------------------------------------------
        test="mean-total-torque";
        if(word==test){
            stats::calculate_system_torque = true;
            output_list.push_back(output_t::mean_total_torque);
            return true;
        }
        //--------------------------------------------------------------------
        test="constraint-phi";
        if(word==test){
            //stats::calculate_torque=true;
            output_list.push_back(output_t::constraint_phi);
            return true;
        }
        //--------------------------------------------------------------------
        test="constraint-theta";
        if(word==test){
            //stats::calculate_torque=true;
            output_list.push_back(output_t::constraint_theta);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-constraint-phi";
        if(word==test){
            //stats::calculate_torque=true;
            output_list.push_back(output_t::material_constraint_phi);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-constraint-theta";
        if(word==test){
            //stats::calculate_torque=true;
            output_list.push_back(output_t::material_constraint_theta);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-torque";
        if(word==test){
            stats::calculate_material_torque = true;
            output_list.push_back(output_t::material_torque);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-torque";
        if(word==test){
            stats::calculate_material_torque = true;
            output_list.push_back(output_t::material_mean_torque);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-standard-deviation";
        if(word==test){
            // Set flags for calculations of standard deviation and magnetization
            stats::calculate_material_standard_deviation=true;
            stats::calculate_material_magnetization=true;
            output_list.push_back(output_t::material_standard_deviation); //AJN
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-susceptibility";
        if(word==test){
            // Set flags for calculations of susceptibility and magnetization
            stats::calculate_system_susceptibility=true;
            stats::calculate_system_magnetization=true;
            output_list.push_back(output_t::mean_susceptibility);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-susceptibility";
        if(word==test){
            // Set flags for calculations of susceptibility and magnetization
            stats::calculate_material_susceptibility=true;
            stats::calculate_material_magnetization=true;
            output_list.push_back(output_t::material_mean_susceptibility);
            return true;
        }
        //-------------------------------------------------------------------
        test="binder-cumulant";
        if(word==test){
            // Set flags for calculations of binder cumulant and magnetization
            stats::calculate_system_binder_cumulant=true;
            stats::calculate_system_magnetization=true;
            output_list.push_back(output_t::binder_cumulant);
            return true;
        }
        //-------------------------------------------------------------------
        test="material-binder-cumulant";
        if(word==test){
            // Set flags for calculations of binder cumulant and magnetization
            stats::calculate_material_binder_cumulant=true;
            stats::calculate_material_magnetization=true;
            output_list.push_back(output_t::material_binder_cumulant);
            return true;
        }
        //-------------------------------------------------------------------
        test="electron-temperature"; // identical to temperature
        if(word==test){
            output_list.push_back(output_t::temperature);
            return true;
        }
        //-------------------------------------------------------------------
        test="phonon-temperature";
        if(word==test){
            output_list.push_back(output_t::phonon_temperature);
            return true;
        }
        //-------------------------------------------------------------------
        test="material-temperature";
        if(word==test){
            output_list.push_back(output_t::material_temperature);
            return true;
        }
        //-------------------------------------------------------------------
        test="material-applied-field-strength";
        if(word==test){
            output_list.push_back(output_t::material_applied_field_strength);
            return true;
        }
        //-------------------------------------------------------------------
        test="material-fmr-field-strength";
        if(word==test){
            output_list.push_back(output_t::material_fmr_field_strength);
            return true;
        }
        //-------------------------------------------------------------------
        test="material-applied-field-alignment";
        if(word==test){
            stats::calculate_material_magnetization=true;
            output_list.push_back(output_t::material_applied_field_alignment);
            return true;
        }
        //-------------------------------------------------------------------
        test="total-energy";
        if(word==test){
            output_list.push_back(output_t::total_energy);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-total-energy";
        if(word==test){
            output_list.push_back(output_t::mean_total_energy);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="anisotropy-energy";
        if(word==test){
            output_list.push_back(output_t::anisotropy_energy);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-anisotropy-energy";
        if(word==test){
            output_list.push_back(output_t::mean_anisotropy_energy);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="exchange-energy";
        if(word==test){
            output_list.push_back(output_t::exchange_energy);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-exchange-energy";
        if(word==test){
            output_list.push_back(output_t::mean_exchange_energy);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="applied-field-energy";
        if(word==test){
            output_list.push_back(output_t::applied_field_energy);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-applied-field-energy";
        if(word==test){
            output_list.push_back(output_t::mean_applied_field_energy);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="magnetostatic-energy";
        if(word==test){
            output_list.push_back(output_t::magnetostatic_energy);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-magnetostatic-energy";
        if(word==test){
            output_list.push_back(output_t::mean_magnetostatic_energy);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="dmi-energy";
        if(word==test){
            output_list.push_back(output_t::dmi_energy);
            stats::calculate_system_energy = true;
            stats::calculate_dmi_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-dmi-energy";
        if(word==test){
            output_list.push_back(output_t::mean_dmi_energy);
            stats::calculate_system_energy = true;
            stats::calculate_dmi_energy = true;
            return true;
        }
        //--------------------------------------------------------------------
        test="height-magnetisation-normalised";
        if(word==test){
            stats::calculate_height_magnetization=true;
            output_list.push_back(output_t::height_magnetisation_normalised);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-height-magnetisation-normalised";
        if(word==test){
            stats::calculate_material_height_magnetization=true;
            output_list.push_back(output_t::material_height_magnetisation_normalised);
            return true;
        }
        //--------------------------------------------------------------------
        test="height-magnetisation";
        if(word==test){
            stats::calculate_height_magnetization=true;
            output_list.push_back(output_t::height_magnetisation);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-height-magnetisation";
        if(word==test){
            stats::calculate_material_height_magnetization=true;
            output_list.push_back(output_t::material_height_magnetisation);
            return true;
        }
        //--------------------------------------------------------------------
        test="fmr-field-strength";
        if(word==test){
            output_list.push_back(output_t::fmr_field_strength);
            return true;
        }
        //--------------------------------------------------------------------
        test="mean-height-magnetisation-length";
        if(word==test){
            stats::calculate_height_magnetization=true;
            output_list.push_back(output_t::mean_height_magnetisation_length);
            return true;
        }
        //--------------------------------------------------------------------
        test="mean-height-magnetisation";
        if(word==test){
            stats::calculate_height_magnetization=true;
            output_list.push_back(output_t::mean_height_magnetisation);
            return true;
        }
        //-------------------------------------------------------------------
        test="mpi-timings";
        if(word==test){
            vmpi::DetailedMPITiming=true;
            output_list.push_back(output_t::mpi_timings);
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-specific-heat";
        if(word==test){
           stats::calculate_system_energy = true;
           stats::calculate_system_specific_heat = true;
           output_list.push_back(output_t::mean_specific_heat);
           return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-specific-heat";
        if(word==test){
           stats::calculate_material_energy = true;
           stats::calculate_material_specific_heat = true;
           output_list.push_back(output_t::material_mean_specific_heat);
           return true;
        }
        //--------------------------------------------------------------------
        test="material-total-energy";
        if(word==test){
           stats::calculate_material_energy = true;
           output_list.push_back(output_t::material_total_energy);
           return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-total-energy";
        if(word==test){
           stats::calculate_material_energy = true;
           output_list.push_back(output_t::material_mean_total_energy);
           return true;
        }
        //--------------------------------------------------------------------
        test="material-dmi-energy";
        if(word==test){
           stats::calculate_material_energy = true;
           stats::calculate_dmi_energy = true;
           output_list.push_back(output_t::material_dmi_energy);
           return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-dmi-energy";
        if(word==test){
           stats::calculate_material_energy = true;
           stats::calculate_dmi_energy = true;
           output_list.push_back(output_t::material_mean_dmi_energy);
           return true;
        }
        //--------------------------------------------------------------------
        test="resistance";
        if(word==test){
           output_list.push_back(output_t::resistance);
           return true;
        }
        //--------------------------------------------------------------------
        test="current";
        if(word==test){
           output_list.push_back(output_t::current);
           return true;
        }
        //--------------------------------------------------------------------
        test="domain-wall-centre";
        if(word==test){
            output_list.push_back(output_t::domain_wall_centre);
            return true;
        }
        //--------------------------------------------------------------------
        test="magneto-resistance";
        if(word==test){
           micromagnetic::enable_resistance = true;
           output_list.push_back(output_t::magneto_resistance);
           return true;
        }
        //--------------------------------------------------------------------
        test="lfa-ms";
        if(word==test){
           output_list.push_back(output_t::lfa_ms);
           return true;
        }
        test="track-pos-x";
        if(word==test){
           output_list.push_back(output_t::track_pos_x);
           return true;
        }
        test="track-pos-z";
        if(word==test){
           output_list.push_back(output_t::track_pos_z);
           return true;
        }
        //--------------------------------------------------------------------
        test="fractional-electric-field-strength";
        if(word==test){
           output_list.push_back(output_t::fractional_electric_field_strength);
           return true;
        }
         //--------------------------------------------------------------------
         test="spin-temperature";
         if(word==test){
            stats::calculate_system_spin_temp = true;
            output_list.push_back(output_t::spin_temperature);
            return true;

         }
         //--------------------------------------------------------------------
         test="mean-spin-temperature";
         if(word==test){
            stats::calculate_system_spin_temp = true;
            output_list.push_back(output_t::mean_spin_temperature);
            return true;
         }
         //--------------------------------------------------------------------
         test="material-spin-temperature";
         if(word==test){
            stats::calculate_material_spin_temp = true;
            output_list.push_back(output_t::material_spin_temperature);
            return true;
         }
         //--------------------------------------------------------------------
         test="material-mean-spin-temperature";
         if(word==test){
            stats::calculate_material_spin_temp = true;
            output_list.push_back(output_t::material_mean_spin_temperature);
            return true;
         }
         //--------------------------------------------------------------------
         test="mean-spin-length";
         if(word==test){
            // Set flags for calculation of spin length
            stats::calculate_system_spin_length=true;
            output_list.push_back(output_t::mean_spin_length);
            return true;
         }
         //--------------------------------------------------------------------
         test="material-mean-spin-length";
         if(word==test){
            // Set flags for calculation of spin length
            stats::calculate_material_spin_length=true;
            output_list.push_back(output_t::material_mean_spin_length);
            return true;
         }
         //--------------------------------------------------------------------
         test="mean-height-spin-length";
         if(word==test){
            // Set flags for calculation of spin length
            stats::calculate_height_spin_length=true;
            output_list.push_back(output_t::mean_height_spin_length);
            return true;
         }
         //--------------------------------------------------------------------
         test="potential-energy";
         if(word==test){
            stats::calculate_system_sld_energy = true;
            output_list.push_back(output_t::potential_energy);
            return true;
         }
         //--------------------------------------------------------------------
         test="kinetic-energy";
         if(word==test){
            stats::calculate_system_sld_energy = true;
            output_list.push_back(output_t::kinetic_energy);
            return true;
         }
         //--------------------------------------------------------------------
         test="sld-exchange-energy";
         if(word==test){
            stats::calculate_system_sld_energy = true;
            output_list.push_back(output_t::sld_exchange_energy);
            return true;
         }
         //--------------------------------------------------------------------
         test="sld-coupling-energy";
         if(word==test){
            stats::calculate_system_sld_energy = true;
            output_list.push_back(output_t::sld_coupling_energy);
            return true;
         }
         //--------------------------------------------------------------------
         test="sld-total-energy";
         if(word==test){
            stats::calculate_system_sld_energy = true;
            output_list.push_back(output_t::sld_total_energy);
            return true;
         }

         //--------------------------------------------------------------------
        test="mean-coupling-field-vec";
        if(word==test){
        output_list.push_back(output_t::mean_coupling_field_vec);
         return true;
      }

         //--------------------------------------------------------------------
         //test="sld-total-spin-energy";
         //if(word==test){
         //   stats::calculate_system_sld_energy = true;
         //   output_list.push_back(86);
         //   return true;
         //}
         //--------------------------------------------------------------------
         //lattice temperature statistics
         //--------------------------------------------------------------------
         test="lattice-temperature";
         if(word==test){
            stats::calculate_system_lattice_temp = true;
            output_list.push_back(output_t::lattice_temperature);
            return true;
         }
         test="mean-lattice-temperature";
         if(word==test){
            stats::calculate_system_lattice_temp = true;
            output_list.push_back(output_t::mean_lattice_temperature);
            return true;
         }
         //--------------------------------------------------------------------
         test="material-lattice-temperature";
         if(word==test){
            stats::calculate_material_lattice_temp = true;
            output_list.push_back(output_t::material_lattice_temperature);
            return true;
         }
         //--------------------------------------------------------------------
         test="material-mean-lattice-temperature";
         if(word==test){
            stats::calculate_material_lattice_temp = true;
            output_list.push_back(output_t::material_mean_lattice_temperature);
            return true;
         }
        //--------------------------------------------------------------------
        test="gnuplot-array-format";
        if(word==test){
           vout::gnuplot_array_format=true;
           return true;
        }
        //--------------------------------------------------------------------
        test="output-rate";
        if(word==test){
            int r=atoi(value.c_str());
            vin::check_for_valid_int(r, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
            vout::output_rate=r;
            return true;
        }

        //--------------------------------------------------------------------
        // Keyword not found
        //--------------------------------------------------------------------
        return false;

   }

} // end of vout namespace
