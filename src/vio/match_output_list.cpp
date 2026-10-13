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
   bool match_vout_list(std::string const word, std::string const value, int const line, std::vector<unsigned int> & output_list){

        std::string prefix="output:";

        std::string test="time-steps";
        if(word==test){
            output_list.push_back(0);
            return true;
        }
        //--------------------------------------------------------------------
        test="real-time";
        if(word==test){
            output_list.push_back(1);
            return true;
        }
        //--------------------------------------------------------------------
        test="temperature";
        if(word==test){
            output_list.push_back(2);
            return true;
        }
        //--------------------------------------------------------------------
        test="applied-field-strength";
        if(word==test){
            output_list.push_back(3);
            return true;
        }
        //--------------------------------------------------------------------
        test="applied-field-unit-vector";
        if(word==test){
            output_list.push_back(4);
            return true;
        }
        //--------------------------------------------------------------------
        test="applied-field-alignment";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(12);
            return true;
        }
        //--------------------------------------------------------------------
        test="magnetisation";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(5);
            return true;
        }
        //-------------------------------------------------------------------
        test="magnetisation-length";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(6);
            return true;
        }
        //--------------------------------------------------------------------
        test="mean-magnetisation-length";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(7);
            return true;
        }
        //--------------------------------------------------------------------
        test="mean-magnetisation";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(48);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-magnetisation";
        if(word==test){
            stats::calculate_material_magnetization=true;
            output_list.push_back(8);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-magnetisation-length";
        if(word==test){
            stats::calculate_material_magnetization=true;
            output_list.push_back(9);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-magnetisation";
        if(word==test){
            stats::calculate_material_magnetization=true;
            output_list.push_back(49);
            return true;
        }
        //--------------------------------------------------------------------
        test="total-torque";
        if(word==test){
            stats::calculate_system_torque = true;
            output_list.push_back(14);
            return true;
        }
        //--------------------------------------------------------------------
        test="mean-total-torque";
        if(word==test){
            stats::calculate_system_torque = true;
            output_list.push_back(15);
            return true;
        }
        //--------------------------------------------------------------------
        test="constraint-phi";
        if(word==test){
            //stats::calculate_torque=true;
            output_list.push_back(16);
            return true;
        }
        //--------------------------------------------------------------------
        test="constraint-theta";
        if(word==test){
            //stats::calculate_torque=true;
            output_list.push_back(17);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-constraint-phi";
        if(word==test){
            //stats::calculate_torque=true;
            output_list.push_back(18);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-constraint-theta";
        if(word==test){
            //stats::calculate_torque=true;
            output_list.push_back(19);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-torque";
        if(word==test){
            stats::calculate_material_torque = true;
            output_list.push_back(31);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-torque";
        if(word==test){
            stats::calculate_material_torque = true;
            output_list.push_back(20);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-standard-deviation";
        if(word==test){
            // Set flags for calculations of standard deviation and magnetization
            stats::calculate_material_standard_deviation=true;
            stats::calculate_material_magnetization=true;
            output_list.push_back(999); //AJN
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-susceptibility";
        if(word==test){
            // Set flags for calculations of susceptibility and magnetization
            stats::calculate_system_susceptibility=true;
            stats::calculate_system_magnetization=true;
            output_list.push_back(21);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-susceptibility";
        if(word==test){
            // Set flags for calculations of susceptibility and magnetization
            stats::calculate_material_susceptibility=true;
            stats::calculate_material_magnetization=true;
            output_list.push_back(50);
            return true;
        }
        //-------------------------------------------------------------------
        test="binder-cumulant";
        if(word==test){
            // Set flags for calculations of binder cumulant and magnetization
            stats::calculate_system_binder_cumulant=true;
            stats::calculate_system_magnetization=true;
            output_list.push_back(998);
            return true;
        }
        //-------------------------------------------------------------------
        test="material-binder-cumulant";
        if(word==test){
            // Set flags for calculations of binder cumulant and magnetization
            stats::calculate_material_binder_cumulant=true;
            stats::calculate_material_magnetization=true;
            output_list.push_back(997);
            return true;
        }
        //-------------------------------------------------------------------
        test="electron-temperature"; // identical to temperature
        if(word==test){
            output_list.push_back(2);
            return true;
        }
        //-------------------------------------------------------------------
        test="phonon-temperature";
        if(word==test){
            output_list.push_back(22);
            return true;
        }
        //-------------------------------------------------------------------
        test="material-temperature";
        if(word==test){
            output_list.push_back(23);
            return true;
        }
        //-------------------------------------------------------------------
        test="material-applied-field-strength";
        if(word==test){
            output_list.push_back(24);
            return true;
        }
        //-------------------------------------------------------------------
        test="material-fmr-field-strength";
        if(word==test){
            output_list.push_back(25);
            return true;
        }
        //-------------------------------------------------------------------
        test="material-applied-field-alignment";
        if(word==test){
            stats::calculate_material_magnetization=true;
            output_list.push_back(26);
            return true;
        }
        //-------------------------------------------------------------------
        test="total-energy";
        if(word==test){
            output_list.push_back(27);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-total-energy";
        if(word==test){
            output_list.push_back(28);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="anisotropy-energy";
        if(word==test){
            output_list.push_back(29);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-anisotropy-energy";
        if(word==test){
            output_list.push_back(30);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="exchange-energy";
        if(word==test){
            output_list.push_back(35);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-exchange-energy";
        if(word==test){
            output_list.push_back(36);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="applied-field-energy";
        if(word==test){
            output_list.push_back(37);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-applied-field-energy";
        if(word==test){
            output_list.push_back(38);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="magnetostatic-energy";
        if(word==test){
            output_list.push_back(39);
            stats::calculate_system_energy = true;
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-magnetostatic-energy";
        if(word==test){
            output_list.push_back(40);
            stats::calculate_system_energy = true;
            return true;
        }
        //--------------------------------------------------------------------
        test="height-magnetisation-normalised";
        if(word==test){
            stats::calculate_height_magnetization=true;
            output_list.push_back(43);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-height-magnetisation-normalised";
        if(word==test){
            stats::calculate_material_height_magnetization=true;
            output_list.push_back(44);
            return true;
        }
        //--------------------------------------------------------------------
        test="height-magnetisation";
        if(word==test){
            stats::calculate_height_magnetization=true;
            output_list.push_back(45);
            return true;
        }
        //--------------------------------------------------------------------
        test="material-height-magnetisation";
        if(word==test){
            stats::calculate_material_height_magnetization=true;
            output_list.push_back(46);
            return true;
        }
        //--------------------------------------------------------------------
        test="fmr-field-strength";
        if(word==test){
            output_list.push_back(47);
            return true;
        }
        //--------------------------------------------------------------------
        test="mean-height-magnetisation-length";
        if(word==test){
            stats::calculate_height_magnetization=true;
            output_list.push_back(51);
            return true;
        }
        //--------------------------------------------------------------------
        test="mean-height-magnetisation";
        if(word==test){
            stats::calculate_height_magnetization=true;
            output_list.push_back(52);
            return true;
        }
        //-------------------------------------------------------------------
        test="mpi-timings";
        if(word==test){
            vmpi::DetailedMPITiming=true;
            output_list.push_back(60);
            return true;
        }
        //-------------------------------------------------------------------
        test="mean-specific-heat";
        if(word==test){
           stats::calculate_system_energy = true;
           stats::calculate_system_specific_heat = true;
           output_list.push_back(61);
           return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-specific-heat";
        if(word==test){
           stats::calculate_material_energy = true;
           stats::calculate_material_specific_heat = true;
           output_list.push_back(62);
           return true;
        }
        //--------------------------------------------------------------------
        test="material-total-energy";
        if(word==test){
           stats::calculate_material_energy = true;
           output_list.push_back(63);
           return true;
        }
        //--------------------------------------------------------------------
        test="material-mean-total-energy";
        if(word==test){
           stats::calculate_material_energy = true;
           output_list.push_back(64);
           return true;
        }
        //--------------------------------------------------------------------
        test="resistance";
        if(word==test){
           output_list.push_back(65);
           return true;
        }
        //--------------------------------------------------------------------
        test="current";
        if(word==test){
           output_list.push_back(66);
           return true;
        }
        //--------------------------------------------------------------------
        test="domain-wall-centre";
        if(word==test){
            output_list.push_back(67);
            return true;
        }
        //--------------------------------------------------------------------
        test="magneto-resistance";
        if(word==test){
           micromagnetic::enable_resistance = true;
           output_list.push_back(68);
           return true;
        }
        //--------------------------------------------------------------------
        test="lfa-ms";
        if(word==test){
           output_list.push_back(69);
           return true;
        }
        test="track-pos-x";
        if(word==test){
           output_list.push_back(70);
           return true;
        }
        test="track-pos-z";
        if(word==test){
           output_list.push_back(71);
           return true;
        }
        //--------------------------------------------------------------------
        test="fractional-electric-field-strength";
        if(word==test){
           output_list.push_back(72);
           return true;
        }
         //--------------------------------------------------------------------
         test="spin-temperature";
         if(word==test){
            stats::calculate_system_spin_temp = true;
            output_list.push_back(74);
            return true;

         }
         //--------------------------------------------------------------------
         test="mean-spin-temperature";
         if(word==test){
            stats::calculate_system_spin_temp = true;
            output_list.push_back(75);
            return true;
         }
         //--------------------------------------------------------------------
         test="material-spin-temperature";
         if(word==test){
            stats::calculate_material_spin_temp = true;
            output_list.push_back(77);
            return true;
         }
         //--------------------------------------------------------------------
         test="material-mean-spin-temperature";
         if(word==test){
            stats::calculate_material_spin_temp = true;
            output_list.push_back(76);
            return true;
         }
         //--------------------------------------------------------------------
         test="mean-spin-length";
         if(word==test){
            // Set flags for calculation of spin length
            stats::calculate_system_spin_length=true;
            output_list.push_back(78);
            return true;
         }
         //--------------------------------------------------------------------
         test="material-mean-spin-length";
         if(word==test){
            // Set flags for calculation of spin length
            stats::calculate_material_spin_length=true;
            output_list.push_back(79);
            return true;
         }
         //--------------------------------------------------------------------
         test="mean-height-spin-length";
         if(word==test){
            // Set flags for calculation of spin length
            stats::calculate_height_spin_length=true;
            output_list.push_back(80);
            return true;
         }
         //--------------------------------------------------------------------
         test="potential-energy";
         if(word==test){
            stats::calculate_system_sld_energy = true;
            output_list.push_back(81);
            return true;
         }
         //--------------------------------------------------------------------
         test="kinetic-energy";
         if(word==test){
            stats::calculate_system_sld_energy = true;
            output_list.push_back(82);
            return true;
         }
         //--------------------------------------------------------------------
         test="sld-exchange-energy";
         if(word==test){
            stats::calculate_system_sld_energy = true;
            output_list.push_back(83);
            return true;
         }
         //--------------------------------------------------------------------
         test="sld-coupling-energy";
         if(word==test){
            stats::calculate_system_sld_energy = true;
            output_list.push_back(84);
            return true;
         }
         //--------------------------------------------------------------------
         test="sld-total-energy";
         if(word==test){
            stats::calculate_system_sld_energy = true;
            output_list.push_back(85);
            return true;
         }

         //--------------------------------------------------------------------
        test="mean-coupling-field-vec";
        if(word==test){
        output_list.push_back(91);
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
            output_list.push_back(87);
            return true;
         }
         test="mean-lattice-temperature";
         if(word==test){
            stats::calculate_system_lattice_temp = true;
            output_list.push_back(88);
            return true;
         }
         //--------------------------------------------------------------------
         test="material-lattice-temperature";
         if(word==test){
            stats::calculate_material_lattice_temp = true;
            output_list.push_back(89);
            return true;
         }
         //--------------------------------------------------------------------
         test="material-mean-lattice-temperature";
         if(word==test){
            stats::calculate_material_lattice_temp = true;
            output_list.push_back(90);
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
