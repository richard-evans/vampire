//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

// Vampire headers
#include "constants.hpp"
#include "errors.hpp"
#include "stats.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

namespace stats{

//------------------------------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------------------------------
//specific_heat_statistic_t::specific_heat_statistic_t (){}

//------------------------------------------------------------------------------------------------------
// Function to initialize data structures
//------------------------------------------------------------------------------------------------------
void specific_heat_statistic_t::initialize(stats::energy_statistic_t& energy_statistic) {

   // Check that magnetization statistic is properly initialized
   if(!energy_statistic.is_initialized()){
      terminaltextcolor(RED);
      std::cerr << "Programmer Error - Uninitialized energy statistic passed to specific_heat statistic - please initialize first." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Programmer Error - Uninitialized energy statistic passed to specific_heat statistic - please initialize first." << std::endl;
      err::vexit();
   }

   // Determine number of magnetization statistics*4
   std::vector<double> temp = energy_statistic.get_total_energy();
   num_elements = temp.size();

   // Now set number of specific_heat values to match
   mean_specific_heat.resize(num_elements,0.0);
   mean_specific_heat_squared.resize(num_elements,0.0);

   // copy normalisation data for each set of atoms (num atoms)
   normalisation = energy_statistic.normalisation;

   // initialize mean counter
   mean_counter = 0.0;

   // Set flag indicating correct initialization
   initialized=true;

}

//------------------------------------------------------------------------------------------------------
// Function to calculate specific heat of the system and retain the mean value
//
//       C_v =    ( <E^2> - <E>^2 )
//               -------------------
//                    k_B T^2
//
//-------------------------------------------------------------------------------------------------------
void specific_heat_statistic_t::calculate(const std::vector<double>& energy){

   // loop over all elements
   for(int id=0; id< num_elements; ++id){

      // copy energy
      const double E = energy[id];

      mean_specific_heat[id] += E;
      mean_specific_heat_squared[id] += E*E;

   }

   mean_counter += 1.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to reset statistical averages
//------------------------------------------------------------------------------------------------------
void specific_heat_statistic_t::reset_averages(){

   // reinitialise mean magnetization to zero
   std::fill(mean_specific_heat.begin(),mean_specific_heat.end(),0.0);
   std::fill(mean_specific_heat_squared.begin(),mean_specific_heat_squared.end(),0.0);

   // reset data counter
   mean_counter = 0.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to output mean specific_heat values as string
//------------------------------------------------------------------------------------------------------
std::string specific_heat_statistic_t::output_mean_specific_heat(const double temperature,bool header){

   // result string stream
   std::ostringstream result;

   // set custom precision if enabled
   if(vout::custom_precision){
      result.precision(vout::precision);
      if(vout::fixed) result.setf( std::ios::fixed, std::ios::floatfield );
   }

   // determine inverse temperature mu_B/(kB T) (flushing to zero for very low temperatures to avoid NaN)
   const double itemp2 = temperature < 1.e-300 ? +0.0 : constants::muB / ( constants::kB * temperature * temperature );

   // determine inverse mean counter and its square
   const double imean_counter = 1.0 / mean_counter;
   const double imean_counter_sq = 1.0 / (mean_counter*mean_counter);

   // loop over all elements
   for(int id=0; id< num_elements - 1; ++id){ // ignore last element as always contains non-magnetic atoms

      const double prefactor = itemp2 / normalisation[id]; // in 1 / Kelvin per spin
      const double sh = prefactor * ( ( mean_specific_heat_squared[id] * imean_counter ) - ( mean_specific_heat[id] * mean_specific_heat[id] * imean_counter_sq) );

      if(header){
          result<<name<<id<<"_spec_heat"<<"\t";
      }else{
          result << sh << "\t";
      }

   }

   return result.str();

}

} // end of namespace stats
