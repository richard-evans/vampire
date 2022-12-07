// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
// Vampire headers
#include "errors.hpp"
#include "stats.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
#include "constants.hpp"
#include "atoms.hpp"

namespace stats{

//------------------------------------------------------------------------------------------------------
// Function to initialize data structures
//------------------------------------------------------------------------------------------------------
void spin_length_statistic_t::initialize(stats::magnetization_statistic_t& mag_stat) {

   // Check that magnetization statistic is properly initialized
   if(!mag_stat.is_initialized()){
      terminaltextcolor(RED);
      std::cerr << "Programmer Error - Uninitialized magnetization statistic passed to spin length statistic - please initialize first." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Programmer Error - Uninitialized magnetization statistic passed to spin length statistic - please initialize first." << std::endl;
      err::vexit();
   }

   // Determine number of magnetization statistics
   std::vector<double> temp = mag_stat.get_magnetization();
   num_elements = temp.size()/4;

   // Now set number of spin length values to match
   mean_spin_length.resize(4*num_elements,0.0);

   // initialize mean counter
   mean_counter = 0.0;

   // Set flag indicating correct initialization
   initialized=true;

}

//------------------------------------------------------------------------------------------------------
// Function to calculate spin length of the system and retain the mean value
//-------------------------------------------------------------------------------------------------------
void spin_length_statistic_t::calculate(const std::vector<double>& magnetization){

   // loop over all elements
   for(int id=0; id< num_elements; ++id){

      // copy reduced magnetisation
      const double mx = magnetization[4*id + 0];
      const double my = magnetization[4*id + 1];
      const double mz = magnetization[4*id + 2];
      const double mm = magnetization[4*id + 3];


      mean_spin_length[4*id + 0] += sqrt(mx*mx + my*my + mz*mz);

      const double temp_sl = sqrt(mx*mx + my*my + mz*mz);

      std::cout << num_elements << " " << temp_sl << " " << mx << " " << my << " " << mz << std::endl;

    }

   mean_counter += 1.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to reset statistical averages
//------------------------------------------------------------------------------------------------------
void spin_length_statistic_t::reset_averages(){

   // reinitialise mean magnetization to zero
   std::fill(mean_spin_length.begin(),mean_spin_length.end(),0.0);
   std::fill(mean_spin_length_squared.begin(),mean_spin_length_squared.end(),0.0);

   // reset data counter
   mean_counter = 0.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to output mean specific_heat values as string
//------------------------------------------------------------------------------------------------------
std::string spin_length_statistic_t::output_mean_spin_length(const double temperature,bool header){

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

      const double sl = mean_spin_length[4*id + 0] * imean_counter;

      if(header){
          result<<name<<id<<"_spin_l"<<"\t";
      }else{
          result << sl << "\t";
      }

   }

   return result.str();

}

} // end of namespace stats 