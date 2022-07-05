//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

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

namespace stats{

//------------------------------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------------------------------
//binder_cumulant_statistic_t::binder_cumulant_statistic_t (){}

//------------------------------------------------------------------------------------------------------
// Function to initialize data structures
//------------------------------------------------------------------------------------------------------
void binder_cumulant_statistic_t::initialize(stats::magnetization_statistic_t& mag_stat) {

   // Check that magnetization statistic is properly initialized
   if(!mag_stat.is_initialized()){
      terminaltextcolor(RED);
      std::cerr << "Programmer Error - Uninitialized magnetization statistic passed to binder cumulant statistic - please initialize first." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Programmer Error - Uninitialized magnetization statistic passed to binder cumulant statistic - please initialize first." << std::endl;
      err::vexit();
   }

   // Determine number of magnetization statistics*4
   std::vector<double> temp = mag_stat.get_magnetization();
   num_elements = temp.size()/4;

   // Now set number of binder cumulant values to match
   binder_cumulant_squared.resize(num_elements,0.0);
   binder_cumulant_fourth_power.resize(num_elements,0.0);

   // initialize mean counter
   mean_counter = 0.0;

   // Set flag indicating correct initialization
   initialized=true;

}

//------------------------------------------------------------------------------------------------------
// Function to calculate binder cumulant of the magnetisation and retain the mean value
//
//       V_L =          <m_l^4>
//             1  -  -------------
//                    3*<m_l^2>^2
//
//       m_l = sum_i mu_i S_i
//             --------------
//               sum_i mu_i
//
//-------------------------------------------------------------------------------------------------------
void binder_cumulant_statistic_t::calculate(const std::vector<double>& magnetization){

   // loop over all elements
   for(int id=0; id< num_elements; ++id){

      // copy reduced magnetisation
      const double mm = magnetization[4*id + 3];

      binder_cumulant_squared[id]+=mm*mm;
      binder_cumulant_fourth_power[id]+=mm*mm*mm*mm;

   }

   mean_counter+=1.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to reset statistical averages
//------------------------------------------------------------------------------------------------------
void binder_cumulant_statistic_t::reset_averages(){

   // reinitialise mean magnetization to zero
   std::fill(binder_cumulant_squared.begin(),binder_cumulant_squared.end(),0.0);
   std::fill(binder_cumulant_fourth_power.begin(),binder_cumulant_fourth_power.end(),0.0);

   // reset data counter
   mean_counter = 0.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to output mean binder cumulant values as string
//------------------------------------------------------------------------------------------------------
std::string binder_cumulant_statistic_t::output_binder_cumulant(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }
   vout::fixed_width_output result(res,vout::fw_size); 
   if(!header){

   // determine inverse mean counter and its square
   const double imean_counter = 1.0/mean_counter;
   const double imean_counter_sq = 1.0/(mean_counter*mean_counter);

   // loop over all elements
   for(int id=0; id< num_elements - 1; ++id){ // ignore last element as always contains non-magnetic atoms

      const double bc = 1 - binder_cumulant_fourth_power[id]*imean_counter/(3.0*binder_cumulant_squared[id]*binder_cumulant_squared[id]*imean_counter_sq);

      result << bc;

   }
   }else{
       for(int id=0; id< num_elements - 1; ++id){ // ignore last element as always contains non-magnetic atoms
          result << name + std::to_string(id) + "_bc";
       }
   }
   return result.str();

}

} // end of namespace stats
