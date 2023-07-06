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

// file scope function in annonymous namespace to expand statistic type string
namespace {

inline std::string expand_str(std::string str){

   if(str == "s") return "system_susceptibility";
   if(str == "g") return "grain_susceptibility";
   if(str == "m") return "material_susceptibility";

   return "";

}

}

namespace stats{

//------------------------------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------------------------------
//susceptibility_statistic_t::susceptibility_statistic_t (){}

//------------------------------------------------------------------------------------------------------
// Function to initialize data structures
//------------------------------------------------------------------------------------------------------
void susceptibility_statistic_t::initialize(stats::magnetization_statistic_t& mag_stat) {

   // Check that magnetization statistic is properly initialized
   if(!mag_stat.is_initialized()){
      terminaltextcolor(RED);
      std::cerr << "Programmer Error - Uninitialized magnetization statistic passed to susceptibility statistic - please initialize first." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Programmer Error - Uninitialized magnetization statistic passed to susceptibility statistic - please initialize first." << std::endl;
      err::vexit();
   }

   // Determine number of magnetization statistics*4
   std::vector<double> temp = mag_stat.get_magnetization();
   num_elements = temp.size()/4;

   // Now set number of susceptibility values to match
   mean_susceptibility.resize(4*num_elements,0.0);
   mean_susceptibility_squared.resize(4*num_elements,0.0);
   mean_absolute_susceptibility.resize(4*num_elements,0.0);
   mean_absolute_susceptibility_squared.resize(4*num_elements,0.0);

   // copy saturation data
   saturation = mag_stat.saturation;

   // initialize mean counter
   mean_counter = 0.0;

   // Set flag indicating correct initialization
   initialized=true;

}

//------------------------------------------------------------------------------------------------------
// Function to calculate susceptibility of the magnetisation and retain the mean value
//
//       chi_l = sum_i mu_i
//               ----------  ( <m_l^2> - <m_l>^2 )
//                 k_B T
//
//       m_l = sum_i mu_i S_i
//             --------------
//               sum_i mu_i
//
//-------------------------------------------------------------------------------------------------------
void susceptibility_statistic_t::calculate(const std::vector<double>& magnetization){

   // loop over all elements
   for(int id=0; id< num_elements; ++id){

      // copy reduced magnetisation
      const double mx = magnetization[4*id + 0];
      const double my = magnetization[4*id + 1];
      const double mz = magnetization[4*id + 2];
      const double mm = magnetization[4*id + 3];

      mean_susceptibility[4*id + 0]+=mx*mm;
      mean_susceptibility[4*id + 1]+=my*mm;
      mean_susceptibility[4*id + 2]+=mz*mm;
      mean_susceptibility[4*id + 3]+=mm;

      mean_susceptibility_squared[4*id + 0]+=mx*mx*mm*mm;
      mean_susceptibility_squared[4*id + 1]+=my*my*mm*mm;
      mean_susceptibility_squared[4*id + 2]+=mz*mz*mm*mm;
      mean_susceptibility_squared[4*id + 3]+=mm*mm;

      mean_absolute_susceptibility[4*id + 0]+=fabs(mx*mm);
      mean_absolute_susceptibility[4*id + 1]+=fabs(my*mm);
      mean_absolute_susceptibility[4*id + 2]+=fabs(mz*mm);
      mean_absolute_susceptibility[4*id + 3]+=mm;

      mean_absolute_susceptibility_squared[4*id + 0]+=fabs(mx*mx*mm*mm);
      mean_absolute_susceptibility_squared[4*id + 1]+=fabs(my*my*mm*mm);
      mean_absolute_susceptibility_squared[4*id + 2]+=fabs(mz*mz*mm*mm);
      mean_absolute_susceptibility_squared[4*id + 3]+=mm*mm;

   }

   mean_counter+=1.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to write mean susceptibility data to a checkpoint file
//------------------------------------------------------------------------------------------------------
void susceptibility_statistic_t::save_checkpoint(std::ofstream& chkfile){

   const uint64_t num_elements = mean_susceptibility.size();

   chkfile.write(reinterpret_cast<const char*>(&num_elements),sizeof(uint64_t));
   chkfile.write(reinterpret_cast<const char*>(&mean_counter),sizeof(double));
   chkfile.write(reinterpret_cast<const char*>(&mean_susceptibility[0]),sizeof(double)*mean_susceptibility.size());
   chkfile.write(reinterpret_cast<const char*>(&mean_susceptibility_squared[0]),sizeof(double)*mean_susceptibility_squared.size());

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to write mean susceptibility data to a checkpoint file
//------------------------------------------------------------------------------------------------------
void susceptibility_statistic_t::load_checkpoint(std::ifstream& chkfile, bool chk_continue){

   // load number of elements to see how much data to read
   uint64_t num_elements = 0;
   chkfile.read((char*)&num_elements,sizeof(uint64_t));

   // set up data storage for reading
   double read_mean_counter = 0.0;
   std::vector<double> read_mean_susceptibility(num_elements, 0.0);
   std::vector<double> read_mean_susceptibility_squared(num_elements, 0.0);

   // read data elements
   chkfile.read((char*)&read_mean_counter,sizeof(double));
   chkfile.read((char*)&read_mean_susceptibility[0],sizeof(double)*num_elements);
   chkfile.read((char*)&read_mean_susceptibility_squared[0],sizeof(double)*num_elements);

   // check that simulation is a continuation (in the case of not continuing do nothing)
   if(chk_continue){

      // check that the number of elements (materials, heights, etc) is the same
      if(num_elements == mean_susceptibility.size()){

         // load mean counter and data into class variables
         mean_counter = read_mean_counter;
         mean_susceptibility = read_mean_susceptibility;
         mean_susceptibility_squared = read_mean_susceptibility_squared;

      }
      // if not, don't load them (allowing changing of stats after checkpoint)
      // but print out warning message to user
      else{
         zlog << zTs() << "Warning - checkpoint loaded for previously unused statistic " << expand_str(name) << std::endl;
      }
   }

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to reset statistical averages
//------------------------------------------------------------------------------------------------------
void susceptibility_statistic_t::reset_averages(){

   // reinitialise mean magnetization to zero
   std::fill(mean_susceptibility.begin(),mean_susceptibility.end(),0.0);
   std::fill(mean_susceptibility_squared.begin(),mean_susceptibility_squared.end(),0.0);
   std::fill(mean_absolute_susceptibility.begin(),mean_absolute_susceptibility.end(),0.0);
   std::fill(mean_absolute_susceptibility_squared.begin(),mean_absolute_susceptibility_squared.end(),0.0);

   // reset data counter
   mean_counter = 0.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to output mean susceptibility values as string
//------------------------------------------------------------------------------------------------------
std::string susceptibility_statistic_t::output_mean_susceptibility(const double temperature, bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }
   vout::fixed_width_output result(res,vout::fw_size);
   if(!header){

      // determine inverse temperature mu_B/(kB T) (flushing to zero for very low temperatures)
      const double itemp = temperature < 1.e-300 ? 0.0 : 9.274e-24/(1.3806503e-23*temperature);

      // determine inverse mean counter and its square
      const double imean_counter = 1.0/mean_counter;
      const double imean_counter_sq = 1.0/(mean_counter*mean_counter);

      // loop over all elements
      for(int id=0; id< num_elements - 1; ++id){ // ignore last element as always contains non-magnetic atoms

         const double prefactor = itemp*saturation[id]; // in mu_B
         const double sus_x = prefactor*(mean_susceptibility_squared[4*id + 0]*imean_counter-mean_susceptibility[4*id + 0]*mean_susceptibility[4*id + 0]*imean_counter_sq);
         const double sus_y = prefactor*(mean_susceptibility_squared[4*id + 1]*imean_counter-mean_susceptibility[4*id + 1]*mean_susceptibility[4*id + 1]*imean_counter_sq);
         const double sus_z = prefactor*(mean_susceptibility_squared[4*id + 2]*imean_counter-mean_susceptibility[4*id + 2]*mean_susceptibility[4*id + 2]*imean_counter_sq);
         const double sus_m = prefactor*(mean_susceptibility_squared[4*id + 3]*imean_counter-mean_susceptibility[4*id + 3]*mean_susceptibility[4*id + 3]*imean_counter_sq);

         result << sus_x << sus_y << sus_z << sus_m;

      }
   }
   else {
      for(int id=0; id< num_elements - 1; ++id){ // ignore last element as always contains non-magnetic atoms
         result << name + std::to_string(id) + "_sus_x"
                << name + std::to_string(id) + "_sus_y"
                << name + std::to_string(id) + "_sus_z"
                << name + std::to_string(id) + "_sus_m";
      }
   }
   return result.str();

}

} // end of namespace stats
