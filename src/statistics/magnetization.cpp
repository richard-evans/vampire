//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------
//
// In principle better to use kahan summation for accumulated statistics to
// avoid accumulated errors, but up to 20% performance loss for unusual edge
// case. Code left here for future reference
//
// for(int idx = 0; idx < msize; ++idx){
//    const double sum = mean_magnetization[idx];  // load magnetization to temp
//    const double acc = kahan_roundoff[idx];      // save accumulated error
//    const double mag = magnetization[idx];       // load magnetisation
//    const double err = mag - acc;                // net error
//    const double cor = sum + y;                  // corrected sum
//    kahan_roundoff[idx] = (corr - sum) - err;    // retain roundoff for next iteration
//    mean_magnetization[idx] = corr;              // save final sum with correction
// }
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
// Vampire headers
#include "errors.hpp"
#include "stats.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
#include "lsf_mc.hpp"

// file scope function in annonymous namespace to expand statistic type string
namespace {

std::string expand_str(std::string str){

   if(str == "s") return "system_magnetization";
   if(str == "g") return "grain_magnetization";
   if(str == "m") return "material_magnetization";
   if(str == "mg") return "material_grain_magnetization";
   if(str == "h") return "height_magnetization";
   if(str == "mh") return "material_height_magnetization";
   if(str == "mgh") return "material_grain_height_magnetization";

   return "";

}

}

namespace stats{

//------------------------------------------------------------------------------------------------------
// Constructor to initialize data structures
//------------------------------------------------------------------------------------------------------
//magnetization_statistic_t::magnetization_statistic_t (std::string n): initialized(false){}

//------------------------------------------------------------------------------------------------------
// Function to determine if class is properly initialized
//------------------------------------------------------------------------------------------------------
bool magnetization_statistic_t::is_initialized(){
   return initialized;
}

//------------------------------------------------------------------------------------------------------
// Function to initialize mask
//------------------------------------------------------------------------------------------------------
void magnetization_statistic_t::set_mask(const int in_mask_size, std::vector<int> in_mask, const std::vector<double>& mm){

   // Check that mask values never exceed mask_size
   for(unsigned int atom=0; atom<in_mask.size(); ++atom){
      if(in_mask[atom] > in_mask_size-1){
         terminaltextcolor(RED);
         std::cerr << "Programmer Error - mask id " << in_mask[atom] << " is greater than number of elements for mask "<< in_mask_size << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Programmer Error - mask id " << in_mask[atom] << " is greater than number of elements for mask "<< in_mask_size << std::endl;
         err::vexit();
      }
   }

   // save mask to internal storage
   num_atoms = in_mask.size();
   mask_size = in_mask_size - 1; // last element contains magnetization for non-magnetic atoms
   mean_counter = 0.0;
   mask=in_mask; // copy contents of vector
   magnetization.resize(4 * in_mask_size, 0.0);
   mean_magnetization.resize(4 * in_mask_size, 0.0);
   saturation.resize(in_mask_size, 0.0);

   // calculate contributions of spins to each magetization category
   for(int atom=0; atom<num_atoms; ++atom){
      const int mask_id = mask[atom]; // get mask id
      saturation[mask_id] += mm[atom];
   }

   // Add saturation for all CPUs
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &saturation[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

   // determine mask id's with no atoms
   std::vector<int> num_atoms_in_mask(in_mask_size,0);
   for(unsigned int atom=0; atom<in_mask.size(); ++atom){
      int mask_id = in_mask[atom];
      // add atoms to mask
      num_atoms_in_mask[mask_id]++;
   }

   // Reduce on all CPUs
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &num_atoms_in_mask[0], mask_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif

   // Check for no atoms in mask on any CPU
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      // if no atoms exist then add to zero list
      if(num_atoms_in_mask[mask_id]==0){
         zero_list.push_back(4*mask_id+0);
         zero_list.push_back(4*mask_id+1);
         zero_list.push_back(4*mask_id+2);
         zero_list.push_back(4*mask_id+3);
      }
   }

   // Set flag indicating correct initialization
   initialized=true;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to get mask needed for gpu acceleration of statistics calculation
//------------------------------------------------------------------------------------------------------
void magnetization_statistic_t::get_mask(std::vector<int>& out_mask, std::vector<double>& out_saturation){

   // copy data to objects
   out_mask = mask;
   out_saturation = saturation;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to calculate magnetisation of spins given a mask and place result in a magnetization array
//------------------------------------------------------------------------------------------------------
void magnetization_statistic_t::calculate_magnetization(const std::vector<double>& sx, // spin unit vector
                                                        const std::vector<double>& sy,
                                                        const std::vector<double>& sz,
                                                        const std::vector<double>& mm){

   // initialise magnetization to zero [.end() seems to be optimised away by the compiler...]
   std::fill(magnetization.begin(),magnetization.end(),0.0);

   // calculate contributions of spins to each magetization category
   for(int atom=0; atom<num_atoms; ++atom){

      const int mask_id = mask[atom]; // get mask id
      magnetization[4*mask_id + 0] += sx[atom]*mm[atom];
      magnetization[4*mask_id + 1] += sy[atom]*mm[atom];
      magnetization[4*mask_id + 2] += sz[atom]*mm[atom];
      magnetization[4*mask_id + 3] += mm[atom];
   }

   // Reduce on all CPUS
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &magnetization[0], 4*mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

   // Calculate magnetisation length and normalize
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      double msat = magnetization[4*mask_id + 3];
      double magm = sqrt(magnetization[4*mask_id + 0]*magnetization[4*mask_id + 0] +
                         magnetization[4*mask_id + 1]*magnetization[4*mask_id + 1] +
                         magnetization[4*mask_id + 2]*magnetization[4*mask_id + 2]);

      // normalize to msat  // this is what we want std_dev of in time - AJN
      magnetization[4*mask_id + 0] = magnetization[4*mask_id + 0]/magm; // unit vector // x - AJN
      magnetization[4*mask_id + 1] = magnetization[4*mask_id + 1]/magm;                // y
      magnetization[4*mask_id + 2] = magnetization[4*mask_id + 2]/magm;                // z
      magnetization[4*mask_id + 3] = magm/msat; // m/m_s                               // m
   }

   // Zero empty mask id's
   for(unsigned int id=0; id<zero_list.size(); ++id) magnetization[zero_list[id]]=0.0;

   // Add magnetisation to mean
   const int msize = magnetization.size();
   for(int idx=0; idx<msize; ++idx) mean_magnetization[idx] += magnetization[idx];
   mean_counter+=1.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to get magnetisation data
//------------------------------------------------------------------------------------------------------
const std::vector<double>& magnetization_statistic_t::get_magnetization(){

   return magnetization;

}

//------------------------------------------------------------------------------------------------------
// Function to write mean magnetisation data to a checkpoint file
//------------------------------------------------------------------------------------------------------
void magnetization_statistic_t::save_checkpoint(std::ofstream& chkfile){

   const uint64_t num_elements = mean_magnetization.size();

   chkfile.write(reinterpret_cast<const char*>(&num_elements),sizeof(uint64_t));
   chkfile.write(reinterpret_cast<const char*>(&mean_counter),sizeof(double));
   chkfile.write(reinterpret_cast<const char*>(&mean_magnetization[0]),sizeof(double)*mean_magnetization.size());
   return;

}

//------------------------------------------------------------------------------------------------------
// Function to write mean magnetisation data to a checkpoint file
//------------------------------------------------------------------------------------------------------
void magnetization_statistic_t::load_checkpoint(std::ifstream& chkfile, bool chk_continue){

   // load number of elements to see how much data to read
   uint64_t num_elements = 0;
   chkfile.read((char*)&num_elements,sizeof(uint64_t));

   // set up data storage for reading
   double read_mean_counter = 0.0;
   std::vector<double> read_mean_magnetization(num_elements, 0.0);

   // read data elements
   chkfile.read((char*)&read_mean_counter,sizeof(double));
   chkfile.read((char*)&read_mean_magnetization[0],sizeof(double)*num_elements);

   // check that simulation is a continuation (in the case of not continuing do nothing)
   if(chk_continue){

      // check that the number of elements (materials, heights, etc) is the same
      if(num_elements == mean_magnetization.size()){

         // load mean counter and magnetization into class variables
         mean_counter = read_mean_counter;
         mean_magnetization = read_mean_magnetization;

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
// Function to set magnetisation data
//------------------------------------------------------------------------------------------------------
void magnetization_statistic_t::set_magnetization(std::vector<double>& new_magnetization, std::vector<double>& new_mean_magnetization, long counter){

   // copy magnetization vector
   magnetization = new_magnetization;
   //magnetisation.swap(new_magnetization); fatsre but too dangerous?

   const size_t array_size = mean_magnetization.size();
   for(size_t i=0; i< array_size; ++i){
      mean_magnetization[i] += new_mean_magnetization[i];
   }

   // update counter
   mean_counter+=double(counter);

}

//------------------------------------------------------------------------------------------------------
// Function to reset magnetization averages
//------------------------------------------------------------------------------------------------------
void magnetization_statistic_t::reset_magnetization_averages(){

   // reinitialise mean magnetization to zero
   std::fill(mean_magnetization.begin(),mean_magnetization.end(),0.0);

   // reset data counter
   mean_counter = 0.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised magnetisation values as string
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_normalized_magnetization(bool header){

   // result string stream
   std::ostringstream res;
   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }
   vout::fixed_width_output result(res,vout::fw_size);
   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header){
         result << name + std::to_string(mask_id) + "M_norm_x"
                << name + std::to_string(mask_id) + "M_norm_y"
                << name + std::to_string(mask_id) + "M_norm_z"
                << name + std::to_string(mask_id) + "M_norm_l";
      }
      else{
         result << magnetization[4*mask_id + 0]
                << magnetization[4*mask_id + 1]
                << magnetization[4*mask_id + 2]
                << magnetization[4*mask_id + 3];
      }
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output actual magnetisation values as string (in Bohr magnetons)
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_magnetization(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }
   vout::fixed_width_output result(res,vout::fw_size);

   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
   if(header){
      result << name + std::to_string(mask_id) + "M_x"
             << name + std::to_string(mask_id) + "M_y"
             << name + std::to_string(mask_id) + "M_z"
             << name + std::to_string(mask_id) + "M_l";
   }
   else{
      result << magnetization[4*mask_id + 0]
             << magnetization[4*mask_id + 1]
             << magnetization[4*mask_id + 2]
             << magnetization[4*mask_id + 3]*saturation[mask_id];
   }
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised mean magnetisation length values as string
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_normalized_magnetization_length(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }

   vout::fixed_width_output result(res,vout::fw_size);
   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
       if(header){
          result << name + std::to_string(mask_id) + "M_norm_l";
       }
       else{
          result << magnetization[4*mask_id + 3];
       }

   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised mean magnetisation values as string
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_normalized_mean_magnetization(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }

   vout::fixed_width_output result(res,vout::fw_size);
   // inverse number of data samples
   const double ic = 1.0/mean_counter;

   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header){
         result << name + std::to_string(mask_id) + "M_norm_mean_x"
                << name + std::to_string(mask_id) + "M_norm_mean_y"
                << name + std::to_string(mask_id) + "M_norm_mean_z"
                << name + std::to_string(mask_id) + "M_norm_mean_l";
      }
      else{
         result << mean_magnetization[4*mask_id + 0]*ic
                << mean_magnetization[4*mask_id + 1]*ic
                << mean_magnetization[4*mask_id + 2]*ic
                << mean_magnetization[4*mask_id + 3]*ic;
      }
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised mean magnetisation length values as string
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_normalized_mean_magnetization_length(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }

   vout::fixed_width_output result(res,vout::fw_size);
   // inverse number of data samples
   const double ic = 1.0/mean_counter;

   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
       if(header){
           result << name + std::to_string(mask_id) + "M_norm_mean_l";
       }else{
           result << mean_magnetization[4*mask_id + 3]*ic;
       }
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised mean magnetisation length values as string
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_normalized_magnetization_dot_product(const std::vector<double>& vec,bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }

   vout::fixed_width_output result(res,vout::fw_size);
   // check vector has correct size
   if(vec.size()!=3){
      terminaltextcolor(RED);
      std::cerr << "Programmer Error - dot product vector in magnetization statistics must have three elements." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Programmer Error - dot product vector in magnetization statistics must have three elements." << std::endl;
      err::vexit();
   }

   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header){
         result << name + std::to_string(mask_id) + "M_norm_dot_B";
      }
      else{
         const double mm  = magnetization[4*mask_id + 3];
         const double mhx = magnetization[4*mask_id + 0]*mm*vec[0];
         const double mhy = magnetization[4*mask_id + 1]*mm*vec[1];
         const double mhz = magnetization[4*mask_id + 2]*mm*vec[2];
         result << mhx + mhy + mhz;
      }
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output actual magnetisation values as string (in Bohr magnetons)
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_mean_magnetization_length(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }

   vout::fixed_width_output result(res,vout::fw_size);
   // inverse number of data samples
   const double ic = 1.0/mean_counter;

   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header) result << name + std::to_string(mask_id) + "M_mean_l";
      else result << mean_magnetization[4*mask_id + 3]*saturation[mask_id]*ic;
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised mean magnetisation values as string
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_mean_magnetization(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }

   vout::fixed_width_output result(res,vout::fw_size);
   // inverse number of data samples
   const double ic = 1.0/mean_counter;

   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header){
         result << name + std::to_string(mask_id) + "M_mean_x"
                << name + std::to_string(mask_id) + "M_mean_y"
                << name + std::to_string(mask_id) + "M_mean_z"
                << name + std::to_string(mask_id) + "M_mean_l";
      }
      else{
         result << mean_magnetization[4*mask_id + 0]*ic
                << mean_magnetization[4*mask_id + 1]*ic
                << mean_magnetization[4*mask_id + 2]*ic
                << mean_magnetization[4*mask_id + 3]*saturation[mask_id]*ic;
      }
   }

   return result.str();

}

} // end of namespace stats
