// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

// Vampire headers
#include "errors.hpp"
#include "stats.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
#include "constants.hpp"
#include "atoms.hpp"
#include "errors.hpp"
#include "sim.hpp"

namespace stats{

//------------------------------------------------------------------------------------------------------
// Function to determine if class is properly initialized
//------------------------------------------------------------------------------------------------------
bool spin_length_statistic_t::is_initialized(){
   return initialized;
}

//------------------------------------------------------------------------------------------------------
// Function to initialize mask
//------------------------------------------------------------------------------------------------------
void spin_length_statistic_t::set_mask(const int in_mask_size, std::vector<int> in_mask){

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
   spin_length.resize(in_mask_size, 0.0);
   mean_spin_length.resize(in_mask_size, 0.0);
   normalisation.resize(in_mask_size, 0.0);

   // calculate number of spins in each mask
   for(int atom=0; atom<num_atoms; ++atom){
      const int mask_id = mask[atom]; // get mask id
      normalisation[mask_id] += 1.0;
   }

   // Calculate normalisation for all CPUs
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &normalisation[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
         zero_list.push_back(mask_id);
      }
   }

   // Set flag indicating correct initialization
   initialized=true;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to get mask needed for gpu acceleration of statistics calculation
//------------------------------------------------------------------------------------------------------
void spin_length_statistic_t::get_mask(std::vector<int>& out_mask){

   // copy data to objects
   out_mask = mask;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to calculate spin length of the system and retain the mean value
//-------------------------------------------------------------------------------------------------------
void spin_length_statistic_t::calculate_spin_length(const std::vector<double>& sx, // full spin vector
                                                    const std::vector<double>& sy,
                                                    const std::vector<double>& sz){

   // initialise spin length to zero [.end() seems to be optimised away by the compiler...]
   std::fill(spin_length.begin(),spin_length.end(),0.0);

   // calculate contributions of spins to spin length
   for(int atom=0; atom<num_atoms; ++atom){
      const int mask_id = mask[atom]; // get mask id
      spin_length[mask_id] += sqrt(sx[atom]*sx[atom] + sy[atom]*sy[atom] + sz[atom]*sz[atom]);
   }

   // Reduce on all CPUS
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &spin_length[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

   // Zero empty mask id's
   for(unsigned int id=0; id<zero_list.size(); ++id) spin_length[zero_list[id]]=0.0;

   // Add spin length to mean
   const int slsize = spin_length.size();
   for(int idx=0; idx<slsize; ++idx){ 
      mean_spin_length[idx] += spin_length[idx];
   }
   mean_counter+=1.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to reset statistical averages
//------------------------------------------------------------------------------------------------------
void spin_length_statistic_t::reset_averages(){

   // reinitialise mean magnetization to zero
   std::fill(mean_spin_length.begin(),mean_spin_length.end(),0.0);

   // reset data counter
   mean_counter = 0.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to output mean spin length value(s) as string
//------------------------------------------------------------------------------------------------------
std::string spin_length_statistic_t::output_mean_spin_length(bool header){

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

   // loop over all spin length value(s)
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header){
         result << name + std::to_string(mask_id) + "|S|";
      }
      else{
         result << mean_spin_length[mask_id]*ic/normalisation[mask_id];
      }
   }

   return result.str();

}

} // end of namespace stats
