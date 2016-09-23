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
#include <iostream>
#include <sstream>

// Vampire headers
#include "errors.hpp"
#include "stats.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

namespace stats{

//------------------------------------------------------------------------------------------------------
// Constructor to initialize data structures
//------------------------------------------------------------------------------------------------------
magnetization_statistic_t::magnetization_statistic_t (): initialized(false){}

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
   mask_size = in_mask_size;
   mean_counter=0.0;
   mask=in_mask; // copy contents of vector
   magnetization.resize(4*mask_size,0.0);
   mean_magnetization.resize(4*mask_size,0.0);
   saturation.resize(mask_size,0.0);

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
   std::vector<int> num_atoms_in_mask(mask_size,0);
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

      // normalize to msat
      magnetization[4*mask_id + 0] = magnetization[4*mask_id + 0]/magm; // unit vector
      magnetization[4*mask_id + 1] = magnetization[4*mask_id + 1]/magm;
      magnetization[4*mask_id + 2] = magnetization[4*mask_id + 2]/magm;
      magnetization[4*mask_id + 3] = magm/msat; // m/m_s
   }

   // Zero empty mask id's
   for(unsigned int id=0; id<zero_list.size(); ++id) magnetization[zero_list[id]]=0.0;

   // Add magnetisation to mean
   const int msize = magnetization.size();
   for(int idx=0; idx<msize; ++idx) mean_magnetization[idx]+=magnetization[idx];
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
std::string magnetization_statistic_t::output_normalized_magnetization(){

   // result string stream
   std::ostringstream result;

   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      result << magnetization[4*mask_id + 0] << "\t" << magnetization[4*mask_id + 1] << "\t" << magnetization[4*mask_id + 2] << "\t" << magnetization[4*mask_id + 3] << "\t";
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output actual magnetisation values as string (in Bohr magnetons)
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_magnetization(){

   // result string stream
   std::ostringstream result;

   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      result << magnetization[4*mask_id + 0] << "\t" << magnetization[4*mask_id + 1] << "\t" << magnetization[4*mask_id + 2] << "\t" << magnetization[4*mask_id + 3]*saturation[mask_id] << "\t";
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised mean magnetisation length values as string
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_normalized_magnetization_length(){

   // result string stream
   std::ostringstream result;

   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      result << magnetization[4*mask_id + 3] << "\t";
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised mean magnetisation values as string
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_normalized_mean_magnetization(){

   // result string stream
   std::ostringstream result;

   // inverse number of data samples
   const double ic = 1.0/mean_counter;

   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      result << mean_magnetization[4*mask_id + 0]*ic << "\t" << mean_magnetization[4*mask_id + 1]*ic << "\t" << mean_magnetization[4*mask_id + 2]*ic << "\t" << mean_magnetization[4*mask_id + 3]*ic << "\t";
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised mean magnetisation length values as string
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_normalized_mean_magnetization_length(){

   // result string stream
   std::ostringstream result;

   // inverse number of data samples
   const double ic = 1.0/mean_counter;

   // loop over all magnetization values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      result << mean_magnetization[4*mask_id + 3]*ic << "\t";
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised mean magnetisation length values as string
//------------------------------------------------------------------------------------------------------
std::string magnetization_statistic_t::output_normalized_magnetization_dot_product(const std::vector<double>& vec){

   // result string stream
   std::ostringstream result;

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
      const double mm = magnetization[4*mask_id + 3];
      const double mhx = magnetization[4*mask_id + 0]*mm*vec[0];
      const double mhy = magnetization[4*mask_id + 1]*mm*vec[1];
      const double mhz = magnetization[4*mask_id + 2]*mm*vec[2];
      result << mhx + mhy + mhz << "\t";
   }

   return result.str();

}

} // end of namespace stats
