//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2022. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

// Vampire headers
#include "constants.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

namespace stats{

//------------------------------------------------------------------------------------------------------
// Constructor to initialize data structures
//------------------------------------------------------------------------------------------------------
//torque_statistic_t::torque_statistic_t (std::string n): initialized(false){}

//------------------------------------------------------------------------------------------------------
// Function to determine if class is properly initialized
//------------------------------------------------------------------------------------------------------
bool torque_statistic_t::is_initialized(){
   return initialized;
}

//------------------------------------------------------------------------------------------------------
// Function to initialize mask
//------------------------------------------------------------------------------------------------------
void torque_statistic_t::set_mask(const int in_mask_size, std::vector<int> in_mask, const std::vector<double>& mm){

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
   mask_size = in_mask_size - 1; // last element contains torques for non-magnetic atoms
   mean_counter = 0.0;
   mask=in_mask; // copy contents of vector
   torque.resize(3 * in_mask_size, 0.0);
   mean_torque.resize(3 * in_mask_size, 0.0);

   // determine mask id's with no atoms
   num_atoms_in_mask.resize(in_mask_size,0);
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
         zero_list.push_back(3*mask_id+0);
         zero_list.push_back(3*mask_id+1);
         zero_list.push_back(3*mask_id+2);
      }
   }

   // Set flag indicating correct initialization
   initialized=true;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to get mask needed for gpu acceleration of statistics calculation
//------------------------------------------------------------------------------------------------------
void torque_statistic_t::get_mask(std::vector<int>& out_mask){

   // copy data to objects
   out_mask = mask;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to calculate torques on spins given a mask and place result in a torque array
//------------------------------------------------------------------------------------------------------
void torque_statistic_t::calculate_torque(const std::vector<double>& sx, // spin unit vector
                                          const std::vector<double>& sy,
                                          const std::vector<double>& sz,
                                          const std::vector<double>& bxs, // spin fields (tesla)
                                          const std::vector<double>& bys,
                                          const std::vector<double>& bzs,
                                          const std::vector<double>& bxe, // external fields (tesla)
                                          const std::vector<double>& bye,
                                          const std::vector<double>& bze,
                                          const std::vector<double>& mm){

   // initialise torques to zero [.end() seems to be optimised away by the compiler...]
   std::fill(torque.begin(),torque.end(),0.0);

   // check for Monte Carlo solvers and recalculate fields
   if(sim::integrator == sim::monte_carlo || sim::integrator == sim::cmc || sim::integrator == sim::hybrid_cmc){
      const int64_t num_atoms = sx.size();
      sim::calculate_spin_fields(0, num_atoms);
      sim::calculate_external_fields(0, num_atoms);
   }

   // calculate contributions of spins to each magetization category
   for(int atom=0; atom < num_atoms; ++atom){

      const int mask_id = mask[atom]; // get mask id

      // get atomic moment
		const double mu = mm[atom];

		// Store local spin in Sand local field in H
		const double S[3] = {sx[atom]*mu,         sy[atom]*mu,         sz[atom]*mu        };
		const double B[3] = {bxs[atom]+bxe[atom], bys[atom]+bye[atom], bzs[atom]+bze[atom]};

		torque[3*mask_id + 0] += S[1]*B[2]-S[2]*B[1];
		torque[3*mask_id + 1] += S[2]*B[0]-S[0]*B[2];
		torque[3*mask_id + 2] += S[0]*B[1]-S[1]*B[0];

	}

   // Reduce on all CPUS
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &torque[0], 3*mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

   // Calculate magnetisation length and normalize
   /*for(int mask_id=0; mask_id < mask_size; ++mask_id){

      // determine inverse number of atoms in mask
      double inv_atoms_in_mask = 1.0 / double(num_atoms_in_mask[mask_id]);

      torque[3*mask_id + 0] *= inv_atoms_in_mask;
      torque[3*mask_id + 1] *= inv_atoms_in_mask;
      torque[3*mask_id + 2] *= inv_atoms_in_mask;

   }*/

   // Zero empty mask id's
   for(unsigned int id=0; id<zero_list.size(); ++id) torque[zero_list[id]]=0.0;

   // Add torque to mean
   const int tsize = torque.size();
   for(int idx = 0; idx < tsize; ++idx) mean_torque[idx] += torque[idx];
   mean_counter+=1.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to get torque data
//------------------------------------------------------------------------------------------------------
const std::vector<double>& torque_statistic_t::get_torque(){

   return torque;

}

//------------------------------------------------------------------------------------------------------
// Function to set torque data (for checkpoints)
//------------------------------------------------------------------------------------------------------
void torque_statistic_t::set_torque(std::vector<double>& new_torque, std::vector<double>& new_mean_torque, long counter){

   // copy torque vector
   torque = new_torque;

   const size_t array_size = mean_torque.size();
   for(size_t i=0; i< array_size; ++i){
      mean_torque[i] += new_mean_torque[i];
   }

   // update counter
   mean_counter += double(counter);

}

//------------------------------------------------------------------------------------------------------
// Function to reset torque averages
//------------------------------------------------------------------------------------------------------
void torque_statistic_t::reset_torque_averages(){

   // reinitialise mean torque to zero
   std::fill(mean_torque.begin(),mean_torque.end(),0.0);

   // reset data counter
   mean_counter = 0.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to output actual torque values as string (in Joules)
//------------------------------------------------------------------------------------------------------
std::string torque_statistic_t::output_torque(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }
   vout::fixed_width_output result(res,vout::fw_size);

   // loop over all torque values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header){
         result << name + std::to_string(mask_id) + "_tx"
                << name + std::to_string(mask_id) + "_ty"
                << name + std::to_string(mask_id) + "_tz";
      }
      else{
         result << constants::muB * torque[3*mask_id + 0]
                << constants::muB * torque[3*mask_id + 1]
                << constants::muB * torque[3*mask_id + 2];
      }
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output mean torque values as string (in Joules)
//------------------------------------------------------------------------------------------------------
std::string torque_statistic_t::output_mean_torque(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }
   vout::fixed_width_output result(res,vout::fw_size);

   // inverse number of data samples * muB
   const double ic = constants::muB / mean_counter;

   // loop over all torque values
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header){
         result << name + std::to_string(mask_id) + "_mean_tx"
                << name + std::to_string(mask_id) + "_mean_ty"
                << name + std::to_string(mask_id) + "_mean_tz";
      }
      else{
         result << mean_torque[3*mask_id + 0]*ic
                << mean_torque[3*mask_id + 1]*ic
                << mean_torque[3*mask_id + 2]*ic;
      }
   }

   return result.str();

}

} // end of namespace stats
