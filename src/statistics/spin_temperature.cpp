//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Mara Strungaru 2022. All rights reserved.
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
#include "sld.hpp"
#include "atoms.hpp"


namespace stats{


//------------------------------------------------------------------------------------------------------
// Function to determine if class is properly initialized
//------------------------------------------------------------------------------------------------------
bool spin_temp_statistic_t::is_initialized(){
   return initialized;
}

//------------------------------------------------------------------------------------------------------
// Function to initialize mask
//------------------------------------------------------------------------------------------------------
void spin_temp_statistic_t::set_mask(const int in_mask_size, std::vector<int> in_mask, const std::vector<double>& mm){

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
   mask_size = in_mask_size - 1; 
   mean_counter = 0.0;
   mask=in_mask; // copy contents of vector
   spin_temp.resize(in_mask_size, 0.0);
   mean_spin_temp.resize(in_mask_size, 0.0);

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
void spin_temp_statistic_t::get_mask(std::vector<int>& out_mask){

   // copy data to objects
   out_mask = mask;

   return;

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void spin_temp_statistic_t::calculate_spin_temp(const std::vector<double>& sx, // spin unit vector
                                          const std::vector<double>& sy,
                                          const std::vector<double>& sz,
                                          const std::vector<double>& bxs, // spin fields (tesla)
                                          const std::vector<double>& bys,
                                          const std::vector<double>& bzs,
                                          const std::vector<double>& bxe, // external fields (tesla)
                                          const std::vector<double>& bye,
                                          const std::vector<double>& bze,
                                          const std::vector<double>& mm){

   std::fill(spin_temp.begin(),spin_temp.end(),0.0);

   
   const int64_t num_tot_atoms = atoms::x_total_spin_field_array.size();
   
   
   std::fill(atoms::x_total_spin_field_array.begin(), atoms::x_total_spin_field_array.end(), 0.0);
   std::fill(atoms::y_total_spin_field_array.begin(), atoms::y_total_spin_field_array.end(), 0.0);
   std::fill(atoms::z_total_spin_field_array.begin(), atoms::z_total_spin_field_array.end(), 0.0);
   sld::compute_fields(0, // first atom for exchange interactions to be calculated
                     num_tot_atoms, // last +1 atom to be calculated
                     atoms::neighbour_list_start_index,
                     atoms::neighbour_list_end_index,
                     atoms::type_array, // type for atom
                     atoms::neighbour_list_array, // list of interactions between atoms
                     atoms::x_coord_array,
                     atoms::y_coord_array,
                     atoms::z_coord_array,
                     atoms::x_spin_array,
                     atoms::y_spin_array,
                     atoms::z_spin_array,
                     atoms::x_total_spin_forces_array,
                     atoms::y_total_spin_forces_array,
                     atoms::z_total_spin_forces_array,
                     atoms::x_total_spin_field_array,
                     atoms::y_total_spin_field_array,
                     atoms::z_total_spin_field_array);

   double SxH2=0.0;
   double SH=0.0;
   // calculate contributions of spins to each magetization category
   for(int atom=0; atom < num_atoms; ++atom){

      const int mask_id = mask[atom]; // get mask id

      // get atomic moment
		const double mu = mm[atom];

		// Store local spin in Sand local field in H
		const double S[3] = {sx[atom],         sy[atom],         sz[atom]        };
		const double B[3] = {bxs[atom], bys[atom], bzs[atom]};

		 double SxHx = S[1]*B[2]-S[2]*B[1];
		 double SxHy = S[2]*B[0]-S[0]*B[2];
		 double SxHz = S[0]*B[1]-S[1]*B[0];
		 SxH2  = SxH2+ SxHx*SxHx + SxHy*SxHy + SxHz*SxHz;
         SH  = SH + S[0]*B[0] + S[1]*B[1] + S[2]*B[2];
         //std::cout<<"s"<<S[0]<<"\t"<<S[1]<<"\t"<<S[2]<<"b"<<B[0]<<"\t"<< B[1]<<"\t"<<B[2]<<"mask id:"<<mask_id<<"\t"<<SxH2<<"\t"<<SH<<std::endl;
         spin_temp[mask_id]=0.5* mu /constants::kB * SxH2 / SH;
         

	}



   // Reduce on all CPUS
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &spin_temp[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif


   // Zero empty mask id's
   for(unsigned int id=0; id<zero_list.size(); ++id) spin_temp[zero_list[id]]=0.0;

   const int tsize = spin_temp.size();
   for(int idx = 0; idx < tsize; ++idx) mean_spin_temp[idx] += spin_temp[idx];
   mean_counter+=1.0;

   return;

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
const std::vector<double>& spin_temp_statistic_t::get_spin_temp(){

   return spin_temp;

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void spin_temp_statistic_t::set_spin_temp(std::vector<double>& new_spin_temp, std::vector<double>& new_mean_spin_temp, long counter){

   spin_temp = new_spin_temp;

   const size_t array_size = mean_spin_temp.size();
   for(size_t i=0; i< array_size; ++i){
      mean_spin_temp[i] += new_mean_spin_temp[i];
   }

   // update counter
   mean_counter += double(counter);

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void spin_temp_statistic_t::reset_spin_temp_averages(){

   std::fill(mean_spin_temp.begin(),mean_spin_temp.end(),0.0);

   // reset data counter
   mean_counter = 0.0;

   return;

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
std::string spin_temp_statistic_t::output_spin_temp(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }
   vout::fixed_width_output result(res,vout::fw_size);

   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header){
         result << name + std::to_string(mask_id) + "_Ts";
      }
      else{
         result << constants::muB * spin_temp[mask_id ];
      }
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
std::string spin_temp_statistic_t::output_mean_spin_temp(bool header){

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

   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header){
         result << name + std::to_string(mask_id) + "_mean_Ts";
      }
      else{
         result << mean_spin_temp[mask_id ]*ic;
      }
   }

   return result.str();

}

} // end of namespace stats
