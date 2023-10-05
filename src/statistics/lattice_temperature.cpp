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
bool lattice_temp_statistic_t::is_initialized(){
   return initialized;
}

//------------------------------------------------------------------------------------------------------
// Function to initialize mask
//------------------------------------------------------------------------------------------------------
void lattice_temp_statistic_t::set_mask(const int in_mask_size, std::vector<int> in_mask, const std::vector<double>& mm){

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
   lattice_temp.resize(in_mask_size, 0.0);
   mean_lattice_temp.resize(in_mask_size, 0.0);

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
void lattice_temp_statistic_t::get_mask(std::vector<int>& out_mask){

   // copy data to objects
   out_mask = mask;

   return;

}


void lattice_temp_statistic_t::calculate_lattice_temp(const std::vector<double>& velo_array_x, // coord vectors for atoms
            const std::vector<double>& velo_array_y,
            const std::vector<double>& velo_array_z){

             
           

           std::fill(lattice_temp.begin(),lattice_temp.end(),0.0);
           const int64_t num_atoms = velo_array_x.size();

    double kinetic=0.0;
   // calculate contributions of spins to each magetization category
    for(int atom=0; atom < num_atoms; ++atom){

       const int mask_id = mask[atom]; // get mask id

      // get atomic moment

		double vx = atoms::x_velo_array[atom];
        double vy = atoms::y_velo_array[atom];
        double vz = atoms::z_velo_array[atom];
        kinetic=kinetic+vx*vx+vy*vy+vz*vz;
         //std::cout<<"s"<<S[0]<<"\t"<<S[1]<<"\t"<<S[2]<<"b"<<B[0]<<"\t"<< B[1]<<"\t"<<B[2]<<"mask id:"<<mask_id<<"\t"<<SxH2<<"\t"<<SH<<std::endl;
        lattice_temp[mask_id]=kinetic*atoms::mass_spin_array[atom];
	    //if(atom==2) std::cout<<"STATS at 2  "<<( 2.0 /(3.0*constants::kB_eV))*kinetic*atoms::mass_spin_array[atom]*0.5<<std::endl;

	}
	//std::cout<<"lattice temp from SLD"<<sld::compute_lattice_temperature(0,atoms::num_atoms,atoms::type_array, atoms::x_velo_array,atoms::y_velo_array,atoms::z_velo_array)<<std::endl;
   int  total_atoms=atoms::num_atoms;
   std::cout<<"tot_at "<<total_atoms<<std::endl;


   // Reduce on all CPUS
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &lattice_temp[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     MPI_Allreduce(MPI_IN_PLACE, &total_atoms,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   #endif


   // Zero empty mask id's
   for(unsigned int id=0; id<zero_list.size(); ++id) lattice_temp[zero_list[id]]=0.0;

   const int tsize = lattice_temp.size();
   for(int idx = 0; idx < tsize; ++idx) mean_lattice_temp[idx] += lattice_temp[idx];
   mean_counter+=1.0;

   return;

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
const std::vector<double>& lattice_temp_statistic_t::get_lattice_temp(){

   return lattice_temp;

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void lattice_temp_statistic_t::set_lattice_temp(std::vector<double>& new_lattice_temp, std::vector<double>& new_mean_lattice_temp, long counter){

   lattice_temp = new_lattice_temp;

   const size_t array_size = mean_lattice_temp.size();
   for(size_t i=0; i< array_size; ++i){
      mean_lattice_temp[i] += new_mean_lattice_temp[i];
   }

   // update counter
   mean_counter += double(counter);

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void lattice_temp_statistic_t::reset_lattice_temp_averages(){

   std::fill(mean_lattice_temp.begin(),mean_lattice_temp.end(),0.0);

   // reset data counter
   mean_counter = 0.0;

   return;

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
std::string lattice_temp_statistic_t::output_lattice_temp(bool header){

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
         result << name + std::to_string(mask_id) + "_TL";
      }
      else{
         result << 0.5*( 2.0 /(3.0*constants::kB_eV)) * lattice_temp[mask_id ];
      }
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
std::string lattice_temp_statistic_t::output_mean_lattice_temp(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }
   vout::fixed_width_output result(res,vout::fw_size);

   // inverse number of data samples * muB
   const double ic = 0.5*( 2.0 /(3.0*constants::kB_eV)) / mean_counter;

   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header){
         result << name + std::to_string(mask_id) + "_mean_Ts";
      }
      else{
         result << mean_lattice_temp[mask_id ]*ic;
      }
   }

   return result.str();

}

} // end of namespace stats
