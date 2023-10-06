//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2018. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>

// Vampire headers
#include "anisotropy.hpp"
#include "constants.hpp"
#include "dipole.hpp"
#include "errors.hpp"
#include "exchange.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
#include "sld.hpp"
#include "atoms.hpp"

namespace stats{

//------------------------------------------------------------------------------------------------------
// Generalised class to calculate energy statistics for a subset of atoms defined by mask
//
// The class sets up data structures to hold statistic variables in order to calculate
// instantaneous and mean values of the energy. Specific values also exist for the exchange,
// anisotropy, applied and magnetostatic energies.
//------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------
// Constructor to initialize data structures
//------------------------------------------------------------------------------------------------------
//energy_statistic_t::energy_statistic_t (): initialized(false){}

//------------------------------------------------------------------------------------------------------
// Function to determine if class is properly initialized
//------------------------------------------------------------------------------------------------------
bool sld_energy_statistic_t::is_initialized(){
   return initialized;
}

//------------------------------------------------------------------------------------------------------
// Function to initialize mask
//------------------------------------------------------------------------------------------------------
void sld_energy_statistic_t::set_mask(const int in_mask_size, const std::vector<int> in_mask){

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
   mask_size = in_mask_size - 1; // last element contains energy for non-magnetic atoms
   mean_counter = 0.0;
   mask = in_mask; // copy contents of vector

   // resize arrays to correct mask size (one value per mask) and set to zero
   sld_total_energy.resize(in_mask_size, 0.0);
   sld_exchange_energy.resize(in_mask_size, 0.0);
   sld_coupling_energy.resize(in_mask_size, 0.0);
   kinetic_energy.resize(in_mask_size, 0.0);
   potential_energy.resize(in_mask_size, 0.0);

   mean_sld_total_energy.resize(in_mask_size, 0.0);
   mean_sld_exchange_energy.resize(in_mask_size, 0.0);
   mean_sld_coupling_energy.resize(in_mask_size, 0.0);
   mean_kinetic_energy.resize(in_mask_size, 0.0);
   mean_potential_energy.resize(in_mask_size, 0.0);

   normalisation.resize(in_mask_size, 0.0);
   num_atoms_in_mask.resize(in_mask_size,0);


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
   initialized = true;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to get mask needed for gpu acceleration of statistics calculation
//------------------------------------------------------------------------------------------------------
void sld_energy_statistic_t::get_mask(std::vector<int>& out_mask, std::vector<double>& out_normalisation){

   // copy data to objects
   out_mask = mask;
   out_normalisation = normalisation;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to calculate spin energy given a mask and place result in energy array
//------------------------------------------------------------------------------------------------------
void sld_energy_statistic_t::calculate(const std::vector<double>& sx,  // spin unit vector
                                   const std::vector<double>& sy,
                                   const std::vector<double>& sz,
                                   const std::vector<double>& mm,  // magnetic moment (Tesla)
                                   const std::vector<int>& mat, // material id
                                   const double temperature){

   // initialise energies to zero
   std::fill(         sld_total_energy.begin(),         sld_total_energy.end(), 0.0 );
   std::fill(      sld_exchange_energy.begin(),      sld_exchange_energy.end(), 0.0 );
   std::fill(    sld_coupling_energy.begin(),    sld_coupling_energy.end(), 0.0 );
   std::fill( kinetic_energy.begin(), kinetic_energy.end(), 0.0 );
   std::fill( potential_energy.begin(), potential_energy.end(), 0.0 );
   
   //test
   double pot_eng=sld::compute_potential_energy(0,atoms::num_atoms, atoms::type_array);
   double kin_eng=sld::compute_kinetic_energy(0,atoms::num_atoms,atoms::type_array,atoms::x_velo_array,atoms::y_velo_array, atoms::z_velo_array);
   double sld_exch_eng=sld::compute_exchange_energy(0,atoms::num_atoms);
   double sld_coupl_eng=sld::compute_coupling_energy(0,atoms::num_atoms);
    /*           
   std::cout<<"kin eng "<<kin_eng<<std::endl;
   std::cout<<"pot eng "<<pot_eng<<std::endl;
   std::cout<<"exch eng "<<sld_exch_eng<<std::endl;
   std::cout<<"coupl eng "<<sld_coupl_eng<<std::endl;
*/
   //---------------------------------------------------------------------------
   // Calculate exchange energy (in Tesla)
   //---------------------------------------------------------------------------

   // loop over all atoms in mask
   for( int atom = 0; atom < num_atoms; ++atom ){
      const int mask_id = mask[atom]; // get mask id
      sld_exchange_energy[mask_id] += sld::compute_exchange_energy(atom,atom+1);
   }

   

   //---------------------------------------------------------------------------
   // Calculate anisotropy energy (in Tesla)
   //---------------------------------------------------------------------------
   for( int atom = 0; atom < num_atoms; ++atom ){
      const int mask_id = mask[atom]; // get mask id
      sld_coupling_energy[mask_id] += sld::compute_coupling_energy(atom,atom+1);
   }

   //---------------------------------------------------------------------------
   // Calculate applied field energy (in Tesla)
   //---------------------------------------------------------------------------
   for( int atom = 0; atom < num_atoms; ++atom ){
      const int mask_id = mask[atom]; // get mask id
      kinetic_energy[mask_id] += sld::compute_kinetic_energy(atom,atom+1,atoms::type_array,atoms::x_velo_array,atoms::y_velo_array, atoms::z_velo_array);
   }

   //---------------------------------------------------------------------------
   // Calculate potential field energy (in Tesla)
   //---------------------------------------------------------------------------
   for( int atom = 0; atom < num_atoms; ++atom ){
      const int mask_id = mask[atom]; // get mask id
      potential_energy[mask_id] += sld::compute_potential_energy(atom,atom+1,atoms::type_array);
   }

   
   //---------------------------------------------------------------------------
   // Calculate total energy (in Tesla)
   //---------------------------------------------------------------------------
   for( int mask_id = 0; mask_id < mask_size; ++mask_id ){
      sld_total_energy[mask_id] = sld_exchange_energy[mask_id] +
                              sld_coupling_energy[mask_id] +
                              kinetic_energy[mask_id] +
                              potential_energy[mask_id];
   }

 //  int  total_atoms=atoms::num_atoms;

   //---------------------------------------------------------------------------
   // Reduce on all CPUS
   //---------------------------------------------------------------------------
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE,      &sld_exchange_energy[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE,    &sld_coupling_energy[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &kinetic_energy[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &potential_energy[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE,         &sld_total_energy[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      //MPI_Allreduce(MPI_IN_PLACE, &total_atoms,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   #endif
   
   // divide by number of atoms
   for(int mask_id=0; mask_id < mask_size; ++mask_id){

   // determine inverse number of atoms in mask
      //std::cout<<num_atoms_in_mask[mask_id]<<std::endl;
      double inv_atoms_in_mask = 1.0 / double(num_atoms_in_mask[mask_id]);
      sld_exchange_energy[mask_id] *= inv_atoms_in_mask;
      sld_coupling_energy[mask_id] *= inv_atoms_in_mask;
      sld_total_energy[mask_id] *= inv_atoms_in_mask;
      kinetic_energy[mask_id] *= inv_atoms_in_mask;
      potential_energy[mask_id] *= inv_atoms_in_mask;
  
  }
  


   //---------------------------------------------------------------------------
   // Add energies to mean energies
   //---------------------------------------------------------------------------
   for(int mask_id=0; mask_id<mask_size; ++mask_id ){
      mean_sld_total_energy[mask_id]         += sld_total_energy[mask_id];
      mean_sld_exchange_energy[mask_id]      += sld_exchange_energy[mask_id];
      mean_sld_coupling_energy[mask_id]    += sld_coupling_energy[mask_id];
      mean_potential_energy[mask_id] += potential_energy[mask_id];
      mean_kinetic_energy[mask_id] += kinetic_energy[mask_id];
   }

   // increment mean counter
   mean_counter += 1.0;

   //---------------------------------------------------------------------------
   // Zero empty mask id's
   //---------------------------------------------------------------------------
   for( unsigned int id=0; id < zero_list.size(); ++id ){
              sld_total_energy[ zero_list[id] ] = 0.0;
           sld_exchange_energy[ zero_list[id] ] = 0.0;
         sld_coupling_energy[ zero_list[id] ] = 0.0;
      kinetic_energy[ zero_list[id] ] = 0.0;
      potential_energy[ zero_list[id] ] = 0.0;
   }

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to get const reference for total energy data
//------------------------------------------------------------------------------------------------------
const std::vector<double>& sld_energy_statistic_t::get_sld_total_energy(){

   return sld_total_energy;

}

//------------------------------------------------------------------------------------------------------
// Function to get const reference for total energy data (in Tesla)
//------------------------------------------------------------------------------------------------------
const std::vector<double>& sld_energy_statistic_t::get_sld_exchange_energy(){

   return sld_exchange_energy;

}

//------------------------------------------------------------------------------------------------------
// Function to get const reference for total energy data
//------------------------------------------------------------------------------------------------------
const std::vector<double>& sld_energy_statistic_t::get_sld_coupling_energy(){

   return sld_coupling_energy;

}

//------------------------------------------------------------------------------------------------------
// Function to get const reference for total energy data
//------------------------------------------------------------------------------------------------------
const std::vector<double>& sld_energy_statistic_t::get_potential_energy(){

   return potential_energy;

}

//------------------------------------------------------------------------------------------------------
// Function to get const reference for total energy data
//------------------------------------------------------------------------------------------------------
const std::vector<double>& sld_energy_statistic_t::get_kinetic_energy(){

   return kinetic_energy;

}

//------------------------------------------------------------------------------------------------------
// Function to set total energy data
//------------------------------------------------------------------------------------------------------
void sld_energy_statistic_t::set_sld_total_energy(std::vector<double>& new_energy, std::vector<double>& new_mean_energy){

   // copy energy vector
   sld_total_energy = new_energy;

   const unsigned int array_size = mean_sld_total_energy.size();
   for( int i=0; i < array_size; ++i ){
      mean_sld_total_energy[i] += new_mean_energy[i];
   }

}

//------------------------------------------------------------------------------------------------------
// Function to set exchange energy data
//------------------------------------------------------------------------------------------------------
void sld_energy_statistic_t::set_sld_exchange_energy(std::vector<double>& new_energy, std::vector<double>& new_mean_energy){

   // copy energy vector
   sld_exchange_energy = new_energy;

   const unsigned int array_size = mean_sld_exchange_energy.size();
   for( int i=0; i < array_size; ++i ){
      mean_sld_exchange_energy[i] += new_mean_energy[i];
   }

}

//------------------------------------------------------------------------------------------------------
// Function to set anisotropy energy data
//------------------------------------------------------------------------------------------------------
void sld_energy_statistic_t::set_sld_coupling_energy(std::vector<double>& new_energy, std::vector<double>& new_mean_energy){

   // copy energy vector
   sld_coupling_energy = new_energy;

   const unsigned int array_size = mean_sld_coupling_energy.size();
   for( int i=0; i < array_size; ++i ){
      mean_sld_coupling_energy[i] += new_mean_energy[i];
   }

}

//------------------------------------------------------------------------------------------------------
// Function to set applied field energy data
//------------------------------------------------------------------------------------------------------
void sld_energy_statistic_t::set_potential_energy(std::vector<double>& new_energy, std::vector<double>& new_mean_energy){

   // copy energy vector
   potential_energy = new_energy;

   const unsigned int array_size = mean_potential_energy.size();
   for( int i=0; i < array_size; ++i ){
      mean_potential_energy[i] += new_mean_energy[i];
   }

}

//------------------------------------------------------------------------------------------------------
// Function to set magnetostatic energy data
//------------------------------------------------------------------------------------------------------
void sld_energy_statistic_t::set_kinetic_energy(std::vector<double>& new_energy, std::vector<double>& new_mean_energy){

   // copy energy vector
   kinetic_energy = new_energy;

   const unsigned int array_size = mean_kinetic_energy.size();
   for( int i=0; i < array_size; ++i ){
      mean_kinetic_energy[i] += new_mean_energy[i];
   }

}

//------------------------------------------------------------------------------------------------------
// Function to update mean counter
//------------------------------------------------------------------------------------------------------
void sld_energy_statistic_t::update_mean_counter(long counter){

   // update counter
   mean_counter += double(counter);

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to reset energy averages
//------------------------------------------------------------------------------------------------------
void sld_energy_statistic_t::reset_averages(){

   // reinitialise mean energies to zero
   std::fill(        mean_sld_total_energy.begin(),         mean_sld_total_energy.end(), 0.0);
   std::fill(     mean_sld_exchange_energy.begin(),      mean_sld_exchange_energy.end(), 0.0);
   std::fill(   mean_sld_coupling_energy.begin(),    mean_sld_coupling_energy.end(), 0.0);
   std::fill(mean_kinetic_energy.begin(), mean_kinetic_energy.end(), 0.0);
   std::fill(mean_potential_energy.begin(), mean_potential_energy.end(), 0.0);

   // reset data counter
   mean_counter = 0.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised magnetisation values as string
//------------------------------------------------------------------------------------------------------
std::string sld_energy_statistic_t::output_sld_energy(enum sld_energy_t sld_energy_type,bool header){

   // result string stream
   std::ostringstream res;
   vout::fixed_width_output result(res,vout::fw_size); 
   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      //if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
      if(vout::fixed) res.setf(std::ios::scientific);
   }

   // could specify an output conversion unit here
   // could run all of the for loops as one outside the switch. -AJN
   if(header){
       for(int mask_id=0;mask_id<mask_size; ++mask_id){
           result<<name + std::to_string(mask_id) + "_E"+std::to_string(sld_energy_type) + "_Energy";
       }
   }else{
   // output correct energy type (in Joules)
       switch(sld_energy_type){
          case sld_total :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result <<         sld_total_energy[mask_id];
             break;

          case sld_exchange :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result <<      sld_exchange_energy[mask_id];
             break;

          case sld_coupling :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result <<   sld_coupling_energy[mask_id] ;
             break;

          case potential :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result << potential_energy[mask_id] ;
             break;

          case kinetic :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result << kinetic_energy[mask_id] ;
             break;

          default :
             break;
       }

   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised magnetisation values as string
//------------------------------------------------------------------------------------------------------
std::string sld_energy_statistic_t::output_mean_sld_energy(enum sld_energy_t sld_energy_type,bool header){

   // result string stream
   std::ostringstream res;
   vout::fixed_width_output result(res,vout::fw_size); 
   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      //if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
      if(vout::fixed) res.setf(std::ios::scientific);

   }

   // could specify an output conversion unit here
   // could run all of the for loops as one outside the switch. -AJN
   if(header){
       for(int mask_id=0;mask_id<mask_size; ++mask_id){
           result << name + std::to_string(mask_id) + "_E" + std::to_string(sld_energy_type) + "_Energy";
       }
   }else{
   // output correct energy type (in Joules)
       switch(sld_energy_type){

          case sld_total :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result <<         mean_sld_total_energy[mask_id] / mean_counter;
             break;

          case sld_exchange :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result <<      mean_sld_exchange_energy[mask_id] / mean_counter;
             break;

          case sld_coupling :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result <<    mean_sld_coupling_energy[mask_id] / mean_counter;
             break;

          case potential :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result << mean_potential_energy[mask_id]  / mean_counter;
             break;

          case kinetic :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result << mean_kinetic_energy[mask_id]  / mean_counter;
             break;

          default :
             break;

       }
   }

   return result.str();

}

} // end of namespace stats
