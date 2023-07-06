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
bool energy_statistic_t::is_initialized(){
   return initialized;
}

//------------------------------------------------------------------------------------------------------
// Function to initialize mask
//------------------------------------------------------------------------------------------------------
void energy_statistic_t::set_mask(const int in_mask_size, const std::vector<int> in_mask){

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
   total_energy.resize(in_mask_size, 0.0);
   exchange_energy.resize(in_mask_size, 0.0);
   anisotropy_energy.resize(in_mask_size, 0.0);
   applied_field_energy.resize(in_mask_size, 0.0);
   magnetostatic_energy.resize(in_mask_size, 0.0);

   mean_total_energy.resize(in_mask_size, 0.0);
   mean_exchange_energy.resize(in_mask_size, 0.0);
   mean_anisotropy_energy.resize(in_mask_size, 0.0);
   mean_applied_field_energy.resize(in_mask_size, 0.0);
   mean_magnetostatic_energy.resize(in_mask_size, 0.0);

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
   initialized = true;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to get mask needed for gpu acceleration of statistics calculation
//------------------------------------------------------------------------------------------------------
void energy_statistic_t::get_mask(std::vector<int>& out_mask, std::vector<double>& out_normalisation){

   // copy data to objects
   out_mask = mask;
   out_normalisation = normalisation;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to calculate spin energy given a mask and place result in energy array
//------------------------------------------------------------------------------------------------------
void energy_statistic_t::calculate(const std::vector<double>& sx,  // spin unit vector
                                   const std::vector<double>& sy,
                                   const std::vector<double>& sz,
                                   const std::vector<double>& mm,  // magnetic moment (Tesla)
                                   const std::vector<int>& mat, // material id
                                   const double temperature){

   // initialise energies to zero
   std::fill(         total_energy.begin(),         total_energy.end(), 0.0 );
   std::fill(      exchange_energy.begin(),      exchange_energy.end(), 0.0 );
   std::fill(    anisotropy_energy.begin(),    anisotropy_energy.end(), 0.0 );
   std::fill( applied_field_energy.begin(), applied_field_energy.end(), 0.0 );
   std::fill( magnetostatic_energy.begin(), magnetostatic_energy.end(), 0.0 );

   //---------------------------------------------------------------------------
   // Calculate exchange energy (in Tesla)
   //---------------------------------------------------------------------------

   // loop over all atoms in mask
   for( int atom = 0; atom < num_atoms; ++atom ){
      const int mask_id = mask[atom]; // get mask id
      exchange_energy[mask_id] += exchange::single_spin_energy(atom, sx[atom], sy[atom], sz[atom]) * mm[atom];
   }

   // Optionally calculate biquadratic exchange energy
   if(exchange::biquadratic){

      // loop over all atoms in mask
      for( int atom = 0; atom < num_atoms; ++atom ){
         const int mask_id = mask[atom]; // get mask id
         exchange_energy[mask_id] += exchange::single_spin_biquadratic_energy(atom, sx[atom], sy[atom], sz[atom]) * mm[atom];
      }

   }

   // save total energy accounting for factor 1/2 in double summation
   for( int mask_id = 0; mask_id < mask_size; ++mask_id ){
      exchange_energy[mask_id] = 0.5 * exchange_energy[mask_id];
   }

   //---------------------------------------------------------------------------
   // Calculate anisotropy energy (in Tesla)
   //---------------------------------------------------------------------------
   for( int atom = 0; atom < num_atoms; ++atom ){
      const int mask_id = mask[atom]; // get mask id
      anisotropy_energy[mask_id] += anisotropy::single_spin_energy(atom, mat[atom], sx[atom], sy[atom], sz[atom], temperature) * mm[atom];
   }

   //---------------------------------------------------------------------------
   // Calculate applied field energy (in Tesla)
   //---------------------------------------------------------------------------
   for( int atom = 0; atom < num_atoms; ++atom ){
      const int mask_id = mask[atom]; // get mask id
      applied_field_energy[mask_id] += sim::spin_applied_field_energy(sx[atom], sy[atom], sz[atom]) * mm[atom];
   }

   //---------------------------------------------------------------------------
   // Calculate magnetostatic field energy (in Tesla)
   //---------------------------------------------------------------------------
   for( int atom = 0; atom < num_atoms; ++atom ){
      const int mask_id = mask[atom]; // get mask id
      magnetostatic_energy[mask_id] += dipole::spin_magnetostatic_energy(atom, sx[atom], sy[atom], sz[atom]) * mm[atom];
   }

   // save energy accounting for factor 1/2 in double summation
   for( int mask_id = 0; mask_id < mask_size; ++mask_id ){
      magnetostatic_energy[mask_id] = 0.5 * magnetostatic_energy[mask_id];
   }

   //---------------------------------------------------------------------------
   // Calculate total energy (in Tesla)
   //---------------------------------------------------------------------------
   for( int mask_id = 0; mask_id < mask_size; ++mask_id ){
      total_energy[mask_id] = exchange_energy[mask_id] +
                              anisotropy_energy[mask_id] +
                              applied_field_energy[mask_id] +
                              magnetostatic_energy[mask_id];
   }

   //---------------------------------------------------------------------------
   // Reduce on all CPUS
   //---------------------------------------------------------------------------
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE,      &exchange_energy[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE,    &anisotropy_energy[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &applied_field_energy[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &magnetostatic_energy[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE,         &total_energy[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

   //---------------------------------------------------------------------------
   // Add energies to mean energies
   //---------------------------------------------------------------------------
   for(int mask_id=0; mask_id<mask_size; ++mask_id ){
      mean_total_energy[mask_id]         += total_energy[mask_id];
      mean_exchange_energy[mask_id]      += exchange_energy[mask_id];
      mean_anisotropy_energy[mask_id]    += anisotropy_energy[mask_id];
      mean_applied_field_energy[mask_id] += applied_field_energy[mask_id];
      mean_magnetostatic_energy[mask_id] += magnetostatic_energy[mask_id];
   }

   // increment mean counter
   mean_counter += 1.0;

   //---------------------------------------------------------------------------
   // Zero empty mask id's
   //---------------------------------------------------------------------------
   for( unsigned int id=0; id < zero_list.size(); ++id ){
              total_energy[ zero_list[id] ] = 0.0;
           exchange_energy[ zero_list[id] ] = 0.0;
         anisotropy_energy[ zero_list[id] ] = 0.0;
      applied_field_energy[ zero_list[id] ] = 0.0;
      magnetostatic_energy[ zero_list[id] ] = 0.0;
   }

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to get const reference for total energy data
//------------------------------------------------------------------------------------------------------
const std::vector<double>& energy_statistic_t::get_total_energy(){

   return total_energy;

}

//------------------------------------------------------------------------------------------------------
// Function to get const reference for total energy data (in Tesla)
//------------------------------------------------------------------------------------------------------
const std::vector<double>& energy_statistic_t::get_exchange_energy(){

   return exchange_energy;

}

//------------------------------------------------------------------------------------------------------
// Function to get const reference for total energy data
//------------------------------------------------------------------------------------------------------
const std::vector<double>& energy_statistic_t::get_anisotropy_energy(){

   return anisotropy_energy;

}

//------------------------------------------------------------------------------------------------------
// Function to get const reference for total energy data
//------------------------------------------------------------------------------------------------------
const std::vector<double>& energy_statistic_t::get_applied_field_energy(){

   return applied_field_energy;

}

//------------------------------------------------------------------------------------------------------
// Function to get const reference for total energy data
//------------------------------------------------------------------------------------------------------
const std::vector<double>& energy_statistic_t::get_magnetostatic_energy(){

   return magnetostatic_energy;

}

//------------------------------------------------------------------------------------------------------
// Function to set total energy data
//------------------------------------------------------------------------------------------------------
void energy_statistic_t::set_total_energy(std::vector<double>& new_energy, std::vector<double>& new_mean_energy){

   // copy energy vector
   total_energy = new_energy;

   const size_t array_size = mean_total_energy.size();
   for( size_t i=0; i < array_size; ++i ){
      mean_total_energy[i] += new_mean_energy[i];
   }

}

//------------------------------------------------------------------------------------------------------
// Function to set exchange energy data
//------------------------------------------------------------------------------------------------------
void energy_statistic_t::set_exchange_energy(std::vector<double>& new_energy, std::vector<double>& new_mean_energy){

   // copy energy vector
   exchange_energy = new_energy;

   const size_t array_size = mean_exchange_energy.size();
   for( size_t i=0; i < array_size; ++i ){
      mean_exchange_energy[i] += new_mean_energy[i];
   }

}

//------------------------------------------------------------------------------------------------------
// Function to set anisotropy energy data
//------------------------------------------------------------------------------------------------------
void energy_statistic_t::set_anisotropy_energy(std::vector<double>& new_energy, std::vector<double>& new_mean_energy){

   // copy energy vector
   anisotropy_energy = new_energy;

   const size_t array_size = mean_anisotropy_energy.size();
   for( size_t i=0; i < array_size; ++i ){
      mean_anisotropy_energy[i] += new_mean_energy[i];
   }

}

//------------------------------------------------------------------------------------------------------
// Function to set applied field energy data
//------------------------------------------------------------------------------------------------------
void energy_statistic_t::set_applied_field_energy(std::vector<double>& new_energy, std::vector<double>& new_mean_energy){

   // copy energy vector
   applied_field_energy = new_energy;

   const size_t array_size = mean_applied_field_energy.size();
   for( size_t i=0; i < array_size; ++i ){
      mean_applied_field_energy[i] += new_mean_energy[i];
   }

}

//------------------------------------------------------------------------------------------------------
// Function to set magnetostatic energy data
//------------------------------------------------------------------------------------------------------
void energy_statistic_t::set_magnetostatic_energy(std::vector<double>& new_energy, std::vector<double>& new_mean_energy){

   // copy energy vector
   magnetostatic_energy = new_energy;

   const size_t array_size = mean_magnetostatic_energy.size();
   for( size_t i=0; i < array_size; ++i ){
      mean_magnetostatic_energy[i] += new_mean_energy[i];
   }

}

//------------------------------------------------------------------------------------------------------
// Function to update mean counter
//------------------------------------------------------------------------------------------------------
void energy_statistic_t::update_mean_counter(long counter){

   // update counter
   mean_counter += double(counter);

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to reset energy averages
//------------------------------------------------------------------------------------------------------
void energy_statistic_t::reset_averages(){

   // reinitialise mean energies to zero
   std::fill(        mean_total_energy.begin(),         mean_total_energy.end(), 0.0);
   std::fill(     mean_exchange_energy.begin(),      mean_exchange_energy.end(), 0.0);
   std::fill(   mean_anisotropy_energy.begin(),    mean_anisotropy_energy.end(), 0.0);
   std::fill(mean_applied_field_energy.begin(), mean_applied_field_energy.end(), 0.0);
   std::fill(mean_magnetostatic_energy.begin(), mean_magnetostatic_energy.end(), 0.0);

   // reset data counter
   mean_counter = 0.0;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to output normalised magnetisation values as string
//------------------------------------------------------------------------------------------------------
std::string energy_statistic_t::output_energy(enum energy_t energy_type,bool header){

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
           result<<name + std::to_string(mask_id) + "_E"+std::to_string(energy_type) + "_Energy";
       }
   }else{
   // output correct energy type (in Joules)
       switch(energy_type){
          case total :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result <<         total_energy[mask_id] * constants::muB;
             break;

          case exchange :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result <<      exchange_energy[mask_id] * constants::muB;
             break;

          case anisotropy :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result <<    anisotropy_energy[mask_id] * constants::muB;
             break;

          case applied_field :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result << applied_field_energy[mask_id] * constants::muB;
             break;

          case magnetostatic :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result << magnetostatic_energy[mask_id] * constants::muB;
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
std::string energy_statistic_t::output_mean_energy(enum energy_t energy_type,bool header){

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
           result << name + std::to_string(mask_id) + "_E" + std::to_string(energy_type) + "_Energy";
       }
   }else{
   // output correct energy type (in Joules)
       switch(energy_type){

          case total :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result <<         mean_total_energy[mask_id] * constants::muB / mean_counter;
             break;

          case exchange :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result <<      mean_exchange_energy[mask_id] * constants::muB / mean_counter;
             break;

          case anisotropy :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result <<    mean_anisotropy_energy[mask_id] * constants::muB / mean_counter;
             break;

          case applied_field :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result << mean_applied_field_energy[mask_id] * constants::muB / mean_counter;
             break;

          case magnetostatic :
             for(int mask_id=0; mask_id<mask_size; ++mask_id) result << mean_magnetostatic_energy[mask_id] * constants::muB / mean_counter;
             break;

          default :
             break;

       }
   }

   return result.str();

}

} // end of namespace stats
