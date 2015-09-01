//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>
#include <algorithm>

// Vampire headers
#include "ltmp.hpp"
#include "random.hpp"

// Local temperature pulse headers
#include "internal.hpp"

namespace ltmp{

   //-----------------------------------------------------------------------------
   // Function for updating local temperature fields
   //-----------------------------------------------------------------------------
   void update_localised_temperature(const double time_from_start){

      // calculate local temperature
      if(ltmp::internal::gradient == false) ltmp::internal::calculate_local_temperature_pulse(time_from_start);
      else ltmp::internal::calculate_local_temperature_gradient();

      // store number of local atoms as local constant for compiler
      const int num_local_atoms = ltmp::internal::num_local_atoms;

      // Initialise thermal field random numbers
      generate (ltmp::internal::x_field_array.begin(),ltmp::internal::x_field_array.begin()+num_local_atoms, mtrandom::gaussian);
      generate (ltmp::internal::y_field_array.begin(),ltmp::internal::y_field_array.begin()+num_local_atoms, mtrandom::gaussian);
      generate (ltmp::internal::z_field_array.begin(),ltmp::internal::z_field_array.begin()+num_local_atoms, mtrandom::gaussian);

      // check for temperature rescaling
      if(ltmp::internal::temperature_rescaling){
         // calculate local thermal field for all atoms with rescaled temperature
         for(int atom=0; atom<ltmp::internal::num_local_atoms; ++atom) {
            const int cell = ltmp::internal::atom_temperature_index[atom]; /// get cell index for atom temperature (Te or Tp)
            const double rootT = ltmp::internal::root_temperature_array[cell]; /// get sqrt(T) for atom
            const double sigma = ltmp::internal::atom_sigma[atom]; /// unrolled list of thermal prefactor

            // Calculate temperature rescaling (using root_T for performance)
            const double alpha = ltmp::internal::atom_rescaling_alpha[atom];
            const double root_Tc = ltmp::internal::atom_rescaling_root_Tc[atom];
            // if T<Tc T/Tc = (T/Tc)^alpha else T = T
            const double rescaled_rootT = rootT < root_Tc ? root_Tc*pow(rootT/root_Tc,alpha) : rootT;

            ltmp::internal::x_field_array[atom]*= sigma*rescaled_rootT;
            ltmp::internal::y_field_array[atom]*= sigma*rescaled_rootT;
            ltmp::internal::z_field_array[atom]*= sigma*rescaled_rootT;
         }
      }
      // otherwise use faster version without rescaling
      else{
         // calculate local thermal field for all atoms
         for(int atom=0; atom<ltmp::internal::num_local_atoms; ++atom) {
            const int cell = ltmp::internal::atom_temperature_index[atom]; /// get cell index for atom temperature (Te or Tp)
            const double rootT = ltmp::internal::root_temperature_array[cell]; /// get sqrt(T) for atom
            const double sigma = ltmp::internal::atom_sigma[atom]; /// unrolled list of thermal prefactor

            ltmp::internal::x_field_array[atom]*= sigma*rootT;
            ltmp::internal::y_field_array[atom]*= sigma*rootT;
            ltmp::internal::z_field_array[atom]*= sigma*rootT;
         }
      }

      return;
   }

   //-----------------------------------------------------------------------------
   // Function for adding local thermal fields to external field array
   //-----------------------------------------------------------------------------
   void get_localised_thermal_fields(std::vector<double>& x_total_external_field_array,
                               std::vector<double>& y_total_external_field_array,
                               std::vector<double>& z_total_external_field_array,
                               const int start_index,
                               const int end_index){

      // Add spin torque fields
      for(int i=start_index; i<end_index; ++i) x_total_external_field_array[i] += ltmp::internal::x_field_array[i];
      for(int i=start_index; i<end_index; ++i) y_total_external_field_array[i] += ltmp::internal::y_field_array[i];
      for(int i=start_index; i<end_index; ++i) z_total_external_field_array[i] += ltmp::internal::z_field_array[i];

      return;
   }

} // end of ltmp namespace
