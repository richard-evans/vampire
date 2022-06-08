//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2021. All rights reserved.
// {J L Ross 2021}
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>


// Vampire headers
#include "errors.hpp"
#include "stats.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
#include "stopwatch.hpp"
#include "atoms.hpp"
namespace stats{

//-------------------------------------------
// Constructor
//---------------------------------------------

//spin_temperature_statistic_t::spin_temperature_statistic_t (){}

//
// Function to initialize the data structures
//--------------------------------------------


void spin_temperature_statistic_t::set_mask(const int mask_length, const std::vector<int>& mag_mask) {
       //check magnetisation statistic
    if(mask_length == 0) {
        terminaltextcolor(RED);
        std::cerr << "Programmer Error - Uninitilized magnetism statistic mask passed to module" << std::endl;
        terminaltextcolor(WHITE);
        zlog << zTs() << "Programmer Error - Uninitilized magnetism statistic mask passed to spin-temp module" << std::endl;
        err::vexit();
    }

    // save mask to internal storage
    num_atoms = mag_mask.size();
    mask_size = mask_length;
    mask = mag_mask;

     // Check that mask values never exceed mask_size
   for(unsigned int atom=0; atom < num_atoms; ++atom){
      if(mask[atom] > mask_size-1) {
         terminaltextcolor(RED);
         std::cerr << "Programmer Error - mask id " << mask[atom] << " is greater than number of elements for mask "<< mask_size << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Programmer Error - mask id " << mask[atom] << " is greater than number of elements for mask "<< mask_size << std::endl;
         err::vexit();
      }
   }


    //resize mask arrays
    total_spin_temperature.resize(mask_size, 0.0);
    spin_temperature_bottom.resize(mask_size, 0.0);
    spin_temperature_top.resize(mask_size, 0.0);
    mean_spin_temperature.resize(mask_size, 0.0);
    mean_spin_counter.resize(mask_size, 0.0);

       // determine mask id's with no atoms
   std::vector<int> num_atoms_in_mask(mask_size,0);
   for(unsigned int atom=0; atom<num_atoms; ++atom){
      int mask_id = mask[atom];
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
    
    initialized = true;
   
   
    if(!initialized) {
        terminaltextcolor(RED);
        std::cerr << "Programmer Error - Uninitilized spin-temp masks" << std::endl;
        terminaltextcolor(WHITE);
        zlog << zTs() << "Programmer Error - Uninitilized spin-temp masks" << std::endl;
        err::vexit();
    }
   return;

}

//Function to calculate the spin temperature according to the equation derived by Ma et al., 2010
// DOI: 10.1103/PhysRevE.82.031111
//---------------------------------
//
//
//
//
//        T_s = sum_i (S_i X H_i)**2
//                 -----------------
//         2 k_B  sum_i S_i [dot] H_i
//
//----------------------------------

void spin_temperature_statistic_t::calculate(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz, const std::vector<double>& sm, 
    const std::vector<double>& Hx_int, const std::vector<double>& Hy_int,const std::vector<double>& Hz_int,
    const std::vector<double>& Hx_ext, const std::vector<double>& Hy_ext, const std::vector<double>& Hz_ext) {
 

    
    //zero spin_temperature array
    for (int mask_num = 0; mask_num < mask_size; ++mask_num ) {
        spin_temperature_top[mask_num] = spin_temperature_bottom[mask_num] = total_spin_temperature[mask_num] = 0;

    }
  
   //create local variables
    double S[3];
    double H[3]; 
    
    
        int mask_num = 0;
        for(int atom = 0; atom < num_atoms; ++atom) {
            mask_num = mask[atom];
            
                //get spin vector and moment per atom
             const double   S[3] = {sx[atom] * sm[atom], 
                        sy[atom] * sm[atom],
                        sz[atom] * sm[atom]};

                //get field value for each atom (note the exclusion of the thermal field)
              const double  H[3] = {Hx_int[atom] + Hx_ext[atom], // + x_thermal[atom];
                        Hy_int[atom] + Hy_ext[atom],// + y_thermal[atom];
                        Hz_int[atom] + Hz_ext[atom]};// + z_thermal[atom];
              
                 //crossproduct**2 can be simplified using Lagrange's Identity
                 // 
                 //  || A x B ||**2 = (a1b2 - a2b1)**2 + (a2b3 - a3b2)**2 + (a3b1 - a1b3)**2
                 
                spin_temperature_top[mask_num]  += (((S[0]*H[1])-(S[1]*H[0])) * ((S[0]*H[1]) - (S[1]*H[0])))  
                + (((S[1]*H[2])-(S[2]*H[1])) * ((S[1]*H[2]) - (S[2]*H[1]))) 
                + (((S[2]*H[0])-(S[0]*H[2])) * ((S[2]*H[0]) - (S[0]*H[2])));

                //dot product
                spin_temperature_bottom[mask_num] += (S[0]*H[0])+(S[1]*H[1])+(S[2]*H[2]);
            
        } 
       //reduction
        #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &total_spin_temperature[0], mask_size,  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        #endif

        //normalize to J then K. 
        const double u_B = 9.274009994e-24;
        const double K_b = 1.38064852e-23;
         for (int mask_num = 0; mask_num < mask_size; ++mask_num) {  
            total_spin_temperature[mask_num] = (u_B * spin_temperature_top[mask_num]) / (2 * K_b * spin_temperature_bottom[mask_num]); 
            mean_spin_counter[mask_num] += 1;
            mean_spin_temperature[mask_num] += total_spin_temperature[mask_num];
         }
          // Zero empty mask id's
         for(unsigned int id=0; id < zero_list.size(); ++id) total_spin_temperature[zero_list[id]]=0.0;

    }  

    
//output function
std::string spin_temperature_statistic_t::output_spin_temperature(bool header) {
    
    std::ostringstream res;

    if(vout::custom_precision) {
        res.precision(vout::precision);
        if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
    }
    vout::fixed_width_output result(res,vout::fw_size);

    for (int mask_id = 0; mask_id < mask_size; ++mask_id) {  

       if (header) {
            result << name + std::to_string(mask_id) + "spinTemp";
        } else {
            result << total_spin_temperature[mask_id];
        }
    
    }return result.str();
}

std::string spin_temperature_statistic_t::output_mean_spin_temperature(bool header) {
    std::ostringstream res;

    if(vout::custom_precision) {
        res.precision(vout::precision);
        if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
    }
    
    vout::fixed_width_output result(res,vout::fw_size); 

   // loop over all masks
   for(int mask_id = 0; mask_id < mask_size; ++mask_id){
       if(header) {
           result << name + std::to_string(mask_id) + "mean_spinTemp";
       } else{
        double ic = 1.0 / mean_spin_counter[mask_id];
        const double temp = mean_spin_temperature[mask_id] * ic;
          result << temp;
        
       }
    }
    return result.str();
}

//reset function at end of statistic
void spin_temperature_statistic_t::reset() {

    fill(total_spin_temperature.begin(), total_spin_temperature.end(), 0.0);
    fill(spin_temperature_bottom.begin(), spin_temperature_bottom.end(), 0.0);
    fill(spin_temperature_top.begin(), spin_temperature_top.end(), 0.0);
    fill(mean_spin_temperature.begin(), mean_spin_temperature.end(), 0.0);
    fill(mean_spin_counter.begin(), mean_spin_counter.end(), 0.0);
    
}

}

