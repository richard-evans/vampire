//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//   implementation based on the paper Phys. Rev. B 103, 024429, (2021) M.Strungaru, M.O.A. Ellis et al
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

// Vampire headers
#include "sld.hpp"
#include "create.hpp"
#include "sim.hpp"
#include "errors.hpp"
#include "material.hpp"
#include <fstream>

// sld module headers
#include "internal.hpp"


namespace sld{

   void compute_forces(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& neighbour_list_start_index,
               const std::vector<int>& neighbour_list_end_index,
               const std::vector<int>& type_array, // type for atom
               const std::vector<int>& neighbour_list_array, // list of interactions between atom
               const std::vector<double>& x0_coord_array, // coord vectors for atoms
               const std::vector<double>& y0_coord_array,
               const std::vector<double>& z0_coord_array,
               const std::vector<double>& x_coord_array, // coord vectors for atoms
               const std::vector<double>& y_coord_array,
               const std::vector<double>& z_coord_array,
               std::vector<double>& forces_array_x, //  vectors for forces
               std::vector<double>& forces_array_y,
               std::vector<double>& forces_array_z,
               std::vector<double>& potential_eng){



      // calculate harmonic potential forces
      if(sld::internal::harmonic){

        //std::cout<<"inside compute function harmonic" <<std::endl;

              internal::compute_forces_harmonic(start_index, end_index,
                                                         neighbour_list_start_index, neighbour_list_end_index,
                                                         type_array, neighbour_list_array,
                                                         x0_coord_array, y0_coord_array, z0_coord_array,
                                                         x_coord_array, y_coord_array, z_coord_array,
                                                         forces_array_x, forces_array_y, forces_array_z, potential_eng);
      }
      // calculate Morse potential forces
      if(sld::internal::morse){

                  internal::compute_forces_morse(start_index, end_index,
                                                             neighbour_list_start_index, neighbour_list_end_index,
                                                             type_array, neighbour_list_array,
                                                             x_coord_array, y_coord_array, z_coord_array,
                                                             forces_array_x, forces_array_y, forces_array_z, potential_eng);
      }
      
      // calculate THz pulse force
      if(sld::internal::linear_pump_enabled){

        internal::compute_thz(start_index, end_index,
                                x_coord_array, y_coord_array, z_coord_array,
                                forces_array_x, forces_array_y, forces_array_z);
      
        }

      return;

}

double PBC_wrap ( double dx, double L, bool bounds){
    return (bounds) ? dx - floor( (dx/L) + 0.5) * L : dx;
}


namespace internal{


void compute_forces_harmonic(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x0_coord_array, // coord vectors for atoms
            const std::vector<double>& y0_coord_array,
            const std::vector<double>& z0_coord_array,
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng){



            double rx, ry, rz;
            double rx0, ry0, rz0;
            double dx, dy, dz;
            double dx0, dy0, dz0;
            double fx = 0.0, fy = 0.0, fz = 0.0;
            double rji_sqr, rji, rji0, inv_rji;
            int j;
            //int total_int;
            double r_sqr_cut=sld::internal::r_cut_pot*sld::internal::r_cut_pot;
            double energy;



            for(int i=start_index;i< end_index; ++i){

               fx = 0.0;
               fy = 0.0;
               fz = 0.0;
               energy=0.0;
               //total_int=0;

               rx = x_coord_array[i];
               ry = y_coord_array[i];
               rz = z_coord_array[i];
               rx0 = x0_coord_array[i];
               ry0 = y0_coord_array[i];
               rz0 = z0_coord_array[i];

                //note for sld_neighbour_list_array
                // int nbr_end = neighbour_list_end_index[i];
                //for(int i=start_index;i<= end_index; ++i)
                //for( int n = nbr_start; n <=nbr_end; ++n)



        	      int nbr_start = neighbour_list_start_index[i];
        	      int nbr_end = neighbour_list_end_index[i]+1;

        	      for( int n = nbr_start; n < nbr_end; ++n){
        	        j = neighbour_list_array[n];

        	        if ( j != i){
        		       dx = -x_coord_array[j] + rx;
        		       dy = -y_coord_array[j] + ry;
        		       dz = -z_coord_array[j] + rz;
                   dx0 = -x0_coord_array[j] + rx0;
                   dy0 = -y0_coord_array[j] + ry0;
                   dz0 = -z0_coord_array[j] + rz0;

                   dx = sld::PBC_wrap( dx, cs::system_dimensions[0], cs::pbc[0]);
                   dy = sld::PBC_wrap( dy, cs::system_dimensions[1], cs::pbc[1]);
                   dz = sld::PBC_wrap( dz, cs::system_dimensions[2], cs::pbc[2]);
                   dx0 = sld::PBC_wrap( dx0, cs::system_dimensions[0], cs::pbc[0]);
                   dy0 = sld::PBC_wrap( dy0, cs::system_dimensions[1], cs::pbc[1]);
                   dz0 = sld::PBC_wrap( dz0, cs::system_dimensions[2], cs::pbc[2]);



        		       rji_sqr = dx*dx + dy*dy + dz*dz;

        		       if( rji_sqr < r_sqr_cut){

        		          //total_int++;


        		           rji = sqrt(rji_sqr);
                           rji0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
        		           inv_rji = 1.0/ rji;

                       energy += (rji-rji0)*(rji-rji0);


		               fx -=  (rji-rji0)*dx*inv_rji ; //2 (rji-rj0)*dx*inv_rji -> 2 went at the end
                       fy -=  (rji-rji0)*dy*inv_rji ;
                       fz -=  (rji-rji0)*dz*inv_rji ;


                     }
        	   	}
        	    }
             const unsigned int imat = type_array[i];

             double V0=sld::internal::mp[imat].V0.get(); //0.15

        	    forces_array_x[i] += V0 * 2.0 * fx;
        	    forces_array_y[i] += V0  * 2.0 * fy;
        	    forces_array_z[i] += V0  * 2.0 * fz;
                potential_eng[i] = 0.5 * V0 * energy;



  }


     return;
            }

void compute_forces_morse(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng){



            double rx, ry, rz;
            double dx, dy, dz;
            double fx = 0.0, fy = 0.0, fz = 0.0;
            double rji_sqr, rji, inv_rji; // rji0,
            int j; //, total_int;
            double r_sqr_cut=sld::internal::r_cut_pot*sld::internal::r_cut_pot;
            double energy;

           double alpha_m= sld::internal::alpha_m; //1.3885;
           //double r0_m= sld::internal::r0_m;//2.845;
           double morse_D=sld::internal::morse_D; //0.4174;;
           double morse_beta=sld::internal::morse_beta;//exp( alpha_m * r0_m);
           double morse_factor = sld::internal::morse_factor; //-2.0 * morse_D * alpha_m;




            for(int i=start_index;i< end_index; ++i){

               fx = 0.0;
               fy = 0.0;
               fz = 0.0;
               energy=0.0;
               //total_int=0;


               rx = x_coord_array[i];
               ry = y_coord_array[i];
               rz = z_coord_array[i];


        	      int nbr_start = neighbour_list_start_index[i];
        	      int nbr_end = neighbour_list_end_index[i]+1;

        	      for( int n = nbr_start; n < nbr_end; ++n){
        	        j = neighbour_list_array[n];

        	        if ( j != i){
        		       dx = x_coord_array[j] -rx;
        		       dy = y_coord_array[j]- ry;
        		       dz = z_coord_array[j]- rz;

                   dx = sld::PBC_wrap( dx, cs::system_dimensions[0], cs::pbc[0]);
                   dy = sld::PBC_wrap( dy, cs::system_dimensions[1], cs::pbc[1]);
                   dz = sld::PBC_wrap( dz, cs::system_dimensions[2], cs::pbc[2]);



        		       rji_sqr = dx*dx + dy*dy + dz*dz;

        		       if( rji_sqr < r_sqr_cut){

        		           rji = sqrt(rji_sqr);
    		               inv_rji = 1.0/ rji;

		                   double y = morse_beta * exp( - alpha_m * rji);
  		                   double f_morse = y * ( y - 1.0);

        		           //std::cout<<"embedded "<<i<<"\t"<<j<<"\t"<<rji<<"\t"<<position_rij<<"\t"<<int(position_rij)<<"\t"<<v_rij<<"\t"<<rho_rij<<std::endl;



		                 fx +=  f_morse*dx*inv_rji;
                         fy +=  f_morse*dy*inv_rji;
                         fz +=  f_morse*dz*inv_rji;

                       //if(i==0) std::cout<<"fxyz "<<i <<"\t"<<j<<"\t"<<rji<<"\t"<<f_morse*dx*inv_rji<<"\t"<<f_morse*dy*inv_rji<<"\t"<<f_morse*dz*inv_rji<< std::endl;//<<(rji-rji0)*dx*inv_rji<<"\t"<<(rji-rji0)*dy*inv_rji<<"\t"<<(rji-rji0)*dz*inv_rji<<std::endl;
                       //if(i==0) std::cout<<"fxyz "<<i <<"\t"<<j<<"\t"<<rji<<"\t"<<x_coord_array[j]<<"\t"<<y_coord_array[j]<<"\t"<<z_coord_array[j]<< "\t"<<dx<<"\t"<<dy<<"\t"<<dz<<"\t"<<f_morse<<std::endl;//<<(rji-rji0)*dx*inv_rji<<"\t"<<(rji-rji0)*dy*inv_rji<<"\t"<<(rji-rji0)*dz*inv_rji<<std::endl;

                       //if(i==100) std::cout<<"fxyz "<<i <<"\t"<<j<<"\t"<<type_array[i]<<"\t"<<inv_rji<<std::endl;//<<(rji-rji0)*dx*inv_rji<<"\t"<<(rji-rji0)*dy*inv_rji<<"\t"<<(rji-rji0)*dz*inv_rji<<std::endl;
                         energy += y * ( y - 2.0);

                     }


        	   	}


        	    }





        	    forces_array_x[i] += fx*morse_factor;
        	    forces_array_y[i] += fy*morse_factor;
        	    forces_array_z[i] += fz*morse_factor;

              potential_eng[i] = morse_D * energy;

  }



   return;
          }

void compute_thz(const int start_index,
                 const int end_index,
                 const std::vector<double>& x_coord_array,
                 const std::vector<double>& y_coord_array,
                 const std::vector<double>& z_coord_array,
                 std::vector<double>& forces_array_x,
                 std::vector<double>& forces_array_y,
                 std::vector<double>& forces_array_z)
{

   static std::ofstream force_debug_file("force_debug_file.txt", std::ios::out);

    // Calculate the Current Physical Time
    // Converts the current step number into the actual time in seconds.
    uint64_t current_step = sim::time; 
    double dt = mp::dt_SI;             // Get the duration of one step in seconds 
    double current_time = static_cast<double>(current_step) * dt; // Total elapsed time (seconds)


    if (current_time < sld::internal::phonon_pulse_start_time || current_time > sld::internal::phonon_pulse_end_time) {
        return; // The pulse is off, so do nothing.
    }
 
    // Calculate the Raw Force for Each Atom
    // Loop through each atom to calculate the
    // force from the THz pulse and add it to a running total for later averaging.
    double sumx = 0.0, sumy = 0.0, sumz = 0.0;
    std::vector<double> f_thz_temp_x(end_index - start_index); // Temporary storage
    std::vector<double> f_thz_temp_y(end_index - start_index);
    std::vector<double> f_thz_temp_z(end_index - start_index);

    // Pre-calculate angular frequency (ω = 2πν) for efficiency
    const double twopi_niu = sld::internal::phonon_frequency * 6.28318530718;
    for (int i = start_index; i < end_index; ++i) {


        // Calculate the position-dependent part of the wave (k·r)
        double kr = sld::internal::phonon_wavevector[0] * x_coord_array[i] +
                    sld::internal::phonon_wavevector[1] * y_coord_array[i] +
                    sld::internal::phonon_wavevector[2] * z_coord_array[i];

        // Calculate the argument for the main physics equation: cos(ωt - k·r)
        double arg = current_time * twopi_niu - kr;
        double cos_factor = std::cos(arg); // The oscillating part of the force
        double sin_factor = std::sin(arg);
        
        // Calculate the raw force and store it temporarily
        int local_index = i - start_index;
        f_thz_temp_x[local_index] = sld::internal::phonon_force_amplitude[0] * cos_factor;
        f_thz_temp_y[local_index] = sld::internal::phonon_force_amplitude[1] * sin_factor;
        f_thz_temp_z[local_index] = sld::internal::phonon_force_amplitude[2] * cos_factor;

        // Add to the running total
        sumx += f_thz_temp_x[local_index];
        sumy += f_thz_temp_y[local_index];
        sumz += f_thz_temp_z[local_index];

        if (i == 0 && current_step % 200 == 0) { // Using 200 to match your output:output-rate
            std::cout << "\n--- THz FORCE DEBUG (Step " << current_step << ") ---" << std::endl;
        std::cout << "  Time (ps)   : " << current_time * 1.0e12 << std::endl;
        std::cout << "  arg (rad)   : " << arg << std::endl;
        std::cout << "  cos(arg)    : " << cos_factor << std::endl;
        std::cout << std::scientific; // Switch to scientific notation for forces
        std::cout << "  Raw Force X : " << f_thz_temp_x[local_index] << " N" << std::endl;
        std::cout << "  Sum X : " << sumx << " N" << std::endl;

      }
      }

    // Calculate the Center-of-Mass Correction
    // Calculate the average force and will subtract it from each atom's individual force.
    int counter = end_index - start_index;
    double avg_fx = 0.0, avg_fy = 0.0, avg_fz = 0.0;
    if (counter > 0) {
        avg_fx = sumx / static_cast<double>(counter);
        avg_fy = sumy / static_cast<double>(counter);
        avg_fz = sumz / static_cast<double>(counter);
    }
  
    
    // Apply the Final, Corrected Force 
    // Loop through the atoms again, retrieve the
    // temporarily stored raw force, subtract the average
    for (int i = start_index; i < end_index; ++i) {
        int local_index = i - start_index;
        double corrected_fx = f_thz_temp_x[local_index] - avg_fx;
        double corrected_fy = f_thz_temp_y[local_index] - avg_fy;
        double corrected_fz = f_thz_temp_z[local_index] - avg_fz;

        // Add the final calculated force to the atom's total force
        forces_array_x[i] += corrected_fx;
        forces_array_y[i] += corrected_fy;
        forces_array_z[i] += corrected_fz;

        // Forced applied to each atom over the simulation period 
        if (i == 0) {
            force_debug_file << current_time << "\t" << i << "\t" << local_index << "\t" << corrected_fx << "\t" << f_thz_temp_x[local_index] << "\t" << x_coord_array[i] << "\t" << y_coord_array[i] << "\t" << z_coord_array[i] << std::endl;
        }

        if (i == 0 && current_step % 200 == 0) {
        std::cout << "  Corrected Force X : " << corrected_fx << " N" << std::endl;
       
      }
    }

    return;
}
}
} 