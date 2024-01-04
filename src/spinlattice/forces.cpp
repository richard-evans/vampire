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

// Vampire headers
#include "sld.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include "create.hpp"
#include <iomanip>
#include <sim.hpp>

#include <errors.hpp>





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
            int j, total_int;
            double r_sqr_cut=sld::internal::r_cut_pot*sld::internal::r_cut_pot;
            double energy;



            for(int i=start_index;i< end_index; ++i){
            
               fx = 0.0;
               fy = 0.0;
               fz = 0.0;
               energy=0.0;
               total_int=0;

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
        		          
        		          total_int++;


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
         } //end of internal
      } // end of sld namespace
