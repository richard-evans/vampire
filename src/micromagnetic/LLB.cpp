//-----------------------------------------------------------------------------
//
//  Vampire - A code for cellistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
///====================================================================================================
///
///       				                    	LLB
///
///  			 Subroutine to simulate an cellistic system with LLB integration scheme
///
///									Version 1.0 R Evans 02/10/2008
///
///====================================================================================================


// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

#include "random.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

namespace mm = micromagnetic::internal;

int LLB_serial_heun(const int);

namespace micromagnetic{
/// Master LLB Function - dispatches code path to desired LLB routine
/// \f$ \frac{\partial S}{\partial t} \f$
int LLB(const int num_steps){
   //----------------------------------------------------------
   // check calling of routine if error checking is activated
   //----------------------------------------------------------


   //if(err::check==true){std::cout << "LLB has been called" << std::endl;}

   #ifdef MPICF
      //LLB_mpi(num_steps);
   #else
      LLB_serial_heun(num_steps);
   #endif

   return 0;
}
}


/// Performs serial Heun integration of the Landau-Lifshitz-Bloch Equation of motion
int LLB_serial_heun(
   const int num_steps,
   const int num_cells,
   const double temperature,
   std::vector<double> x_cell_array,
   std::vector<double> y_cell_array,
   std::vector<double> z_cell_array,
   std::vector<double> x_total_cell_field_array,
   std::vector<double> y_total_cell_field_array,
   std::vector<double> z_total_cell_field_array,
   std::vector<double> x_total_external_field_array,
   std::vector<double> y_total_external_field_array,
   std::vector<double> z_total_external_field_array,
   double Hx,
   double Hy,
   double Hz
){

   std::vector<double> alpha_para(num_cells,0.0);
   std::vector<double> alpha_perp(num_cells,0.0);
   std::vector<double> m_e(num_cells,0.0);
   std::vector<double> m_e_squared(num_cells,0.0);
   std::vector<double> one_o_chi_perp(num_cells,0.0);
   std::vector<double> one_o_2_chi_para(num_cells,0.0);
   std::vector<double> reduced_temperature(num_cells,0.0);
   std::vector<double> sigma_para(num_cells,0.0);
   std::vector<double> sigma_perp(num_cells,0.0);
   std::vector<double> Tc_o_Tc_m_T(num_cells,0.0);
   std::vector<double> dt(num_cells,0.0);

	const double kB = 1.3806503e-23;
   const double dt_SI = 1.0E-15;    //why?!?!?!?!?!?

   for (int cell =0; cell <num_cells; cell++)
   {
      reduced_temperature[cell] = temperature/mm::Tc[cell];
      Tc_o_Tc_m_T[cell] = mm::Tc[cell]/(temperature - mm::Tc[cell]);
   }

   for (int cell =0; cell <num_cells; cell++)
   {
      if (temperature<=mm::Tc[cell])
      {
         m_e[cell] = pow((mm::Tc[cell]-temperature)/(mm::Tc[cell]),0.365);
         alpha_para[cell] = (2.0/3.0)*mm::alpha[cell]*reduced_temperature[cell];
         alpha_perp[cell] = mm::alpha[cell]*(1.0-temperature/(3.0*mm::Tc[cell]));
      }
      else
      {
         m_e[cell] = 0.0;
         alpha_para[cell] = mm::alpha[cell]*(2.0/3.0)*reduced_temperature[cell];
         alpha_perp[cell] = alpha_para[cell];
      }
      m_e_squared[cell] = m_e[cell]*m_e[cell];
      one_o_chi_perp[cell] = 1.0/mm::chi_perp[cell];
      one_o_2_chi_para[cell] = 1.0/(2.0*mm::chi_para[cell]);
      dt[cell] = dt_SI*mm::gamma[cell];
   }

   for (int cell =0; cell <num_cells; cell++)
   {
	     if(temperature<0.1)   sigma_para[cell] = 1.0;
        else                  sigma_para[cell] = sqrt(2.0*kB*temperature/(mm::ms[cell]*mm::gamma[cell]*alpha_para[cell]*dt_SI));
        sigma_perp[cell] = sqrt(2.0*kB*temperature/(mm::ms[cell]*mm::gamma[cell]*alpha_perp[cell]*dt_SI));
   }

   //does mu_s*n_spins --> Ms ??
   //what are alpha/gamma/dt_SI


   //replace with an array over the cells?
   std::vector <double> Htx_perp(num_cells);
   std::vector <double> Hty_perp(num_cells);
   std::vector <double> Htz_perp(num_cells);
   std::vector <double> Htx_para(num_cells);
   std::vector <double> Hty_para(num_cells);
   std::vector <double> Htz_para(num_cells);

   // precalculate thermal fields
   generate (Htx_perp.begin(),Htx_perp.end(), mtrandom::gaussian);
   generate (Hty_perp.begin(),Hty_perp.end(), mtrandom::gaussian);
   generate (Htz_perp.begin(),Htz_perp.end(), mtrandom::gaussian);
   generate (Htx_para.begin(),Htx_para.end(), mtrandom::gaussian);
   generate (Hty_para.begin(),Hty_para.end(), mtrandom::gaussian);
   generate (Htz_para.begin(),Htz_para.end(), mtrandom::gaussian);


   std::vector<double> x_initial_cell_array(num_cells,0.0);
   std::vector<double> y_initial_cell_array(num_cells,0.0);
   std::vector<double> z_initial_cell_array(num_cells,0.0);
   std::vector<double> x_euler_array(num_cells,0.0);
   std::vector<double> y_euler_array(num_cells,0.0);
   std::vector<double> z_euler_array(num_cells,0.0);
   std::vector<double> x_heun_array(num_cells,0.0);
   std::vector<double> y_heun_array(num_cells,0.0);
   std::vector<double> z_heun_array(num_cells,0.0);
   std::vector<double> x_cell_storage_array(num_cells,0.0);
   std::vector<double> y_cell_storage_array(num_cells,0.0);
   std::vector<double> z_cell_storage_array(num_cells,0.0);
   std::vector<double> xyz(3,0.0);

   for(unsigned int cell=0;cell<num_cells;cell++){
		Htx_perp[cell] *= sigma_perp[cell];
		Hty_perp[cell] *= sigma_perp[cell];
		Htz_perp[cell] *= sigma_perp[cell];
		Htx_para[cell] *= sigma_para[cell];
		Hty_para[cell] *= sigma_para[cell];
		Htz_para[cell] *= sigma_para[cell];
	}

   for(int t=0;t<num_steps;t++){

      // Store initial spin positions for each cell
      for(unsigned int cell=0;cell<num_cells;cell++){
         x_initial_cell_array[cell] = x_cell_array[cell];
         y_initial_cell_array[cell] = y_cell_array[cell];
         z_initial_cell_array[cell] = z_cell_array[cell];
      }



      fill (x_total_cell_field_array.begin(),x_total_cell_field_array.end(),0.0);
      fill (y_total_cell_field_array.begin(),y_total_cell_field_array.end(),0.0);
      fill (z_total_cell_field_array.begin(),z_total_cell_field_array.end(),0.0);
      fill (x_total_external_field_array.begin(),x_total_external_field_array.end(),Hx);
      fill (y_total_external_field_array.begin(),y_total_external_field_array.end(),Hy);
      fill (z_total_external_field_array.begin(),z_total_external_field_array.end(),Hz);


      for(unsigned int cell=0;cell<num_cells;cell++){
			double m[3] = {x_cell_array[cell],y_cell_array[cell],z_cell_array[cell]};
			double m_squared = m[1]*m[1]+m[2]*m[2]+m[0]*m[0];
			double pf;
			if(temperature<=mm::Tc[cell]){
				pf = one_o_2_chi_para[cell]*(1.0 - m_squared/m_e_squared[cell]);
			}
			else{
				pf = -2.0*one_o_2_chi_para[cell]*(1.0 + Tc_o_Tc_m_T[cell]*3.0*m_squared/5.0);
			}

			x_total_cell_field_array[cell] = (pf-one_o_chi_perp[cell])*m[0];
			y_total_cell_field_array[cell] = (pf-one_o_chi_perp[cell])*m[1];
			z_total_cell_field_array[cell] = (pf-0.0				 )*m[2];
		}


      		// Calculate Euler Step
      		for(unsigned int cell=0;cell<num_cells;cell++){

      			// Store local spin in Sand local field in H
      			const double S[3] = {x_cell_array[cell],y_cell_array[cell],z_cell_array[cell]};

      			const double H[3] = {x_total_cell_field_array[cell]+x_total_external_field_array[cell],
      										y_total_cell_field_array[cell]+y_total_external_field_array[cell],
      										z_total_cell_field_array[cell]+z_total_external_field_array[cell]};

      			const double H_perp[3]={H[0]+Htx_perp[cell], H[1]+Hty_perp[cell], H[2]+Htz_perp[cell]};
      			const double H_para[3]={H[0]+Htx_para[cell], H[1]+Hty_para[cell], H[2]+Htz_para[cell]};

      			const double one_o_m_squared = 1.0/(S[1]*S[1]+S[2]*S[2]+S[0]*S[0]);

      			// Calculate Delta S
      			xyz[0]= 	-(S[1]*H[2]-S[2]*H[1])
      						+ alpha_para[cell]*S[0]*S[0]*H_para[0]*one_o_m_squared
      						-alpha_perp[cell]*(S[1]*(S[0]*H_perp[1]-S[1]*H_perp[0])-S[2]*(S[2]*H_perp[0]-S[0]*H_perp[2]))*one_o_m_squared;

      			xyz[1]= 	-(S[2]*H[0]-S[0]*H[2])
      						+ alpha_para[cell]*S[1]*S[1]*H_para[1]*one_o_m_squared
      						-alpha_perp[cell]*(S[2]*(S[1]*H_perp[2]-S[2]*H_perp[1])-S[0]*(S[0]*H_perp[1]-S[1]*H_perp[0]))*one_o_m_squared;

      			xyz[2]=	-(S[0]*H[1]-S[1]*H[0])
      						+ alpha_para[cell]*S[2]*S[2]*H_para[2]*one_o_m_squared
      						-alpha_perp[cell]*(S[0]*(S[2]*H_perp[0]-S[0]*H_perp[2])-S[1]*(S[1]*H_perp[2]-S[2]*H_perp[1]))*one_o_m_squared;

      			// Store dS in euler array
      			x_euler_array[cell]=xyz[0];
      			y_euler_array[cell]=xyz[1];
      			z_euler_array[cell]=xyz[2];

      			// Calculate Euler Step
      			x_cell_storage_array[cell]=S[0]+xyz[0]*dt[cell];
      			y_cell_storage_array[cell]=S[1]+xyz[1]*dt[cell];
      			z_cell_storage_array[cell]=S[2]+xyz[2]*dt[cell];

      		}
            for(unsigned int cell=0;cell<num_cells;cell++){
               x_cell_array[cell]=x_cell_storage_array[cell];
               y_cell_array[cell]=y_cell_storage_array[cell];
               z_cell_array[cell]=z_cell_storage_array[cell];
            }

            // Recalculate spin dependent fields
            fill (x_total_cell_field_array.begin(),x_total_cell_field_array.end(),0.0);
            fill (y_total_cell_field_array.begin(),y_total_cell_field_array.end(),0.0);
            fill (z_total_cell_field_array.begin(),z_total_cell_field_array.end(),0.0);


            for(unsigned int cell=0;cell<num_cells;cell++){
      			double m[3] = {x_cell_array[cell],y_cell_array[cell],z_cell_array[cell]};
      			double m_squared = m[1]*m[1]+m[2]*m[2]+m[0]*m[0];
      			double pf;
      			if(temperature<=mm::Tc[cell]){
      				pf = one_o_2_chi_para[cell]*(1.0 - m_squared/m_e_squared[cell]);
      			}
      			else{
      				pf = -2.0*one_o_2_chi_para[cell]*(1.0 + Tc_o_Tc_m_T[cell]*3.0*m_squared/5.0);
      			}

      			x_total_cell_field_array[cell] = (pf-one_o_chi_perp[cell])*m[0];
      			y_total_cell_field_array[cell] = (pf-one_o_chi_perp[cell])*m[1];
      			z_total_cell_field_array[cell] = (pf-0.0				 )*m[2];
      		}

            for(unsigned int cell=0;cell<num_cells;cell++){

               // Store local spin in Sand local field in H
               const double S[3] = {x_cell_array[cell],y_cell_array[cell],z_cell_array[cell]};
               const double H[3] = {x_total_cell_field_array[cell]+x_total_external_field_array[cell],
                                    y_total_cell_field_array[cell]+y_total_external_field_array[cell],
                                    z_total_cell_field_array[cell]+z_total_external_field_array[cell]};

               const double H_perp[3]={H[0]+Htx_perp[cell], H[1]+Hty_perp[cell], H[2]+Htz_perp[cell]};
               const double H_para[3]={H[0]+Htx_para[cell], H[1]+Hty_para[cell], H[2]+Htz_para[cell]};
               const double one_o_m_squared = 1.0/(S[1]*S[1]+S[2]*S[2]+S[0]*S[0]);

               // Calculate Delta S
               xyz[0]= 	-(S[1]*H[2]-S[2]*H[1])
                        + alpha_para[cell]*S[0]*S[0]*H_para[0]*one_o_m_squared
                        -alpha_perp[cell]*(S[1]*(S[0]*H_perp[1]-S[1]*H_perp[0])-S[2]*(S[2]*H_perp[0]-S[0]*H_perp[2]))*one_o_m_squared;

               xyz[1]= 	-(S[2]*H[0]-S[0]*H[2])
                        + alpha_para[cell]*S[1]*S[1]*H_para[1]*one_o_m_squared
                        -alpha_perp[cell]*(S[2]*(S[1]*H_perp[2]-S[2]*H_perp[1])-S[0]*(S[0]*H_perp[1]-S[1]*H_perp[0]))*one_o_m_squared;

               xyz[2]=	-(S[0]*H[1]-S[1]*H[0])
                        + alpha_para[cell]*S[2]*S[2]*H_para[2]*one_o_m_squared
                        -alpha_perp[cell]*(S[0]*(S[2]*H_perp[0]-S[0]*H_perp[2])-S[1]*(S[1]*H_perp[2]-S[2]*H_perp[1]))*one_o_m_squared;

               // Store dS in heun array
               x_heun_array[cell]=xyz[0];
               y_heun_array[cell]=xyz[1];
               z_heun_array[cell]=xyz[2];
            }
            for(unsigned int cell=0;cell<num_cells;cell++){
               x_cell_array[cell]=x_initial_cell_array[cell]+0.5*dt[cell]*(x_euler_array[cell]+x_heun_array[cell]);
               y_cell_array[cell]=y_initial_cell_array[cell]+0.5*dt[cell]*(y_euler_array[cell]+y_heun_array[cell]);
               z_cell_array[cell]=z_initial_cell_array[cell]+0.5*dt[cell]*(z_euler_array[cell]+z_heun_array[cell]);
            }

   }


	return EXIT_SUCCESS;
}
