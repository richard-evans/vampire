
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

#include "random.hpp"
#include "sim.hpp"
#include "cells.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>


namespace mm = micromagnetic::internal;

namespace micromagnetic{

   void mm::step( int num_cells,
                  double temperature,
                  std::vector<double> x_array,
                  std::vector<double> y_array,
                  std::vector<double> z_array,
                  std::vector<double> ext_field,
                  double dt,
                  std::vector<double>& new_x_array,
                  std::vector<double>& new_y_array,
                  std::vector<double>& new_z_array){

      const double kB = 1.3806503e-23;
      double one_o_chi_perp, one_o_2_chi_para;
      double reduced_temperature, Tc_o_Tc_m_T;
      double m_e, alpha_para, alpha_perp, m_e_squared;

      std::vector<double> m(3,0.0);
      std::vector<double> spin_field(3,0.0);

      //6 arrays of gaussian random numbers to store the stochastic noise terms for x,y,z parallel and perperdicular
      std::vector <double> GW1x(num_cells);
      std::vector <double> GW1y(num_cells);
      std::vector <double> GW1z(num_cells);
      std::vector <double> GW2x(num_cells);
      std::vector <double> GW2y(num_cells);
      std::vector <double> GW2z(num_cells);

      //fill the noise terms
      generate (GW1x.begin(),GW1x.end(), mtrandom::gaussian);
      generate (GW1y.begin(),GW1y.end(), mtrandom::gaussian);
      generate (GW1z.begin(),GW1z.end(), mtrandom::gaussian);
      generate (GW2x.begin(),GW2x.end(), mtrandom::gaussian);
      generate (GW2y.begin(),GW2y.end(), mtrandom::gaussian);
      generate (GW2z.begin(),GW2z.end(), mtrandom::gaussian);


   //loops over all the cells to calculate the spin terms per cell - only filled cells where MS>0
   for (int cell =0; cell <num_cells; cell++){
      if (mm::ms[cell] > 1e-100){
         m[0] = x_array[cell];
         m[1] = y_array[cell];
         m[2] = z_array[cell];
         const double m_squared = m[1]*m[1]+m[2]*m[2]+m[0]*m[0];


         //chi is usually sued as 2/chi
         one_o_chi_perp = 1.0/mm::chi_perp[cell];
         one_o_2_chi_para = 1.0/(2.0*mm::chi_para[cell]);

         //the temperature is usually used as a reduced temperature.
         reduced_temperature = temperature/mm::Tc[cell];
         Tc_o_Tc_m_T = mm::Tc[cell]/(temperature - mm::Tc[cell]);

         //calculates me and alpha depending on the temperature.
         if (temperature<=mm::Tc[cell]){
            m_e = pow((mm::Tc[cell]-temperature)/(mm::Tc[cell]),0.365);
            alpha_para = (2.0/3.0)*mm::alpha[cell]*reduced_temperature;
            alpha_perp = mm::alpha[cell]*(1.0-temperature/(3.0*mm::Tc[cell]));
         }
         else{
            m_e = 0.0;
            alpha_para = mm::alpha[cell]*(2.0/3.0)*reduced_temperature;
            alpha_perp = alpha_para;
         }
         //   std::cout << temperature << '\t' << mm::chi_perp[cell] << '\t' << mm::chi_para[cell] << alpha_para << '\t' << alpha_perp << "\t" << m_e <<std::endl;
         m_e_squared = m_e*m_e;


         //calculates the final term of the field pf - this is depnedent on temperature.
         double pf;
         if(temperature<=mm::Tc[cell]) pf = one_o_2_chi_para*(1.0 - m_squared/m_e_squared);
         else pf = -2.0*one_o_2_chi_para*(1.0 + Tc_o_Tc_m_T*3.0*m_squared/5.0);


         //calculates the exchage fields as me^1.66 *A*(xi-xj)/m_e^2
         double exchange_field[3]={0.0,0.0,0.0};
         //is T < TC the exchange field = 0
         double sumx =0;
         double sumy = 0;
         double sumz = 0;
         if (num_cells > 1){
            int j2 = cell*num_cells;
            //loops over all other cells to sum the interaction
            double mi = pow(m_e_squared,0.5);

            for(int j = mm::macro_neighbour_list_start_index[cell];j<mm::macro_neighbour_list_end_index[cell] +1;j++){
               // calculate reduced exchange constant factor
               const int cellj = mm::macro_neighbour_list_array[j];
               const double mj = sqrt(x_array[cellj]*x_array[cellj] + y_array[cellj]*y_array[cellj] + z_array[cellj]*z_array[cellj]);
               const double A = mm::A[cellj]*pow(mj,1.66);
               exchange_field[0] -= A*(x_array[cellj] - x_array[cell]);
               exchange_field[1] -= A*(y_array[cellj] - y_array[cell]);
               exchange_field[2] -= A*(z_array[cellj] - z_array[cell]);
               //   std::cout << mm::macro_neighbour_list_start_index[cell] << '\t' << mm::macro_neighbour_list_end_index[cell] << '\t' << mm::macro_neighbour_list_array[j] <<std::endl;
            }
         }

         //Sum H = H_exch + H_A +H_exch_grains +H_App + H+dip
         spin_field[0] = pf*m[0] - one_o_chi_perp*m[0] + exchange_field[0] + mm::ext_field[0] + cells::x_field_array[cell]*cells::num_atoms_in_cell[cell];
         spin_field[1] = pf*m[1] - one_o_chi_perp*m[1] + exchange_field[1] + mm::ext_field[1] + cells::y_field_array[cell]*cells::num_atoms_in_cell[cell];
         spin_field[2] = pf*m[2]                       + exchange_field[2] + mm::ext_field[2] + cells::z_field_array[cell]*cells::num_atoms_in_cell[cell];

         //calculates the stochatic parallel and perpendicular terms
         double a;
         if (micromagnetic::stochastic == true) a = 1.0;
         else if (micromagnetic::stochastic == false) a = 0.0;
         double sigma_para = a*sqrt(2*kB*temperature*alpha_para/(mm::ms[cell]*dt)); //why 1e-27
         double sigma_perp = a*sqrt(2*kB*temperature*(alpha_perp-alpha_para)/(dt*mm::ms[cell]*alpha_perp*alpha_perp));

         const double H[3] = {spin_field[0], spin_field[1], spin_field[2]};

         //saves the noise terms to an array
         const double GW2t[3] = {GW2x[cell],GW2y[cell],GW2z[cell]};
         const double one_o_m_squared = 1.0/(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
         const double SdotH = m[0]*H[0] + m[1]*H[1] + m[2]*H[2];


         double xyz[3];
         //calculates the LLB equation
         xyz[0]= 	- (m[1]*H[2]-m[2]*H[1])
                  + alpha_para*m[0]*SdotH*one_o_m_squared
                  - alpha_perp*(m[1]*(m[0]*H[1]-m[1]*H[0])-m[2]*(m[2]*H[0]-m[0]*H[2]))*one_o_m_squared
                  + GW1x[cell]*sigma_para
                  - alpha_perp*(m[1]*(m[0]*GW2t[1]-m[1]*GW2t[0])-m[2]*(m[2]*GW2t[0]-m[0]*GW2t[2]))*one_o_m_squared*sigma_perp;

        xyz[1]= 	- (m[2]*H[0]-m[0]*H[2])
                  + alpha_para*m[1]*SdotH*one_o_m_squared
                  - alpha_perp*(m[2]*(m[1]*H[2]-m[2]*H[1])-m[0]*(m[0]*H[1]-m[1]*H[0]))*one_o_m_squared
                  + GW1y[cell]*sigma_para
                  - alpha_perp*(m[2]*(m[1]*GW2t[2]-m[2]*GW2t[1])-m[0]*(m[0]*GW2t[1]-m[1]*GW2t[0]))*one_o_m_squared*sigma_perp;

        xyz[2]=	- (m[0]*H[1]-m[1]*H[0])
                  + alpha_para*m[2]*SdotH*one_o_m_squared
                  - alpha_perp*(m[0]*(m[2]*H[0]-m[0]*H[2])-m[1]*(m[1]*H[2]-m[2]*H[1]))*one_o_m_squared
                  + GW1z[cell]*sigma_para
                  - alpha_perp*(m[0]*(m[2]*GW2t[0]-m[0]*GW2t[2])-m[1]*(m[1]*GW2t[2]-m[2]*GW2t[1]))*one_o_m_squared*sigma_perp;

                  //returns the values of the LLB as the new
                  new_x_array[cell] = xyz[0];
                  new_y_array[cell] = xyz[1];
                  new_z_array[cell] = xyz[2];
         }

      }
   }
}
