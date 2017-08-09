

// Vampire headers
#include "micromagnetic.hpp"
#include "cells.hpp"
#include "internal.hpp"

// micromagnetic module headers
#include <stdlib.h>
#include <vector>
#include "errors.hpp"
#include "vio.hpp"
#include "random.hpp"

namespace micromagnetic
{

   namespace internal
   {

      std::vector<double> calculate_llg_fields(std::vector <double > m,
                                               double temperature,
                                               int num_cells,
                                               int cell,
                                               std::vector<double> x_array,
                                               std::vector<double> y_array,
                                               std::vector<double> z_array){


      std::vector<double> spin_field(3,0.0);

      const double kB = 1.3806503e-23;

      //chi is usually used as 2/chi
      const double one_o_2_chi_para = (one_o_chi_para[cell]/2.0);

      //the temperature is usually used as a reduced temperature.
      const double reduced_temperature = temperature/Tc[cell];
      const double Tc_o_Tc_m_T = Tc[cell]/(temperature - Tc[cell]);

      //calcualted m_e and alpha temperature dependant
      if (temperature<=Tc[cell]){
         m_e[cell] = pow((Tc[cell]-temperature)/(Tc[cell]),0.365);
         alpha_para[cell] = (2.0/3.0)*alpha[cell]*reduced_temperature;
         alpha_perp[cell] = alpha[cell]*(1.0-temperature/(3.0*Tc[cell]));
      }
      else{
         m_e[cell] = 0.01;
         alpha_para[cell] = alpha[cell]*(2.0/3.0)*reduced_temperature;
         alpha_perp[cell] = alpha_para[cell];
      }

      //saved me_2
      const double m_e_squared = m_e[cell]*m_e[cell];

      //calculates the intercell exchange field (pf) - this is dependent on temperature.
      double pf = 0;
      if(temperature<=Tc[cell]) pf = 0;
      else pf = -2.0*one_o_2_chi_para*(1.0 + Tc_o_Tc_m_T*3.0*m_e_squared/5.0);

      //calculates the exchage fields as me^1.66 *A*(xi-xj)/m_e^2
      double exchange_field[3]={0.0,0.0,0.0};
      //is T < TC the exchange field = 0
      if (num_cells > 1){

        const int start = macro_neighbour_list_start_index[cell];
        const int end = macro_neighbour_list_end_index[cell] +1;
        for(int j = start;j< end;j++){
          // calculate reduced exchange constant factor
          const int cellj = macro_neighbour_list_array[j];
          const double mj = m_e[cellj];
          const double Ac = A[cellj]*pow(mj,1.66);
          exchange_field[0] -= Ac*(x_array[cellj]*m_e[cellj] - x_array[cell]*m_e[cell]);
          exchange_field[1] -= Ac*(y_array[cellj]*m_e[cellj] - y_array[cell]*m_e[cell]);
          exchange_field[2] -= Ac*(z_array[cellj]*m_e[cellj] - z_array[cell]*m_e[cell]);
        }
      }
      //calcualtes thesigma values
      double sigma_para = sqrt(2*kB*temperature*alpha_para[cell]/(ms[cell]*mp::dt));
      double sigma_perp = sqrt(2*kB*temperature*(alpha_perp[cell]-alpha_para[cell])/(mp::dt*ms[cell]*alpha_perp[cell]*alpha_perp[cell]));



      //Sum H = H_exch + H_A +H_exch_grains +H_App + H+dip
      spin_field[0] = one_o_chi_perp[cell]*m[0]*m_e[cell] + ext_field[0] + cells::field_array_x[cell] + exchange_field[0] + sigma_para*mtrandom::gaussian();
      spin_field[1] = one_o_chi_perp[cell]*m[1]*m_e[cell] + ext_field[1] + cells::field_array_y[cell] + exchange_field[1] + sigma_para*mtrandom::gaussian();
      spin_field[2] =                                     + ext_field[2] + cells::field_array_z[cell] + exchange_field[2] + sigma_para*mtrandom::gaussian();

      return spin_field;
     }
   }
 }
