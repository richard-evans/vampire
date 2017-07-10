

// Vampire headers
#include "micromagnetic.hpp"
#include "cells.hpp"
#include "internal.hpp"

// micromagnetic module headers
#include <stdlib.h>
#include <vector>
#include "errors.hpp"
#include "cells.hpp"
#include "vio.hpp"
#include "environment.hpp"
namespace micromagnetic
{

   namespace internal
   {

      std::vector<double> calculate_llb_fields(std::vector <double > m,
                                               double temperature,
                                               int num_cells,
                                               int cell,
                                               std::vector<double> x_array,
                                               std::vector<double> y_array,
                                               std::vector<double> z_array){

                                                // std::cout <<'a' << std::endl;
      std::vector<double> spin_field(3,0.0);

      //chi is usually used as 1/2*chi
      const double one_o_2_chi_para = (one_o_chi_para[cell]/2.0);

      //the temperature is usually used as a reduced temperature or Tc/(Tc-T).
      const double reduced_temperature = temperature/Tc[cell];
      const double Tc_o_Tc_m_T = Tc[cell]/(temperature - Tc[cell]);

      if (temperature<=Tc[cell]){
         m_e[cell] = pow((Tc[cell]-temperature)/(Tc[cell]),0.365);
         alpha_para[cell] = (2.0/3.0)*alpha[cell]*reduced_temperature;
         alpha_perp[cell] = alpha[cell]*(1.0-temperature/(3.0*Tc[cell]));
      }
      else{
         m_e[cell] = 0.0;
         alpha_para[cell] = alpha[cell]*(2.0/3.0)*reduced_temperature;
         alpha_perp[cell] = alpha_para[cell];
      }

      const double m_e_squared = m_e[cell]*m_e[cell];
      const double m_squared = m[1]*m[1]+m[2]*m[2]+m[0]*m[0];

      //calculates the intercell exchange field (pf) - this is dependent on temperature.
      double pf = 0;
      if(temperature<=Tc[cell]) pf = one_o_2_chi_para*(1.0 - m_squared/m_e_squared);
      else pf = -2.0*one_o_2_chi_para*(1.0 + Tc_o_Tc_m_T*3.0*m_squared/5.0);


      //calculates the exchage fields as me^1.66 *A*(xi-xj)/m_e^2
      double exchange_field[3]={0.0,0.0,0.0};
      //is T < TC the exchange field = 0
      int cellj;
      double mj;


      if (num_cells > 1){
       int j2 = cell*num_cells;
       //loops over all other cells to sum the interaction
       for(int j = macro_neighbour_list_start_index[cell];j<macro_neighbour_list_end_index[cell] +1;j++){
          // calculate reduced exchange constant factor
          cellj = macro_neighbour_list_array[j];
          mj = sqrt(x_array[cellj]*x_array[cellj] + y_array[cellj]*y_array[cellj] + z_array[cellj]*z_array[cellj]);
          const double Ac = A[cellj]*pow(mj,1.66);
          exchange_field[0] -= Ac*(x_array[cellj] - x_array[cell]);
          exchange_field[1] -= Ac*(y_array[cellj] - y_array[cell]);
          exchange_field[2] -= Ac*(z_array[cellj] - z_array[cell]);
        }
      }

      //Sum H = H_exch + H_A +H_exch_grains +H_App + H+dip
      spin_field[0] = pf*m[0] - one_o_chi_perp[cell]*m[0] + ext_field[0] + cells::field_array_x[cell] + exchange_field[0];
      spin_field[1] = pf*m[1] - one_o_chi_perp[cell]*m[1] + ext_field[1] + cells::field_array_y[cell] + exchange_field[1];
      spin_field[2] = pf*m[2]                             + ext_field[2] + cells::field_array_z[cell] + exchange_field[2];

     if (environment::enabled){
       spin_field[0] = spin_field[0] + environment::environment_field_x[cell];
       spin_field[1] = spin_field[1] + environment::environment_field_y[cell];
       spin_field[2] = spin_field[2] + environment::environment_field_z[cell];
      // std::cout << cell << '\t' << environment::environment_field_x[cell]  << '\t' << environment::environment_field_y[cell]  << '\t' << environment::environment_field_z[cell] << "\t" << spin_field[0] << '\t' << spin_field[1] << '\t' << spin_field[2] << std::endl;
     }
    //  std::cout << spin_field[0] << '\t' << environment::environment_field_x[cell] <<std::endl;
      return spin_field;
     }
   }
 }
