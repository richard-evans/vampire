

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
		//std::cout << "f3" <<std::endl;
      const double m_e_squared = m_e[cell]*m_e[cell];

      //calculates the intercell exchange field (pf) - this is dependent on temperature.
      double pf = 0;
      if(temperature<=Tc[cell]) pf = 0;
      else pf = -2.0*one_o_2_chi_para*(1.0 + Tc_o_Tc_m_T*3.0*m_e_squared/5.0);
		//std::cout << "f4" <<std::endl;

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
          mj = m_e[cell];
          const double Ac = A[cellj]*pow(mj,1.66);
          exchange_field[0] -= Ac*(x_array[cellj] - x_array[cell]);
          exchange_field[1] -= Ac*(y_array[cellj] - y_array[cell]);
          exchange_field[2] -= Ac*(z_array[cellj] - z_array[cell]);
        }
      }
      double sigma_para = sqrt(2*kB*temperature*alpha_para[cell]/(ms[cell]*mp::dt));
      double sigma_perp = sqrt(2*kB*temperature*(alpha_perp[cell]-alpha_para[cell])/(mp::dt*ms[cell]*alpha_perp[cell]*alpha_perp[cell]));



      //Sum H = H_exch + H_A +H_exch_grains +H_App + H+dip
      spin_field[0] = pf*m[0] - one_o_chi_perp[cell]*m[0] + ext_field[0] + cells::field_array_x[cell]*cells::num_atoms_in_cell[cell] + exchange_field[0] + sigma_perp*mtrandom::gaussian();
      spin_field[1] = pf*m[1] - one_o_chi_perp[cell]*m[1] + ext_field[1] + cells::field_array_y[cell]*cells::num_atoms_in_cell[cell] + exchange_field[1] + sigma_perp*mtrandom::gaussian();
      spin_field[2] = pf*m[2]                       + ext_field[2] + cells::field_array_z[cell]*cells::num_atoms_in_cell[cell] + exchange_field[2] + sigma_para*mtrandom::gaussian();

    //  std::cout <<temperature << '\t' << pf << "\t" << spin_field[0] << '\t' << std::endl;

      return spin_field;
     }
   }
 }
