//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Vampire headers
#include "micromagnetic.hpp"
#include "cells.hpp"
#include "internal.hpp"
#include "../cells/internal.hpp"
// micromagnetic module headers
#include <stdlib.h>
#include <vector>
#include "errors.hpp"
#include "cells.hpp"
#include "vio.hpp"
#include "environment.hpp"

namespace micromagnetic{

   namespace internal{

      std::vector<double> calculate_llb_fields(std::vector <double > m,
                                               double temperature,
                                               int num_cells,
                                               int cell,
                                               std::vector<double> x_array,
                                               std::vector<double> y_array,
                                               std::vector<double> z_array){


      //array to save the field to for each cell
      std::vector<double> spin_field(3,0.0);

      //chi is usually used as 1/2*chi
      const double one_o_2_chi_para = (one_o_chi_para[cell]/2.0);


      //the temperature is usually used as a reduced temperature or Tc/(Tc-T).
      const double reduced_temperature = temperature/Tc[cell];
      const double Tc_o_Tc_m_T = Tc[cell]/(temperature - Tc[cell]);

      //sets m_e and alpha temperature dependant parameters
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

      //m and me are usually used squared
      const double m_e_squared = m_e[cell]*m_e[cell];
      const double m_squared = m[1]*m[1]+m[2]*m[2]+m[0]*m[0];

      //calculates the intercell exchange field (pf) - this is dependent on temperature.
      double pf = 0;
      if(temperature<=Tc[cell]) pf = one_o_2_chi_para*(1.0 - m_squared/m_e_squared);
      else pf = -2.0*one_o_2_chi_para*(1.0 + Tc_o_Tc_m_T*3.0*m_squared/5.0);

      //calculates the exchage fields as me^1.66 *A*(xi-xj)/m_e^2
      //array to store the exchanege field

      double exchange_field[3]={0.0,0.0,0.0};
      if (num_cells > 1){

         //loops over all other cells with interactions to this cell
         const int start = macro_neighbour_list_start_index[cell];
         const int end = macro_neighbour_list_end_index[cell] +1;

         for(int j = start;j< end;j++){
            // calculate reduced exchange constant factor
            const int cellj = macro_neighbour_list_array[j];
            const double mj = sqrt(x_array[cellj]*x_array[cellj] + y_array[cellj]*y_array[cellj] + z_array[cellj]*z_array[cellj]);
            double Ac = A[j]*pow(mj,1.66);

      //if (zj < 140 && zi > 150)   std::cout << cell << '\t' << cellj << '\t' << Ac << '\t' << x_array[cellj] - x_array[cell] << "\t" << exchange_field[0] <<"\t" << exchange_field[1] <<"\t" << exchange_field[2] << std::endl;
      //else std::cout << "A" << cell << '\t' << cellj << '\t' << Ac << '\t' << x_array[cellj] - x_array[cell] << "\t" << exchange_field[0] <<"\t" << exchange_field[1] <<"\t" << exchange_field[2] << std::endl;
            exchange_field[0] -= Ac*(x_array[cellj] - x_array[cell]);
            exchange_field[1] -= Ac*(y_array[cellj] - y_array[cell]);
            exchange_field[2] -= Ac*(z_array[cellj] - z_array[cell]);


         }
      }
   //   std::cin.get();

//if (exchange_field[0] != exchange_field[0]) std::cin.get();

      //Sum H = H_exch + H_A +H_exch_grains +H_App + H+dip
      spin_field[0] = pf*m[0] + exchange_field[0] + pinning_field_x[cell] + ext_field[0] - one_o_chi_perp[cell]*m[0] + cells::field_array_x[cell];// + ext_field[0] + exchange_field[0];// + pinning_field_x[cell];// + cells::field_array_x[cell];
      spin_field[1] = pf*m[1] + exchange_field[1] + pinning_field_y[cell] + ext_field[1] - one_o_chi_perp[cell]*m[1] + cells::field_array_y[cell];// + ext_field[1] + exchange_field[1];// + pinning_field_y[cell];// + cells::field_array_y[cell];
      spin_field[2] = pf*m[2] + exchange_field[2] + pinning_field_z[cell] + ext_field[2]                             + cells::field_array_z[cell];// + ext_field[2] + exchange_field[2];// + pinning_field_z[cell];// + cells::field_array_z[cell];

      //if environment is enabled add the environment field.


      if (environment::enabled){
         spin_field[0] = spin_field[0] + environment::environment_field_x[cell];
         spin_field[1] = spin_field[1] + environment::environment_field_y[cell];
         spin_field[2] = spin_field[2] + environment::environment_field_z[cell];
      }

      //std::cout << environment::environment_field_x[cell] << '\t' << environment::environment_field_x[cell]  << '\t' << environment::environment_field_x[cell] << std::endl;
      return spin_field;
   }

} // end of internal namespace

} // end of micromagnetic namespace
