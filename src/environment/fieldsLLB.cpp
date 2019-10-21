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
#include "environment.hpp"
#include "internal.hpp"

// environment module headers
#include <stdlib.h>
#include <vector>
#include "errors.hpp"
#include "vio.hpp"
#include "sim.hpp"
namespace environment
{

   namespace internal
   {

      std::vector<double> calculate_llb_fields(std::vector <double>& m,
         double t,
         int cell,
         std::vector<double>& x_array,
         std::vector<double>& y_array,
         std::vector<double>& z_array){

            double temperature;
            //vector to store fields
            std::vector<double> spin_field(3,0.0);
            //stops T ==Tc - breaks LLB
            if (t == Tc) temperature = t+0.1;
            else temperature = t;

            const double one_o_2_chi_para = (one_o_chi_para/2.0);

            //the temperature is usually used as a reduced temperature.
            const double reduced_temperature = temperature/Tc;
            const double Tc_o_Tc_m_T = Tc/(temperature - Tc);

            //calcualtes m_e, and alpha temeprature dependant
            if (temperature<=Tc){
               m_e = pow((Tc-temperature)/(Tc),0.365);
               alpha_para = (2.0/3.0)*alpha*reduced_temperature;
               alpha_perp = alpha*(1.0-temperature/(3.0*Tc));
            }
            else{
               m_e = 0.0;
               alpha_para = alpha*(2.0/3.0)*reduced_temperature;
               alpha_perp = alpha_para;
            }
            //m and m_e are squared
            const double m_e_squared = m_e*m_e;
            const double m_squared = m[1]*m[1]+m[2]*m[2]+m[0]*m[0];

            //calculates the intercell exchange field (pf) - this is dependent on temperature.
            double pf;
            if(temperature<=Tc) pf = one_o_2_chi_para*(1.0 - m_squared/m_e_squared);
            else pf = -2.0*one_o_2_chi_para*(1.0 + Tc_o_Tc_m_T*3.0*m_squared/5.0);

            //calculates the exchage fields as me^1.66 *A*(xi-xj)/m_e^2
            double exchange_field[3]={0.0,0.0,0.0};
            const double Ar = A*2.0/(Ms);
            //saves the x,y,z prefactors for the exchange constant
            //const double Acx = A*2/(Ms)*cell_size[2]*cell_size[1];
            //const double Acy = A*2/(Ms)*cell_size[2]*cell_size[0];
            //const double Acz = A*2/(Ms)*cell_size[0]*cell_size[1];
            const double Acx = Ar*cell_size[2]*cell_size[1];
            const double Acy = Ar*cell_size[2]*cell_size[0];
            const double Acz = Ar*cell_size[0]*cell_size[1];

            if (num_cells > 1){
               const int start = neighbour_list_start_index[cell];
               const int end =   neighbour_list_end_index[cell] +1;
               //loops over cells with an interaction from the neighbour lists
               for(int j = start;j<end;j++){

                  const int cellj = neighbour_list_array[j];
                  //calculate |mj|
                  const double mj = sqrt(x_array[cellj]*x_array[cellj] + y_array[cellj]*y_array[cellj] + z_array[cellj]*z_array[cellj]);
                  //calcaultes the temperature dependant terms
                  const double AT = pow(mj,1.66);
                  const double ECx = AT*Acx;
                  const double ECy = AT*Acy;
                  const double ECz = AT*Acz;
                  //calcualtes the exchange field from i to j
                  exchange_field[0] -= ECx*(x_array[cellj] - x_array[cell]);
                  exchange_field[1] -= ECy*(y_array[cellj] - y_array[cell]);
                  exchange_field[2] -= ECz*(z_array[cellj] - z_array[cell]);

               }
            }


            //Sum H = H_exch + H_A +H_exch_grains +H_App + H+dip
            spin_field[0] = pf*m[0] + exchange_field[0] - one_o_chi_perp*m[0] + ext_field[0] + dipole_field_x[cell] + env_field_uv[0] + bias_field_x[cell];
            spin_field[1] = pf*m[1] + exchange_field[1] - one_o_chi_perp*m[1] + ext_field[1] + dipole_field_y[cell] + env_field_uv[1] + bias_field_y[cell];
            spin_field[2] = pf*m[2] + exchange_field[2]                       + ext_field[2] + dipole_field_z[cell] + env_field_uv[2] + bias_field_z[cell];
          //    std::cout << "x" << '\t' << pf  << '\t' << m[0] << '\t' << exchange_field[0] << "\t" << dipole_field_x[cell] << "\t" << spin_field[0] <<std::endl;
          //    std::cout << "y" << '\t' << pf  << '\t' << m[1] << '\t' << exchange_field[1] << "\t" << dipole_field_y[cell] << "\t" << spin_field[1] <<std::endl;
          //    std::cout << "z" << '\t' << pf  << '\t' << m[2] << '\t' << exchange_field[2] << "\t" << dipole_field_z[cell] << "\t" << spin_field[2] <<std::endl;
              if (spin_field[0] != spin_field[0] ){
                std::cin.get();
                std::cout << "Error in environment"<<std::endl;
            }
            return spin_field;

         }
      }
   }
