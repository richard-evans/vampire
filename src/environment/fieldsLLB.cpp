

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

      std::vector<double> calculate_llb_fields(std::vector <double > m,
                                               double t,
                                               int cell,
                                               std::vector<double> x_array,
                                               std::vector<double> y_array,
                                               std::vector<double> z_array){



      double temperature;
      std::vector<double> spin_field(3,0.0);
      if (t == Tc) temperature = t+0.1;
      else temperature = t;
      const double one_o_2_chi_para = (one_o_chi_para/2.0);

      //the temperature is usually used as a reduced temperature.
      const double reduced_temperature = temperature/Tc;
      const double Tc_o_Tc_m_T = Tc/(temperature - Tc);

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


        const double m_e_squared = m_e*m_e;
        const double m_squared = m[1]*m[1]+m[2]*m[2]+m[0]*m[0];

        //calculates the intercell exchange field (pf) - this is dependent on temperature.
        double pf;
        if(temperature<=Tc) pf = one_o_2_chi_para*(1.0 - m_squared/m_e_squared);
        else pf = -2.0*one_o_2_chi_para*(1.0 + Tc_o_Tc_m_T*3.0*m_squared/5.0);

        //calculates the exchage fields as me^1.66 *A*(xi-xj)/m_e^2
        double exchange_field[3]={0.0,0.0,0.0};
        //is T < TC the exchange field = 0
        int cellj;
        double mj;
        double Ac = A*2/(m_e_squared*Ms*cell_size[0]*cell_size[1]);
        if (num_cells > 1){
         int start = neighbour_list_start_index[cell];
         int end =   neighbour_list_end_index[cell] +1;
         //loops over all other cells to sum the interaction
         for(int j = start;j<end;j++){

            cellj = neighbour_list_array[j];

          // std::cout << cell << '\t' << start << '\t' << end << '\t' << j << '\t' << cellj << std::endl;
            mj = sqrt(x_array[cellj]*x_array[cellj] + y_array[cellj]*y_array[cellj] + z_array[cellj]*z_array[cellj]);
            const double EC = pow(mj,1.66)*Ac;

            exchange_field[0] -= EC*(x_array[cellj] - x_array[cell]);
            exchange_field[1] -= EC*(y_array[cellj] - y_array[cell]);
            exchange_field[2] -= EC*(z_array[cellj] - z_array[cell]);
            //std::cout <<EC << std::endl;
            //std::cout << temperature << '\t' << cell << '\t' << start << '\t' << end << '\t' << j << '\t' << cellj << "\t" << A << '\t' << Ac << '\t' << exchange_field[0] << '\t' << x_array[cell] << '\t' << x_array[cellj]<<  std::endl;

          }
        }

        //Sum H = H_exch + H_A +H_exch_grains +H_App + H+dip
        spin_field[0] = pf*m[0] - one_o_chi_perp*m[0] + ext_field[0] + exchange_field[0] + dipole_field_x[cell];
        spin_field[1] = pf*m[1] - one_o_chi_perp*m[1] + ext_field[1] + exchange_field[1] + dipole_field_y[cell];
        spin_field[2] = pf*m[2]                       + ext_field[2] + exchange_field[2] + dipole_field_z[cell];
      // std::cout << m[0] << '\t' << m[1] << '\t' << m[2] << '\t' << spin_field[0] << '\t' << spin_field[1] << '\t' << spin_field[2] << '\t' << dipole_field_x[cell]<< '\t' << dipole_field_y[cell]<< '\t' << dipole_field_z[cell] <<std::endl;
        if (m[0] != m[0]) std::cin.get();
      return spin_field;

      }
    }
  }
