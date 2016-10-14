
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

#include "random.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>

namespace mm = micromagnetic::internal;


int LLB_serial_heun(int num_steps,int num_cells,double temperature,std::vector<double>& x_mag_array,std::vector<double>& y_mag_array,std::vector<double>& z_mag_array,double Hx,double Hy,double Hz, double H,double dt,std::vector <double> volume_array);

namespace micromagnetic{
/// Master LLB Function - dispatches code path to desired LLB routine
/// \f$ \frac{\partial S}{\partial t} \f$
int LLB(int num_steps,int num_cells,double temperature,std::vector<double>& x_mag_array,std::vector<double>& y_mag_array,std::vector<double>& z_mag_array,double Hx,double Hy,double Hz, double H,double dt, std::vector <double> volume_array){

   //----------------------------------------------------------
   // check calling of routine if error checking is activated
   //----------------------------------------------------------


   //if(err::check==true){std::cout << "LLB has been called" << std::endl;}

   #ifdef MPICF
      //LLB_mpi(num_steps);
   #else
   LLB_serial_heun(num_steps,num_cells,temperature,x_mag_array,y_mag_array,z_mag_array,Hx,Hy,Hz,H,dt, volume_array);
   #endif

   return 0;
}
}

int LLB_serial_heun(
   int num_steps,
   int num_cells,
   double temperature,
   std::vector<double>& x_mag_array,
   std::vector<double>& y_mag_array,
   std::vector<double>& z_mag_array,
   double Hx,
   double Hy,
   double Hz,
   double H,
   double dt,
   std::vector <double> volume_array
){

	const double kB = 1.3806503e-23;
   double one_o_chi_perp, one_o_2_chi_para, reduced_temperature, Tc_o_Tc_m_T, m_e, alpha_para, alpha_perp, m_e_squared;
   std::vector<double> m(3,0.0);
   std::vector<double> spin_field(3,0.0);
   std::vector<double> ext_field(3,0.0);
   std::vector<double> dip_field(3,0.0);

   std::vector<double> x_array(num_cells,0.0);
   std::vector<double> y_array(num_cells,0.0);
   std::vector<double> z_array(num_cells,0.0);

   std::vector<double> x_euler_array(num_cells,0.0);
   std::vector<double> y_euler_array(num_cells,0.0);
   std::vector<double> z_euler_array(num_cells,0.0);

   std::vector<double> x_heun_array(num_cells,0.0);
   std::vector<double> y_heun_array(num_cells,0.0);
   std::vector<double> z_heun_array(num_cells,0.0);

   std::vector<double> mx_store(num_cells,0.0);
   std::vector<double> my_store(num_cells,0.0);
   std::vector<double> mz_store(num_cells,0.0);

   std::vector<double> mx_init(num_cells,0.0);
   std::vector<double> my_init(num_cells,0.0);
   std::vector<double> mz_init(num_cells,0.0);


   std::vector <double> GW1x(num_cells);
   std::vector <double> GW1y(num_cells);
   std::vector <double> GW1z(num_cells);
   std::vector <double> GW2x(num_cells);
   std::vector <double> GW2y(num_cells);
   std::vector <double> GW2z(num_cells);

   //--------------------------------------------------------------------------------------------------------------------------------

   //----------------------------------------------------Loop over N steps-------------- --------------------------------------------

   //--------------------------------------------------------------------------------------------------------------------------------

   ext_field[0] = H*Hx;
   ext_field[1] = H*Hy;
   ext_field[2] = H*Hz;

   for(int t=0;t<num_steps;t++){

      for (int cell = 0; cell < num_cells; cell ++)
      {
         x_array[cell] = x_mag_array[cell]/mm::ms[cell];
         y_array[cell] = x_mag_array[cell]/mm::ms[cell];
         z_array[cell] = x_mag_array[cell]/mm::ms[cell];
      }

      for(unsigned int cell=0;cell<num_cells;cell++){
         mx_init[cell] = x_array[cell];
         my_init[cell] = y_array[cell];
         mz_init[cell] = z_array[cell];
      }

      mm::chi_para =  mm::calculate_chi_para(num_cells, temperature);
      mm::chi_perp =  mm::calculate_chi_perp(num_cells, temperature);

      generate (GW1x.begin(),GW1x.end(), mtrandom::gaussian);
      generate (GW1y.begin(),GW1y.end(), mtrandom::gaussian);
      generate (GW1z.begin(),GW1z.end(), mtrandom::gaussian);
      generate (GW2x.begin(),GW2x.end(), mtrandom::gaussian);
      generate (GW2y.begin(),GW2y.end(), mtrandom::gaussian);
      generate (GW2z.begin(),GW2z.end(), mtrandom::gaussian);

      //--------------------------------------------------------------------------------------------------------------------------------

      //----------------------------------------------------Calculation of the euler vectors --------------------------------------------

      //--------------------------------------------------------------------------------------------------------------------------------


      for (int cell =0; cell <num_cells; cell++)
      {
         one_o_chi_perp = 1.0/mm::chi_perp[cell];
         one_o_2_chi_para = 1.0/(2.0*mm::chi_para[cell]);
         reduced_temperature = temperature/mm::Tc[cell];
         Tc_o_Tc_m_T = mm::Tc[cell]/(temperature - mm::Tc[cell]);

         if (temperature<=mm::Tc[cell])
         {
            m_e = pow((mm::Tc[cell]-temperature)/(mm::Tc[cell]),0.365);
            alpha_para = (2.0/3.0)*mm::alpha[cell]*reduced_temperature;
            alpha_perp = mm::alpha[cell]*(1.0-temperature/(3.0*mm::Tc[cell]));
         }
         else
         {
            m_e = 0.0;
            alpha_para = mm::alpha[cell]*(2.0/3.0)*reduced_temperature;
            alpha_perp = alpha_para;
         }

         m_e_squared = m_e*m_e;
         m[0] = x_array[cell];
         m[1] = y_array[cell];
         m[2] = z_array[cell];

         const double m_squared = m[1]*m[1]+m[2]*m[2]+m[0]*m[0];

         double pf;
         if(temperature<=mm::Tc[cell]) pf = one_o_2_chi_para*(1.0 - m_squared/m_e_squared);
         else pf = -2.0*one_o_2_chi_para*(1.0 + Tc_o_Tc_m_T*3.0*m_squared/5.0);



         double exchange_field[3]={0.0,0.0,0.0};
         int j2 = cell*num_cells;
         double mcell = pow(m[0]*m[0] +m[1]*m[1] + m[2]*m[2],0.5);
         double mj;
			for(int j=0;j<num_cells;j++){
            mj = pow(x_array[j]*x_array[j] +y_array[j]*y_array[j] + z_array[j]*z_array[j],0.5);
				exchange_field[0]+=pow(m_e,1.66)*(x_array[j]/mj -m[0]/mcell)*mm::Ax[j2];
				exchange_field[1]+=pow(m_e,1.66)*(y_array[j]/mj -m[1]/mcell)*mm::Ay[j2];
				exchange_field[2]+=pow(m_e,1.66)*(z_array[j]/mj -m[2]/mcell)*mm::Az[j2];
            j2++;
   		}
         spin_field[0] = (pf)*m[0];// - (one_o_chi_perp)*m[0];// + exchange_field[0];
         spin_field[1] = (pf)*m[1];// - (one_o_chi_perp)*m[1];// + exchange_field[1];
         spin_field[2] = (pf - 0             )*m[2];// + exchange_field[2];

         double g = sqrt(mm::gamma[cell]*mm::gamma[cell]);
         double sigma_para = sqrt(2*kB*temperature*alpha_para*g/mm::ms[cell])*dt; //why 1e-27
         double sigma_perp = sqrt(2*kB*temperature*(alpha_perp-alpha_para)/(g*mm::ms[cell]*alpha_perp*alpha_perp))*dt;

         const double H[3] = {spin_field[0],// + ext_field[0] + dip_field[0],
                              spin_field[1],// + ext_field[1] + dip_field[1],
                              spin_field[2]};// + ext_field[2] + dip_field[2]};

         const double GW2t[3] = {GW2x[cell],GW2y[cell],GW2z[cell]};
         const double one_o_m_squared = 1.0/(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
         const double SdotH = m[0]*H[0] + m[1]*H[1] + m[2]*H[2];

         double xyz[3];

         xyz[0]= 	- mm::gamma[cell]*(m[1]*H[2]-m[2]*H[1])
                  + g*alpha_para*m[0]*SdotH*one_o_m_squared
                  - g*alpha_perp*(m[1]*(m[0]*H[1]-m[1]*H[0])-m[2]*(m[2]*H[0]-m[0]*H[2]))*one_o_m_squared
                  + GW1x[cell]*sigma_para
                  - g*alpha_perp*(m[1]*(m[0]*GW2t[1]-m[1]*GW2t[0])-m[2]*(m[2]*GW2t[0]-m[0]*GW2t[2]))*one_o_m_squared*GW2x[cell]*sigma_perp;

         xyz[1]= 	- mm::gamma[cell]*(m[2]*H[0]-m[0]*H[2])
                  + g*alpha_para*m[1]*SdotH*one_o_m_squared
                  - g*alpha_perp*(m[2]*(m[1]*H[2]-m[2]*H[1])-m[0]*(m[0]*H[1]-m[1]*H[0]))*one_o_m_squared
                  + GW1y[cell]*sigma_para
                  - g*alpha_perp*(m[2]*(m[1]*GW2t[2]-m[2]*GW2t[1])-m[0]*(m[0]*GW2t[1]-m[1]*GW2t[0]))*one_o_m_squared*GW2y[cell]*sigma_perp;

         xyz[2]=	- mm::gamma[cell]*(m[0]*H[1]-m[1]*H[0])
                  + g*alpha_para*m[2]*SdotH*one_o_m_squared
                  - g*alpha_perp*(m[0]*(m[2]*H[0]-m[0]*H[2])-m[1]*(m[1]*H[2]-m[2]*H[1]))*one_o_m_squared
                  + GW1z[cell]*sigma_para
                  - g*alpha_perp*(m[0]*(m[2]*GW2t[0]-m[0]*GW2t[2])-m[1]*(m[1]*GW2t[2]-m[2]*GW2t[1]))*one_o_m_squared*GW2z[cell]*sigma_perp;


         x_euler_array[cell]=xyz[0];
         y_euler_array[cell]=xyz[1];
         z_euler_array[cell]=xyz[2];

         x_array[cell] = m[0]+xyz[0]*dt;
         y_array[cell] = m[1]+xyz[1]*dt;
         z_array[cell] = m[2]+xyz[2]*dt;

         x_mag_array[cell] = x_array[cell]*mm::ms[cell];
         y_mag_array[cell] = y_array[cell]*mm::ms[cell];
         z_mag_array[cell] = z_array[cell]*mm::ms[cell];

std::cout << temperature << '\t' << x_array[cell] << '\t' << y_array[cell] << "\t" << z_array[cell] << "\t" << m_squared << '\t' << pf << "\t" << one_o_chi_perp << "\t" << xyz[0] <<"\t"<< exchange_field[0] <<std::endl;
      }

   }
}
