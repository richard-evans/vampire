
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include "cells.hpp"
#include "atoms.hpp"
#include "vio.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>


namespace mm = micromagnetic::internal;


int LLB_serial_heun(int num_steps,int num_cells,double temperature,std::vector<double>& x_mag_array,std::vector<double>& y_mag_array,std::vector<double>& z_mag_array,double Hx,double Hy,double Hz, double H,double dt,std::vector <double> volume_array, int N);

namespace micromagnetic{
/// Master LLB Function - dispatches code path to desired LLB routine
/// \f$ \frac{\partial S}{\partial t} \f$
int LLB(int num_steps,int num_cells,double temperature,std::vector<double>& x_mag_array,std::vector<double>& y_mag_array,std::vector<double>& z_mag_array,double Hx,double Hy,double Hz, double H,double dt, std::vector <double> volume_array, int N){

   //----------------------------------------------------------
   // check calling of routine if error checking is activated
   //----------------------------------------------------------


   //if(err::check==true){std::cout << "LLB has been called" << std::endl;}

   #ifdef MPICF
      //LLB_mpi(num_steps);
   #else
   LLB_serial_heun(num_steps,num_cells,temperature,x_mag_array,y_mag_array,z_mag_array,Hx,Hy,Hz,H,dt, volume_array, N);
   #endif

   return 0;
}
}

int LLB_serial_heun( int num_steps,
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
                     std::vector <double> volume_array,
                     int N
                  ){




   const double kB = 1.3806503e-23;
   //The external fields equal the length of the field times the applied field vector.
   //This is saved to an array.
   mm::ext_field[0] = H*Hx;
   mm::ext_field[1] = H*Hy;
   mm::ext_field[2] = H*Hz;


   //loop over all cells and calculate m = M/Ms and save it to x,y,z_array
   //there are lots of blank cells so only the cells where ms!=0 are calcualted.
   //save this new m as the initial value, so it can be saved and used in the final equation.
   for (int cell = 0; cell < num_cells; cell++){
         mm::x_array[cell] = x_mag_array[cell]/mm::ms[cell];
         mm::y_array[cell] = y_mag_array[cell]/mm::ms[cell];
         mm::z_array[cell] = z_mag_array[cell]/mm::ms[cell];
         mm::mx_init[cell] = mm::x_array[cell];
         mm::my_init[cell] = mm::y_array[cell];
         mm::mz_init[cell] = mm::z_array[cell];
   }

   //calculte chi(T).
   mm::chi_para =  mm::calculate_chi_para(num_cells, temperature);
   mm::chi_perp =  mm::calculate_chi_perp(num_cells, temperature);

   //calls a function (step) which calculates 1 intergration step - the euler step and saves it to an euler array
   mm::step(num_cells, temperature, mm::x_array,mm::y_array,mm::z_array, mm::ext_field, dt, mm::x_euler_array, mm::y_euler_array, mm::z_euler_array);

   //these new x postiion are stored in an array (store)
   //x = x+step*dt
   for (int cell = 0; cell < num_cells; cell++){
      if (mm::ms[cell] > 0){
         mm::mx_store[cell] = mm::x_array[cell] + mm::x_euler_array[cell]*dt;
         mm::my_store[cell] = mm::y_array[cell] + mm::y_euler_array[cell]*dt;
         mm::mz_store[cell] = mm::z_array[cell] + mm::z_euler_array[cell]*dt;
      }
   }

   //calls a function (step) which calculates a heun step based on the position after the euler step and save it to a heun array
   mm::step(num_cells, temperature, mm::mx_store, mm::my_store, mm::mz_store, mm::ext_field, dt, mm::x_heun_array, mm::y_heun_array, mm::z_heun_array);


   //calcualtes the final position as x = xinital + 1/2(euler+heun)*dt
   for (int cell = 0; cell < num_cells; cell++){
      if (mm::ms[cell] > 0){
      mm::x_array[cell] = mm::mx_init[cell] + 0.5*dt*(mm::x_euler_array[cell] + mm::x_heun_array[cell]);
      mm::y_array[cell] = mm::my_init[cell] + 0.5*dt*(mm::y_euler_array[cell] + mm::y_heun_array[cell]);
      mm::z_array[cell] = mm::mz_init[cell] + 0.5*dt*(mm::z_euler_array[cell] + mm::z_heun_array[cell]);
      }
   }

   for (int cell = 0; cell < num_cells; cell ++)
   {
         if (mm::ms[cell] > 1e-100){
         cells::x_mag_array[cell] = mm::x_array[cell]*mm::ms[cell];
         cells::y_mag_array[cell] = mm::y_array[cell]*mm::ms[cell];
         cells::z_mag_array[cell] = mm::z_array[cell]*mm::ms[cell];
      }
   }


if((sim::time +1)%vout::output_rate==0){

   for (int atom = 0; atom < atoms::num_atoms; atom ++)
   {
      int cell = atoms::cell_array[atom];
      atoms::x_spin_array[atom] = mm::x_array[cell];
      atoms::y_spin_array[atom] = mm::y_array[cell];
      atoms::z_spin_array[atom] = mm::z_array[cell];
   }


   //calcualtes the average magnetisation of the system per step.
/*
   double m_x = 0;
   double m_y = 0;
   double m_z = 0;
   double m_l = 0;

   for (int cell = 0; cell < num_cells; cell++)
   {
      m_x = m_x + cells::x_mag_array[cell];
      m_y = m_y + cells::y_mag_array[cell];
      m_z = m_z + cells::z_mag_array[cell];
      m_l = m_l + mm::ms[cell];
   }

   m_x = m_x/m_l;
   m_y = m_y/m_l;
   m_z = m_z/m_l;

   m_l = sqrt(m_x*m_x + m_y*m_y + m_z*m_z);

   std::cout << sim::time << '\t' << temperature << '\t' << m_x << '\t' << m_y << '\t' << m_z  << "\t" << m_l <<std::endl;
*/
}
   if (micromagnetic::enable_boltzman_distribution == true){
      if (sim::time > 1000){
         for (int cell = 0; cell < num_cells; cell ++)
         {
            if (mm::ms[cell] > 1e-100)
            {
            double mz=sqrt(mm::z_array[cell]*mm::z_array[cell]);
            double mx=sqrt(mm::x_array[cell]*mm::x_array[cell] + mm::y_array[cell]*mm::y_array[cell]); //sqrt(S[cell]*S[cell]+S[1]*S[1]);
            double mag_m =sqrt(mm::x_array[cell]*mm::x_array[cell] + mm::y_array[cell]*mm::y_array[cell]+ mm::z_array[cell]*mm::z_array[cell]);
            int para = int(mz*100.0+0.5);
            int perp = int(mx*100.0+0.5);
            if (para >100) para =100;
            if (perp >100) perp =100;
            int para1D = int(mag_m*1000.0+0.5);
            micromagnetic::P[para][perp]++;
            micromagnetic::P1D[para1D]++;
            micromagnetic::mean_M+=mag_m;
            micromagnetic::counter++;
         }
         }
      }
   }
}
