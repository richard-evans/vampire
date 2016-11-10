
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

   void mm::step(int num_cells, double temperature, std::vector<double> x_array,std::vector<double> y_array,std::vector<double> z_array, std::vector<double> ext_field, double dt,std::vector<double>& new_x_array,std::vector<double>& new_y_array,std::vector<double>& new_z_array){

      const double kB = 1.3806503e-23;
      double one_o_chi_perp, one_o_2_chi_para, reduced_temperature, Tc_o_Tc_m_T, m_e, alpha_para, alpha_perp, m_e_squared;

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
         //if (temperature < mm::Tc[cell]){
            int j2 = cell*num_cells;
            //loops over all other cells to sum the interaction
            double mi = pow(m_e_squared,0.5);
            for(int j=0;j<num_cells;j++){
               // calculate magnetization length
               const double mj = sqrt(x_array[j]*x_array[j] +y_array[j]*y_array[j] + z_array[j]*z_array[j]);
               // calculate reduced exchange constant factor
               const double A = mm::Ax[j2]*pow(mj,1.66);

               //if(cell == 63 && mm::Ax[j2] < 0.0){
               //   std::cout << cell << "\t" << j2-cell*num_cells << "\t" << mm::Ax[j2] << std::endl;
               //   std::cin.get();
               //}
               //if (mj > 0){
                  //H_exch = sum_j A(T)
                  //exchange_field[0] = exchange_field[0] - (pow(mi,1.66/2.)*pow(mj,1.66/2.))*(x_array[j] - x_array[cell])*(x_array[j])*mm::Ax[j2]/(mj*mj);
                  //exchange_field[1] = exchange_field[1] - (pow(mi,1.66/2.)*pow(mj,1.66/2.))*(y_array[j] - y_array[cell])*(y_array[j])*mm::Ay[j2]/(mj*mj);
                  //exchange_field[2] = exchange_field[2] - (pow(mi,1.66/2.)*pow(mj,1.66/2.))*(z_array[j] - z_array[cell])*(z_array[j])*mm::Az[j2]/(mj*mj);
                  exchange_field[0] -= A*(x_array[j] - x_array[cell]);
                  exchange_field[1] -= A*(y_array[j] - y_array[cell]); //- (pow(mi,1.66/2.)*pow(mj,1.66/2.))*(y_array[j] - y_array[cell])*(y_array[j])*mm::Ay[j2]/(mj*mj);
                  exchange_field[2] -= A*(z_array[j] - z_array[cell]); //- (pow(mi,1.66/2.)*pow(mj,1.66/2.))*(z_array[j] - z_array[cell])*(z_array[j])*mm::Az[j2]/(mj*mj);
                  //sumx = sumx - mm::Ax[j2];
                  //sumy = sumy - mm::Ay[j2];
                  //sumz = sumz - mm::Az[j2];
   //             exchange_field[0] = exchange_field[0] - (pow(mi,1.66/2.)*pow(mj,1.66/2.))*(x_array[j])*mm::Ax[j2]/(mj*mj);
      //          exchange_field[1] = exchange_field[1] - (pow(mi,1.66/2.)*pow(mj,1.66/2.))*(y_array[j])*mm::Ay[j2]/(mj*mj);
         //       exchange_field[2] = exchange_field[2] - (pow(mi,1.66/2.)*pow(mj,1.66/2.))*(z_array[j])*mm::Az[j2]/(mj*mj);

               //}
               j2++;
            }
         //}


         //Sum H = H_exch + H_A +H_exch_grains +H_App + H+dip
         spin_field[0] = pf*m[0] + one_o_chi_perp*m[0] + 1*exchange_field[0] + mm::ext_field[0];
         spin_field[1] = pf*m[1] + one_o_chi_perp*m[1] + 1*exchange_field[1] + mm::ext_field[1];
         spin_field[2] = pf*m[2]                       + 1*exchange_field[2] + mm::ext_field[2];

         //std::cout << sim::time << '\t' << temperature << '\t' << m[0] << '\t' << m[1] << '\t' << m[2] << '\t' << pf*m[0] << '\t' << sumx << '\t' << sumy << '\t' << sumz << '\t' << pf*m[2] << "\t" << one_o_chi_perp*m[0] << '\t' << one_o_chi_perp*m[1] << "\t" << exchange_field[0] << "\t" << exchange_field[1] << "\t" << exchange_field[2] << std::endl;
         //calculates the stochatic parallel and perpendicular terms
         double sigma_para = sqrt(2*kB*temperature*alpha_para/(mm::ms[cell]*dt)); //why 1e-27
         double sigma_perp = sqrt(2*kB*temperature*(alpha_perp-alpha_para)/(dt*mm::ms[cell]*alpha_perp*alpha_perp));

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

/*
void demag_FFT(int num_cells,std::vector<double> x_array,std::vector<double> y_array,std::vector<double> z_array,std::vector<double>& dip_field_x,std::vector<double>& dip_field_y,std::vector<double>& dip_field_z)
{

   Array3D<fftw_complex> Mx; //3D Array for magneetisation
   Array3D<fftw_complex> My;
   Array3D<fftw_complex> Mz;

   Array3D<fftw_complex> Hx; //3D Array for dipolar field
   Array3D<fftw_complex> Hy;
   Array3D<fftw_complex> Hz;

   Hx.resize(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z);
   Hy.resize(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z);
   Hz.resize(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z);

   Mx.resize(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z);
   My.resize(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z);
   Mz.resize(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z);

   for (int cell = 0; cell < num_cells; cell++)

      Mx.IFillReal(x_array[cell]);
      Mx.IFillComplex(0.0);
      My.IFillReal(y_array[cell]);
      My.IFillComplex(0.0);
      Mz.IFillReal(z_array[cell]);
      Mz.IFillComplex(0.0);

      Hx.IFill(0.0);
      Hy.IFill(0.0);
      Hz.IFill(0.0);


      // fft calculations
      fftw_plan NxxP,NxyP,NxzP,NyxP,NyyP,NyzP,NzxP,NzyP,NzzP;
      fftw_plan MxP,MyP,MzP;

      //deterines the forward transform for the N arrays
      NxxP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Nxx.ptr(),Nxx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(NxxP);
      NyxP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Nyx.ptr(),Nyx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(NyxP);
      NzxP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Nzx.ptr(),Nzx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(NzxP);

      NxyP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Nxy.ptr(),Nxy.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(NxyP);
      NyyP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Nyy.ptr(),Nyy.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(NyyP);
      NzyP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Nzy.ptr(),Nzy.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(NzyP);

      NxzP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Nxz.ptr(),Nxz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(NxzP);
      NyzP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Nyz.ptr(),Nyz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(NyzP);
      NzzP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Nzz.ptr(),Nzz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(NzzP);

      MxP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Mx.ptr(),Mx.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(MxP);
      MyP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,My.ptr(),My.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(MyP);
      MzP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Mz.ptr(),Mz.ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(MzP);

      // performs the converlusion between Nk and Mk
      for (int i=0 ; i<2*mm::num_macro_cells_x ; i++){
      for (int j=0 ; j<2*mm::num_macro_cells_y ; j++){
          for (int k=0 ; k<2*mm::num_macro_cells_z ; k++){
           // [Nreal+ iNimag]*(Mreal+iMimag)
           Hx(i,j,k)[0] = Nxx(i,j,k)[0]*Mx(i,j,k)[0] + Nxy(i,j,k)[0]*My(i,j,k)[0] + Nxz(i,j,k)[0]*Mz(i,j,k)[0]; //summing the real part
           Hx(i,j,k)[0] -= (Nxx(i,j,k)[1]*Mx(i,j,k)[1] + Nxy(i,j,k)[1]*My(i,j,k)[1] + Nxz(i,j,k)[1]*Mz(i,j,k)[1]);
           Hx(i,j,k)[1] = Nxx(i,j,k)[0]*Mx(i,j,k)[1] + Nxy(i,j,k)[0]*My(i,j,k)[1] + Nxz(i,j,k)[0]*Mz(i,j,k)[1];
           Hx(i,j,k)[1] += (Nxx(i,j,k)[1]*Mx(i,j,k)[0] + Nxy(i,j,k)[1]*My(i,j,k)[0] + Nxz(i,j,k)[1]*Mz(i,j,k)[0]);

           Hy(i,j,k)[0] = Nyx(i,j,k)[0]*Mx(i,j,k)[0] + Nyy(i,j,k)[0]*My(i,j,k)[0] + Nyz(i,j,k)[0]*Mz(i,j,k)[0];
           Hy(i,j,k)[0] -= (Nyx(i,j,k)[1]*Mx(i,j,k)[1] + Nyy(i,j,k)[1]*My(i,j,k)[1] + Nyz(i,j,k)[1]*Mz(i,j,k)[1]);
           Hy(i,j,k)[1] = Nyx(i,j,k)[0]*Mx(i,j,k)[1] + Nyy(i,j,k)[0]*My(i,j,k)[1] + Nyz(i,j,k)[0]*Mz(i,j,k)[1];
           Hy(i,j,k)[1] += (Nyx(i,j,k)[1]*Mx(i,j,k)[0] + Nyy(i,j,k)[1]*My(i,j,k)[0] + Nyz(i,j,k)[1]*Mz(i,j,k)[0]);

           Hz(i,j,k)[0] = Nzx(i,j,k)[0]*Mx(i,j,k)[0] + Nzy(i,j,k)[0]*My(i,j,k)[0] + Nzz(i,j,k)[0]*Mz(i,j,k)[0]; //summing the real part
           Hz(i,j,k)[0] -= (Nzx(i,j,k)[1]*Mx(i,j,k)[1] + Nzy(i,j,k)[1]*My(i,j,k)[1] + Nzz(i,j,k)[1]*Mz(i,j,k)[1]);
           Hz(i,j,k)[1] = Nzx(i,j,k)[0]*Mx(i,j,k)[1] + Nzy(i,j,k)[0]*My(i,j,k)[1] + Nzz(i,j,k)[0]*Mz(i,j,k)[1];
           Hz(i,j,k)[1] += (Nzx(i,j,k)[1]*Mx(i,j,k)[0] + Nzy(i,j,k)[1]*My(i,j,k)[0] + Nzz(i,j,k)[1]*Mz(i,j,k)[0]);
          }
      }
   }

   // performs the backward transform to give the dipole field, Hx, Hy, Hz
   fftw_plan HxP,HyP,HzP;

   HxP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Hx.ptr(),Hx.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
   fftw_execute(HxP);
   HyP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Hy.ptr(),Hy.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
   fftw_execute(HyP);
   HzP = fftw_plan_dft_3d(2*mm::num_macro_cells_x,2*mm::num_macro_cells_y,2*mm::num_macro_cells_z,Hz.ptr(),Hz.ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
   fftw_execute(HzP);

   for (int i=0 ; i<2*mm::num_macro_cells_x ; i++){
   for (int j=0 ; j<2*mm::num_macro_cells_y ; j++){
      for (int k=0 ; k<2*mm::num_macro_cells_z ; k++){
          if(i==j && i==k){
           Hx(i,j,k)[0] += Mx(i,j,k)[0]*8.0*pi/3.0;
           Hy(i,j,k)[0] += My(i,j,k)[0]*8.0*pi/3.0;
           Hz(i,j,k)[0] += Mz(i,j,k)[0]*8.0*pi/3.0;
           }
        myfile << i  << "\t" << j  << "\t" << k << "\t" << 1.0e-7*Ms*Hx(i,j,k)[0]/((3.0*mm::num_macro_cells_x)*(3.0*mm::num_macro_cells_y)*(3.0*mm::num_macro_cells_z)) << "\t" << 1.0e-7*Ms*Hy(i,j,k)[0]/((3.0*mm::num_macro_cells_x)*(3.0*mm::num_macro_cells_y)*(3.0*mm::num_macro_cells_z)) << "\t" << 1.0e-7*Ms*Hz(i,j,k)[0]/((1.0*mm::num_macro_cells_x)*(1.0*mm::num_macro_cells_y)*(1.0*mm::num_macro_cells_z)) << "\n";
      }
   }
   }

   return 0;


}*/
}
