//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//   implementation based on the paper Phys. Rev. B 103, 024429, (2021) M.Strungaru, M.O.A. Ellis et al
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include <iostream>
#include <vector>
#include <math.h>
#include "atoms.hpp"
#include "create.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "random.hpp"
#include "errors.hpp"
#include "internal.hpp"


namespace sim{
              


   int STDspin(){
      const int num_atoms=atoms::num_atoms;
      double cay_dt=-mp::dt/2.0;//-dt4*consts::gyro - mp::dt contains gamma;
      double dt2=0.5*mp::dt_SI*1e12;
    

      std::vector <double> Hx_th(atoms::x_spin_array.size());
   	  std::vector <double> Hy_th(atoms::x_spin_array.size());
   	  std::vector <double> Hz_th(atoms::x_spin_array.size());

      generate (Hx_th.begin(),Hx_th.end(), mtrandom::gaussian);
      generate (Hy_th.begin(),Hy_th.end(), mtrandom::gaussian);
      generate (Hz_th.begin(),Hz_th.end(), mtrandom::gaussian);


      for(int atom=0;atom<num_atoms;atom++){
      
      calculate_spin_fields(atom, atom+1);
  	  calculate_external_fields(atom, atom+1);

     /* atoms::x_total_spin_field_array[atom]=+atoms::x_total_external_field_array[atom];
      atoms::y_total_spin_field_array[atom]=+atoms::y_total_external_field_array[atom];
      atoms::z_total_spin_field_array[atom]=+atoms::z_total_external_field_array[atom];*/

      
      add_spin_noise(atom,
                  atom+1,
                  mp::dt_SI*1e12,
                  atoms::type_array, // type for atom
                  atoms::x_spin_array,
                  atoms::y_spin_array,
                  atoms::z_spin_array,
                  atoms::x_total_spin_field_array,
                  atoms::y_total_spin_field_array,
                  atoms::z_total_spin_field_array,
                  Hx_th, //  vectors for fields
                  Hy_th,
                  Hz_th);
                  



     cayley_update(atom,
                  atom+1,
                  cay_dt,
                  atoms::x_spin_array,
                  atoms::y_spin_array,
                  atoms::z_spin_array,
                  atoms::x_total_spin_field_array,
                  atoms::y_total_spin_field_array,
                  atoms::z_total_spin_field_array);
                  


      }
      
       for(int atom=num_atoms-1;atom>=0;atom--){
       
       
       calculate_spin_fields(atom, atom+1);
   	   calculate_external_fields(atom,atom+1);
       /*atoms::x_total_spin_field_array[atom]=+atoms::x_total_external_field_array[atom];
       atoms::y_total_spin_field_array[atom]=+atoms::y_total_external_field_array[atom];
       atoms::z_total_spin_field_array[atom]=+atoms::z_total_external_field_array[atom];*/

   	   
       add_spin_noise(atom,
                   atom+1,
                   mp::dt_SI*1e12,
                   atoms::type_array, // type for atom
                   atoms::x_spin_array,
                   atoms::y_spin_array,
                   atoms::z_spin_array,
                   atoms::x_total_spin_field_array,
                   atoms::y_total_spin_field_array,
                   atoms::z_total_spin_field_array,
                   Hx_th, 
                   Hy_th,
                   Hz_th);


      cayley_update(atom,
                   atom+1,
                   cay_dt,
                   atoms::x_spin_array,
                   atoms::y_spin_array,
                   atoms::z_spin_array,
                   atoms::x_total_spin_field_array,
                   atoms::y_total_spin_field_array,
                   atoms::z_total_spin_field_array);

       }

      
      return EXIT_SUCCESS;
  }


void cayley_update(const int start_index,
            const int end_index,
            double dt,
            std::vector<double>& x_spin_array, // coord vectors for atoms
            std::vector<double>& y_spin_array,
            std::vector<double>& z_spin_array,
            std::vector<double>& fields_array_x, //  vectors for fields
            std::vector<double>& fields_array_y,
            std::vector<double>& fields_array_z){

      for( int i = start_index; i<end_index; i++)
      {
          double Sx = x_spin_array[i];
          double Sy = y_spin_array[i];
          double Sz = z_spin_array[i];

          double Ax = fields_array_x[i] * dt;
          double Ay = fields_array_y[i] * dt;
          double Az = fields_array_z[i] * dt;

          double AS = Ax*Sx + Ay*Sy + Az*Sz;
          double A2 = Ax * Ax + Ay* Ay + Az * Az;

          double AxSx = Ay * Sz - Az * Sy;
          double AxSy = Az * Sx - Ax * Sz;
          double AxSz = Ax * Sy - Ay * Sx;

          double factor = 1.0 / (1.0 + 0.25 * A2);

          x_spin_array[i] = (Sx * ( 1.0 - 0.25 * A2) + AxSx + 0.5 * Ax * AS) * factor;
          y_spin_array[i] = (Sy * ( 1.0 - 0.25 * A2) + AxSy + 0.5 * Ay * AS) * factor;
          z_spin_array[i] = (Sz * ( 1.0 - 0.25 * A2) + AxSz + 0.5 * Az * AS) * factor;
      }

  return;
}

void add_spin_noise(const int start_index,
            const int end_index,
            double dt,
            const std::vector<int>& type_array, // type for atom
            std::vector<double>& x_spin_array, // coord vectors for atoms
            std::vector<double>& y_spin_array,
            std::vector<double>& z_spin_array,
            std::vector<double>& fields_array_x, //  vectors for fields
            std::vector<double>& fields_array_y,
            std::vector<double>& fields_array_z,
            std::vector<double>& Hx_th, //  vectors for fields
            std::vector<double>& Hy_th,
            std::vector<double>& Hz_th){

     //std::cout<<"lambda= "<<lambda<<std::endl;

     for( int i = start_index; i<end_index; i++)

    {  
       const unsigned int imat = atoms::type_array[i];
       double lambda=mp::material[imat].alpha;
       double spin_noise=mp::material[imat].H_th_sigma*sqrt(sim::temperature);
      
        double Sx = x_spin_array[i];
        double Sy = y_spin_array[i];
        double Sz = z_spin_array[i];

        double Fx = fields_array_x[i] + spin_noise * Hx_th[i];
        double Fy = fields_array_y[i] + spin_noise * Hy_th[i];
        double Fz = fields_array_z[i] + spin_noise * Hz_th[i];

        double FxSx = Fy * Sz - Fz * Sy;
        double FxSy = Fz * Sx - Fx * Sz;
        double FxSz = Fx * Sy - Fy * Sx;

        double inv_l2 = 1.0 / (1.0 + lambda*lambda);

        fields_array_x[i] = (Fx + lambda * FxSx) * inv_l2;
        fields_array_y[i] = (Fy + lambda * FxSy) * inv_l2;
        fields_array_z[i] = (Fz + lambda * FxSz) * inv_l2;
    }

return;
}//end of add_spin_noise

} // end of sim namespace
