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

#include "internal.hpp"

#include "cells.hpp"
#include "../cells/internal.hpp"

// micromagnetic module headers
#include <vector>

// typesafe sign function
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}


namespace micromagnetic{

   namespace internal{

int calculate_bias_magnets(double system_dimensions_x,double system_dimensions_y,double system_dimensions_z){

  //THIS IS CORRECT I KNOW THE X AND Y ARE MISSED UP! THIS IS BECAUE MS is assumed to be along y but i need it along x so i have switched x and y then switched them back at the end

  double shield_Ms = bias_magnet_ms_input;
  double x_size = system_dimensions_y*bias_magnets_max_width - system_dimensions_y*bias_magnets_min_width;
  double y_size = 1000000;
  double z_size = system_dimensions_z*bias_magnets_max_height - system_dimensions_z*bias_magnets_min_height;

  double x_pos = x_size/2 + system_dimensions_y*bias_magnets_min_width;
  double y_pos;
  double z_pos = z_size/2.0 + system_dimensions_z*bias_magnets_min_height;

  double y_pos_1 = - y_size/2 - bias_magnets_gap;
  double y_pos_2 =   y_size/2 + system_dimensions_x + bias_magnets_gap;
//

//  std::cout << y_pos_1 << '\t' << y_pos_2 << '\t' << x_pos << '\t' << z_pos << std::endl;

   double prefactor = shield_Ms/(4.0*M_PI);
  //save this new m as the initial value, so it can be saved and used in the final equation.

  for (int lc = 0; lc < cells::num_local_cells; lc++){

    int cell = cells::local_cell_array[lc];



    bias_field_x[cell] = 0;
    bias_field_y[cell] = 0;
    bias_field_z[cell] = 0;

     //cell position in Angstrom
     double x_cell = cells::pos_and_mom_array[4*cell + 0];
     double y_cell = cells::pos_and_mom_array[4*cell + 1];
     double z_cell = cells::pos_and_mom_array[4*cell + 2];


     const double xb = x_size * 0.5;
     const double yb = y_size * 0.5;
     const double zb = z_size * 0.5;


     for (int shield = 0; shield < 2; shield++){

       if (shield == 0) y_pos = y_pos_1;
       if (shield == 1) y_pos = y_pos_2;
       //calculates the vector in A from the cell to the shields
       double x = sqrt((y_cell - x_pos)*(y_cell - x_pos));
       double y = sqrt((x_cell - y_pos)*(x_cell - y_pos));
       double z = sqrt((z_cell - z_pos)*(z_cell - z_pos));

       double Bx = 0.0;
       double By = 0.0;
       double Bz = 0.0;

       for(int k=1; k<3; k++){

           // predefine power as fixed for loop iteration
           const double m1k = pow(-1,k);

           for(int l=1; l<3; l++){

              // predefine power as fixed for loop iteration
              const double m1l = pow(-1,l);

              for(int m=1; m<3; m++){

                 const double m1m = pow(-1,m);
                 const double m1klm = pow(-1,k+l+m);

                 const double xp = x + xb*m1k;
                 const double yp = y + yb*m1l;
                 const double zp = z + zb*m1m;

                 const double xabs = fabs(xp);
                 const double yabs = fabs(yp);

                 double r = sqrt(xp*xp + yp*yp + zp*zp);

                 Bx = Bx + m1klm * log(zp + r);
                 By = By - m1klm * sign(yp) * sign(xp) * atan(xabs * zp / (yabs * r));
                 Bz = Bz + m1klm * log(xp + r);


              }
           }
       }
   //     std::cout << cell << '\t' << shield << '\t' <<  bias_field_x[cell] << '\t' << By*prefactor << std::endl;

      // std::cout << shield << '\t' << cell << '\t' << Bx << "\t" << By << '\t' << Bz <<  '\t' << prefactor << std::endl;
       bias_field_x[cell] = bias_field_x[cell] + By*prefactor;
       bias_field_y[cell] = bias_field_y[cell] + Bx*prefactor;
       bias_field_z[cell] = bias_field_z[cell] + Bz*prefactor;
   //    std::cout << cell << '\t' << shield << '\t' <<  bias_field_x[cell] << '\t' << By*prefactor << std::endl;

     }
    // std::cout << cell << '\t' << cells::pos_and_mom_array[4*cell+0] << '\t' << cells::pos_and_mom_array[4*cell+1] << '\t' << cells::pos_and_mom_array[4*cell+2] << '\t'  << bias_field_x[cell] << '\t' << bias_field_y[cell] << '\t' << bias_field_z[cell] <<  "\t" << z_cell << std::endl;

  }


return 0;
 }
}
}
