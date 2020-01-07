//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
// Vampire headers
#include "environment.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "create.hpp"

// C++ headers
#include <math.h>


// typesafe sign function
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}


namespace environment{

namespace internal{

bool  in_x(double x, double z){

   const double xr = x - 500.0; // shift to zero
   const double zr = z - 300.0; // shift to base of side shield

   const double zmin = 200.0 * exp(-((fabs(xr)-190)*0.01));
   //std::cout << x << '\t' << z << '\t'  << xr << '\t' << zr << '\t' << zmin << std::endl;
   if(zr > zmin && z < 500.0) return true;
   return false;

}

//------------------------------------------------------------------------------
// Function to calculate basic shield geometry for reader
//------------------------------------------------------------------------------
int in_shield(double x, double y, double z){


//  if (square_shields){
  //  std::cout << dim[1] << "\t" << cs::system_dimensions[1] << '\t' << gap << std::endl;
  //   int min = dim[0]/2.0 - cs::system_dimensions[0]/2.0 - gap;
  //   int max = min + cs::system_dimensions[0] + 2.0*gap;
  // //  std::cout << min << '\t' << max << std::endl;
  //   if (y < min || y > max) return true;
  //   else return false;

//   }
//
//
//   else if (expoential_shields){
//
//    // height of inner sensor region
//    const double stack_height = 200; // Angstroms
//
//    const double xr = x;
//    const double yr = y;
//    const double zr = z; // reduced height
//
// //   Bottom shield
   if(z <= 302 && z >= -2)    return 1;
  //
  //if (zr < 310  && zr > 290) return true;
  //
  // // Top shield (++z) 31-51 nm
  if(z >= 518 && z <= 722.9) return 2;
  //
  //  //Top shield (+z) 52-72 nm
  if(z >= 759.8 && z <= 942.4) return 3;
  // //
  // //  // side shields
  if(in_x(x, z)) return 4;

  //if (sqrt(x*x + y*y + z*z)  < 90) return 1;
  else return 0;

}

int bias_shields(){

  double shield_Ms = 1;
  double x_size = dim[1];
  double y_size = 1000000;
  double z_size = dim[2];

  double x_pos = x_size/2.0;
  double y_pos;
  double z_pos = dim[2]/2.0;

  double y_pos_1 = -y_size/2.0;
  double y_pos_2 =  y_size/2.0 +dim[0];


   double prefactor = shield_Ms/(4.0*M_PI);
  //save this new m as the initial value, so it can be saved and used in the final equation.
    for (int cell = 0; cell < num_cells; cell ++){

    std::vector <double > B(3,0.0);
    bias_field_x[cell] = 0;
    bias_field_y[cell] = 0;
    bias_field_z[cell] = 0;

     //cell position in Angstrom
     double x_cell = cell_coords_array_x[cell];
     double y_cell = cell_coords_array_y[cell];
     double z_cell = cell_coords_array_z[cell];

     const double xb = x_size * 0.5;
     const double yb = y_size * 0.5;
     const double zb = z_size * 0.5;

     for (int shield = 0; shield < 2; shield++){

       if (shield == 0) y_pos = y_pos_1;
       if (shield == 1) y_pos = y_pos_2;
       //calculates the vector in A from the cell to the shields
       double x = sqrt((x_cell - x_pos)*(x_cell - x_pos));
       double y = sqrt((y_cell - y_pos)*(y_cell - y_pos));
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
                 By = By + m1klm * sign(yp) * sign(xp) * atan(xabs * zp / (yabs * r));
                 Bz = Bz + m1klm * log(xp + r);


              }
           }
       }
       bias_field_x[cell] = bias_field_x[cell] + By*prefactor;
       bias_field_y[cell] = bias_field_y[cell] + Bx*prefactor;
       bias_field_z[cell] = bias_field_z[cell] + Bz*prefactor;

     }
  //  std::cout << bias_field_x[cell] << '\t' << bias_field_y[cell] << '\t' << bias_field_z[cell] << std::endl;

  }
//std::cin.get();



  return 0;

}


}

}
