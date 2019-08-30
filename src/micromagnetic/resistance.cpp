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

// micromagnetic module headers
#include <vector>



namespace micromagnetic{

   namespace internal{

     double calculate_resistance(){

       double x_i,y_i,z_i;
       double x_j,y_j,z_j;
       double mod_i, mod_j;

       double GMR = 1.0;
       double GMR_o_2 = GMR/2.0;

       double Rmin = 1.0;
       double one_o_Rmin = 1.0/Rmin;
       double sum_one_o_R = 0.0;

       for (int lc = 0; lc < cells::num_local_cells; lc++){

         int cell = cells::local_cell_array[lc];
         int mat  = cell_material_array[cell];
         if (mat == resistance_layer_1 ){

           const int start = macro_neighbour_list_start_index[cell];
           const int end = macro_neighbour_list_end_index[cell] +1;

           x_i = cells::mag_array_x[cell];
           y_i = cells::mag_array_y[cell];
           z_i = cells::mag_array_z[cell];

           mod_i =sqrt(sqrt(x_i*x_i+y_i*y_i)*sqrt(x_i*x_i+y_i*y_i) +z_i*z_i);

           for(int j = start;j< end;j++){

            const int cellj = macro_neighbour_list_array[j];

            int matj =cell_material_array[cellj];

            if (mat == resistance_layer_1 && matj == resistance_layer_2){
          //    std::cout << cell << '\t' << cellj << '\t' << dot_product << '\t' << costheta << '\t' << R << '\t' <<sum_one_o_R <<std::endl;
              x_j = cells::mag_array_x[cellj];
              y_j = cells::mag_array_y[cellj];
              z_j = cells::mag_array_z[cellj];

              mod_j =sqrt(sqrt(x_j*x_j+y_j*y_j)*sqrt(x_j*x_j+y_j*y_j) +z_j*z_j);

              double dot_product = x_i*x_j + y_i*y_j + z_i*z_j;

              double costheta = dot_product/(mod_i*mod_j);
              double change  = 1- GMR_o_2*costheta;
            //  std::cout << change <<std::endl;
              double R = Rmin*(change);
          //    std::cout << costheta << '\t' << change << std::endl;
              sum_one_o_R = sum_one_o_R + 1.0/R;
            //  std::cout <<"x=" << x_i << '\t' << "y=" <<  y_i << '\t' << "z=" << z_i << "\t" << "x=" << x_j << '\t'<< "y=" << y_j << '\t' << "z=" << z_j << "\t" << mod_i*mod_j << '\t' << dot_product << '\t' << costheta << '\t' << R << '\t' <<sum_one_o_R <<std::endl;
          //  std::cout << "r" <<resistance_layer_1 << '\t' << resistance_layer_2 << "\t" << cell << '\t' << std::endl;
            }
          }
        }
      }
    //  std::cout << 1.0/sum_one_o_R <<std::endl;
      return 1.0/sum_one_o_R;
    }
   }
 }
