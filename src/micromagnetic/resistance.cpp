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

       double mx_i,my_i,mz_i;
       double mx_j,my_j,mz_j;
       double mod_i, mod_j;

       double GMR = res_GMR;
       double GMR_o_2 = GMR/2.0;

       double Ra = res_RA;
       double area = overlap_area*1e-8;

       double Rmin = Ra/area;
       double one_o_Rmin = 1.0/Rmin;
       double sum_one_o_R = 0.0;

       int i = 0;
       for (int lc = 0; lc < cells::num_local_cells; lc++){

         int cell = cells::local_cell_array[lc];
         int mat  = cell_material_array[cell];
         if (mat == resistance_layer_1 ){

           const int start = macro_neighbour_list_start_index[cell];
           const int end = macro_neighbour_list_end_index[cell] +1;

           mx_i = cells::mag_array_x[cell];
           my_i = cells::mag_array_y[cell];
           mz_i = cells::mag_array_z[cell];

           mod_i =sqrt(sqrt(mx_i*mx_i+my_i*my_i)*sqrt(mx_i*mx_i+my_i*my_i) +mz_i*mz_i);

           for(int j = start;j< end;j++){

            const int cellj = macro_neighbour_list_array[j];

            int matj =cell_material_array[cellj];

            if (mat == resistance_layer_1 && matj == resistance_layer_2){
          //    std::cout << cell << '\t' << cellj << '\t' << dot_product << '\t' << costheta << '\t' << R << '\t' <<sum_one_o_R <<std::endl;
              mx_j = cells::mag_array_x[cellj];
              my_j = cells::mag_array_y[cellj];
              mz_j = cells::mag_array_z[cellj];

              mod_j =sqrt(sqrt(mx_j*mx_j+my_j*my_j)*sqrt(mx_j*mx_j+my_j*my_j) +mz_j*mz_j);

              double dot_product = mx_i*mx_j + my_i*my_j + mz_i*mz_j;
               double costheta = dot_product/(mod_i*mod_j);
              double change  = 1- GMR_o_2*costheta;
            //  std::cout << change <<std::endl;
              double R = Rmin*(change);
            //  std::cout << i << '\t' << costheta << '\t' << change << "\t" << R<<  "\t" << 1.0/R << "\t" << sum_one_o_R << '\t' << 1.0/sum_one_o_R << std::endl;

              sum_one_o_R = sum_one_o_R + 1.0/R;
              i ++ ;
          //    std::cout << dot_product << '\t' << costheta << '\t' << change << '\t' << R << sum_one_o_R << "\t" << 1/sum_one_o_R << std::endl;
            //  std::cout <<"x=" << x_i << '\t' << "y=" <<  y_i << '\t' << "z=" << z_i << "\t" << "x=" << x_j << '\t'<< "y=" << y_j << '\t' << "z=" << z_j << "\t" << mod_i*mod_j << '\t' << dot_product << '\t' << costheta << '\t' << R << '\t' <<sum_one_o_R <<std::endl;
          //  std::cout << "r" <<resistance_layer_1 << '\t' << resistance_layer_2 << "\t" << cell << '\t' << std::endl;
            }
          }
        }
      }
    //  std::cin.get();
    //  std::cout << 1.0/sum_one_o_R <<std::endl;
      return i*1.0/sum_one_o_R;
    }
   }
 }
