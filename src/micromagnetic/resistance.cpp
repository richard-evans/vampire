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
       //double one_o_Rmin = 1.0/Rmin;
       double sum_one_o_R = 0.0;

 //   std::cout <<res_RA << '\t' << Rmin <<"\t" <<  overlap_area <<'\t' << cells::num_cells << '\t' << Ra << '\t' << area << std::endl;
       int i = 0;
      for (int  cell = 0; cell < cells::num_cells; cell++){

         int mat  = cell_material_array[cell];
         if (mat == resistance_layer_1 ){

           const int start = macro_neighbour_list_start_index[cell];
           const int end = macro_neighbour_list_end_index[cell] +1;

           mx_i = cells::mag_array_x[cell];
           my_i = cells::mag_array_y[cell];
           mz_i = cells::mag_array_z[cell];

           mod_i =sqrt(mx_i*mx_i+my_i*my_i +mz_i*mz_i);

           for(int j = start;j< end;j++){

            const int cellj = macro_neighbour_list_array[j];

            int matj =cell_material_array[cellj];

            if (mat == resistance_layer_1 && matj == resistance_layer_2){

               // calculate distance between cells in x,y
               double dx = cells::pos_array[cell*3 +0] - cells::pos_array[cellj*3 +0];
               double dy = cells::pos_array[cell*3 +1] - cells::pos_array[cellj*3 +1];

               // check that cells are on top of each other
               if (dx*dx < cs::unit_cell.dimensions[0]*cs::unit_cell.dimensions[0] && dy*dy < cs::unit_cell.dimensions[1]*cs::unit_cell.dimensions[1]){


          //    std::cout << cell << '\t' << cellj << '\t' << dot_product << '\t' << costheta << '\t' << R << '\t' <<sum_one_o_R <<std::endl;
              mx_j = cells::mag_array_x[cellj];
              my_j = cells::mag_array_y[cellj];
              mz_j = cells::mag_array_z[cellj];

              mod_j =sqrt(mx_j*mx_j+my_j*my_j +mz_j*mz_j);

              double dot_product = mx_i*mx_j + my_i*my_j + mz_i*mz_j;
              double costheta = dot_product/(mod_i*mod_j);
              double change  = 1- GMR_o_2*costheta;
              double R = Rmin*(change);

              sum_one_o_R = sum_one_o_R + 1.0/R;
              //std::cout << cell << '\t' << cellj << '\t' << dot_product << '\t' << costheta << '\t' << R << '\t' <<sum_one_o_R <<std::endl;
              i ++ ;
           }
            }
          }
        }
      }
//      std::cout << "a" <<std::endl;
      return i*1.0/sum_one_o_R;
    }
   }
 }
