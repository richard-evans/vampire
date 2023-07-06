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

// micromagnetic module headers
#include "internal.hpp"
#include "material.hpp"
#include "cells.hpp"
#include "sim.hpp"
#include <string>

namespace mm = micromagnetic::internal;

namespace micromagnetic{


      void outputs(){
         std::string str;
         std::string str_m;
         std::string str_time;
         std::string str_temp;
         std::string str_H;
         std::string str_R;
         for (size_t i = 0; i < mm::output_list.size(); i ++){
         //   std::cout << i << '\t' << "\t" << output_list.size() << "\t" << output_list[i] << std::endl;


      if (internal::output_list[i] == 3){

         str_time = std::to_string(sim::time) + "\t";
            str.append(str_time);
      }
      if (internal::output_list[i] == 2){

         str_temp = std::to_string(sim::temperature) + "\t";
         str.append(str_temp);
      }

      if (internal::output_list[i] == 1){

         str_H = std::to_string(sim::H_applied) + "\t";
         str.append(str_H);
      }

      if (internal::output_list[i] == 4){

         str_R = std::to_string(internal::calculate_resistance()) + "\t";
         str.append(str_R);
      }

      if (internal::output_list[i] == 0){


      std::vector <double> mx(mp::num_materials,0);
      std::vector <double> my(mp::num_materials,0);
      std::vector <double> mz(mp::num_materials,0);
      std::vector <double> ml(mp::num_materials,0);

      for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){

        int cell = list_of_micromagnetic_cells[lc];
        int mat = mm::cell_material_array[cell];
      //   std::cout << mat << "\t" << std::endl;

         double mx1 = cells::mag_array_x[cell];
         double my1 = cells::mag_array_y[cell];
         double mz1 = cells::mag_array_z[cell];
         //const double im_squared = 1.0/sqrt(mx1*mx1 + my1*my1 +mz1*mz1);
         mx[mat] = mx[mat] + mx1;//*im_squared;
         my[mat] = my[mat] + my1;//*im_squared;
         mz[mat] = mz[mat] + mz1;//*im_squared;
         // const double magn = sqrt(cells::mag_array_x[cell]*cells::mag_array_x[cell] +
         //                       cells::mag_array_y[cell]*cells::mag_array_y[cell] +
         //                       cells::mag_array_z[cell]*cells::mag_array_z[cell]);
         ml[mat] = ml[mat] + mm::ms[cell];


   }
  // std::cout << mp::num_materials << std::endl;
   for (int mat = 0; mat < mp::num_materials; mat++){
      str_m = std::to_string((mx[mat]/ml[mat])) + "\t" +  std::to_string(my[mat]/ml[mat]) + "\t" +  std::to_string(mz[mat]/ml[mat]) + "\t" +  std::to_string(sqrt(mx[mat]*mx[mat] + my[mat]*my[mat] +mz[mat]*mz[mat])/ml[mat]) + "\t" ;

   str.append(str_m);
      }



        }


      }



   internal::mm_output << str << std::endl;

      }
}
