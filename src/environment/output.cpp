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
#include "environment.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "sim.hpp"

#include<math.h>
#include<sstream>
#include<iostream>

// set up local namespace nickname
namespace env = environment::internal;

namespace environment{

   namespace internal{

      int output(){

        std::vector <double> mx(env::num_shields+1,0.0);
        std::vector <double> my(env::num_shields+1,0.0);
        std::vector <double> mz(env::num_shields+1,0.0);
         std::vector <double> ml(env::num_shields+1,0.0);

         //calcualtes the mean mx,my,mz,ml for all cells.
         for (int cell = 0; cell < num_cells; cell++){
           int shield = env::shield_number[cell];

            mx[shield] =  mx[shield] + x_mag_array[cell];
            my[shield] =  my[shield] + y_mag_array[cell];
            mz[shield] =  mz[shield] + z_mag_array[cell];
            ml[shield] =  ml[shield] + Ms[cell];
         }

         env::o_file <<sim::time << '\t' << sim::temperature << "\t";

         for (int shield = 0; shield < env::num_shields; shield++){
         double msat = ml[shield];
         double magm = sqrt(mx[shield]*mx[shield] + my[shield]*my[shield] + mz[shield]*mz[shield]);
         mx[shield] = mx[shield]/magm;
         my[shield] = my[shield]/magm;
         mz[shield] = mz[shield]/magm;
      //   std::cout << env::num_shields << '\t' << shield << "\t" << mx[shield] << '\t' << my[shield]<< '\t' << mz[shield] << '\t' <<std::endl;
        // outputs to th  e file environment_output
        env::o_file << mx[shield] << '\t' << my[shield]<< '\t' << mz[shield] << '\t' <<  magm/msat << '\t';

      }
      env::o_file << std::endl;
        if (env::env_output_info){
         std::stringstream filename_sstr;
         filename_sstr << "env_cell_config" << sim::time << ".txt";
         std::ofstream pfile;
         pfile.open(filename_sstr.str());

         //
         for (int i = 0; i < env::num_env_cells; i++){
            int cell = env::none_atomistic_cells[i];
         //for(int cell = 0; cell < num_cells; cell++){
         //
         	pfile << cell_coords_array_x[cell] << '\t' << cell_coords_array_y[cell] << '\t' << cell_coords_array_z[cell] << '\t' <<x_mag_array[cell] << '\t' << y_mag_array[cell] << '\t' << z_mag_array[cell] << '\t' << std::endl;
         }
       }

         return 0;
      }
   }
}
