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

// set up local namespace nickname
namespace env = environment::internal;

namespace environment{

   namespace internal{

      int output(){

         double mx = 0;
         double my = 0;
         double mz = 0;
         double ml = 0;

         //calcualtes the mean mx,my,mz,ml for all cells.
         for (int cell = 0; cell < num_cells; cell++){

            mx =  mx + x_mag_array[cell];
            my =  my + y_mag_array[cell];
            mz =  mz + z_mag_array[cell];
            ml =  ml + Ms;
         }

         double msat = ml;
         double magm = sqrt(mx*mx + my*my + mz*mz);
         mx = mx/magm;
         my = my/magm;
         mz = mz/magm;

         //outputs to the file environment_output
         o_file <<sim::time << '\t' << sim::temperature << "\t" << mx << '\t' << my<< '\t' << mz << '\t' <<  magm/msat << std::endl;


        //  std::stringstream filename_sstr;
        //  filename_sstr << "env_cell_config" << sim::time << ".txt";
        //  std::ofstream pfile;
        //  pfile.open(filename_sstr.str());
         //
        //  //
        //  for (int i = 0; i < env::num_env_cells; i++){
        //     int cell = env::none_atomistic_cells[i];
        //  //for(int cell = 0; cell < num_cells; cell++){
        //  //
        //  	pfile << cell_coords_array_x[cell] << '\t' << cell_coords_array_y[cell] << '\t' << cell_coords_array_z[cell] << '\t' <<x_mag_array[cell] << '\t' << y_mag_array[cell] << '\t' << z_mag_array[cell] << '\t' << std::endl;
        //  }

         return 0;
      }
   }
}
