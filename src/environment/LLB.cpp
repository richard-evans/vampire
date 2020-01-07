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
#include "random.hpp"
#include "sim.hpp"
#include "cells.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "micromagnetic.hpp"


#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>

namespace env = environment::internal;

namespace LLB_arrays{

   // Local arrays for LLG integration
   std::vector <double> x_euler_array;
   std::vector <double> y_euler_array;
   std::vector <double> z_euler_array;

   std::vector <double> x_array;
   std::vector <double> y_array;
   std::vector <double> z_array;

   std::vector <double> x_heun_array;
   std::vector <double> y_heun_array;
   std::vector <double> z_heun_array;

   std::vector <double> x_spin_storage_array;
   std::vector <double> y_spin_storage_array;
   std::vector <double> z_spin_storage_array;

   std::vector <double> x_initial_spin_array;
   std::vector <double> y_initial_spin_array;
   std::vector <double> z_initial_spin_array;

   //stores the gaussian width components of the random thermal fields
   std::vector <double> GW1x;
   std::vector <double> GW1y;
   std::vector <double> GW1z;
   std::vector <double> GW2x;
   std::vector <double> GW2y;
   std::vector <double> GW2z;


   bool LLG_set=false; ///< Flag to define state of LLG arrays (initialised/uninitialised)

}


namespace environment{

   int LLB_init(int num_cells){

      using namespace LLB_arrays;

      //resize arrays
      x_spin_storage_array.resize(num_cells,0.0);
      y_spin_storage_array.resize(num_cells,0.0);
      z_spin_storage_array.resize(num_cells,0.0);

      x_array.resize(num_cells,0.0);
      y_array.resize(num_cells,0.0);
      z_array.resize(num_cells,0.0);

      x_initial_spin_array.resize(num_cells,0.0);
      y_initial_spin_array.resize(num_cells,0.0);
      z_initial_spin_array.resize(num_cells,0.0);

      x_euler_array.resize(num_cells,0.0);
      y_euler_array.resize(num_cells,0.0);
      z_euler_array.resize(num_cells,0.0);

      x_heun_array.resize(num_cells,0.0);
      y_heun_array.resize(num_cells,0.0);
      z_heun_array.resize(num_cells,0.0);

      GW1x.resize(num_cells,0.0);
      GW1y.resize(num_cells,0.0);
      GW1z.resize(num_cells,0.0);
      GW2x.resize(num_cells,0.0);
      GW2y.resize(num_cells,0.0);
      GW2z.resize(num_cells,0.0);


      LLG_set=true;

      return EXIT_SUCCESS;
   }

   int LLB(double temperature,
      double H,
      double Hx,
      double Hy,
      double Hz,
      double dt){

         using namespace LLB_arrays;
               //updating the cell magnetisation (in parallel)
               if (micromagnetic::discretisation_type == 0) cells::mag();

               int num_cells = env::num_cells;

               // set up tranges for different processors
               int num_env_cells = env::num_env_cells;
               int my_num_env_cells = num_env_cells/vmpi::num_processors;
               int my_env_start_index = my_num_env_cells*vmpi::my_rank; // first cell to intergrate on local (my) cpu
               int my_env_end_index = my_num_env_cells*(vmpi::my_rank+1);  // last cell +1 to intergrate on local (my) cpu
               if (vmpi::my_rank == vmpi::num_processors - 1 ) my_env_end_index = num_env_cells;


               // Check for initialisation of LLG integration arrays
               if(LLG_set== false) environment::LLB_init(num_cells);
               // Local variables for system integration

               //sets to 0 for the parallel processing
               for (int cell = 0; cell < num_cells; cell++){
                  x_spin_storage_array[cell] = 0.0;
                  y_spin_storage_array[cell] = 0.0;
                  z_spin_storage_array[cell] = 0.0;
               }

               //save this new m as the initial value, so it can be saved and used in the final equation.
               for (int i = my_env_start_index; i < my_env_end_index; i++){
                  int cell = env::none_atomistic_cells[i];
                  x_array[cell] = env::x_mag_array[cell]*1.0/env::Ms[cell];
                  y_array[cell] = env::y_mag_array[cell]*1.0/env::Ms[cell];
                  z_array[cell] = env::z_mag_array[cell]*1.0/env::Ms[cell];
                  x_initial_spin_array[cell] = x_array[cell];
                  y_initial_spin_array[cell] = y_array[cell];
                  z_initial_spin_array[cell] = z_array[cell];
               }


               const double kB = 1.3806503e-23;
               std::vector<double> m(3,0.0);
               std::vector<double> spin_field(3,0.0);

               env::ext_field[0] = H*Hx;
               env::ext_field[1] = H*Hy;
               env::ext_field[2] = H*Hz;


               //calculte chi as a function of temperature
               env::one_o_chi_para =  env::calculate_chi_para(temperature);
               env::one_o_chi_perp =  env::calculate_chi_perp(temperature);


               //fill the noise terms
               for (int i = my_env_start_index; i < my_env_end_index; i++){
                  int cell = env::none_atomistic_cells[i];
                  GW1x[cell] = mtrandom::gaussian();
                  GW1y[cell] = mtrandom::gaussian();
                  GW1z[cell] = mtrandom::gaussian();
                  GW2x[cell] = mtrandom::gaussian();
                  GW2y[cell] = mtrandom::gaussian();
                  GW2z[cell] = mtrandom::gaussian();
               }
               //iff FFt is enabled calculate the demag fields.
            //   #ifdef FFT
               if (sim::time %demag_update_rate == 0) env::calculate_demag_fields();
            //   #endif
               //std::vector < double > mm_env_exchange(num_cells,0.0);
               //std::vector < double > env_mm_exchange(mm::num_cells,0.0);

               for (int i = 0; i < environment::num_interactions; i ++ ){
                  int mm_cell = environment::list_of_mm_cells_with_neighbours[i];
                  int env_cell = environment::list_of_env_cells_with_neighbours[i];
                  double overlap = list_of_overlap_area[i];

               }
               if (sim::time%1000) env::output();

         return 0;

      }


   }
