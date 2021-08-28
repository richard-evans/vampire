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
                //std::cout << "LL" << "\t" <<env::num_cells << "\t" << my_env_start_index << '\t' << my_env_end_index << std::endl;

               // Check for initialisation of LLG integration arrays
               if(LLG_set== false) environment::LLB_init(num_cells);
               // Local variables for system integration

//std::cout << "HEREa" << "\t" << env::num_cells <<std::endl;

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
               //                     std::cout << "HEREb" << "\t" << env::x_mag_array[cell] << '\t' << env::Ms[cell] << "\t" << x_array[cell] << std::endl;
               }


               //const double kB = 1.3806503e-23;
               std::vector<double> m(3,0.0);
               std::vector<double> spin_field(3,0.0);

               env::ext_field[0] = H*Hx;
               env::ext_field[1] = H*Hy;
               env::ext_field[2] = H*Hz;



            //   std::cout << "HERE" <<"\t" << my_env_start_index << '\t' << my_env_end_index << std::endl;
               //fill the noise terms
               for (int i = my_env_start_index; i < my_env_end_index; i++){
                  int cell = env::none_atomistic_cells[i];
               //   std::cout << cell << '\t' <<  env::none_atomistic_cells.size()<< std::endl;
                  //calculte chi as a function of temperature
                  env::one_o_chi_para[cell] =  env::calculate_chi_para(temperature,cell);
                  env::one_o_chi_perp[cell] =  env::calculate_chi_perp(temperature,cell);
                  GW1x[cell] = mtrandom::gaussian();
                  GW1y[cell] = mtrandom::gaussian();
                  GW1z[cell] = mtrandom::gaussian();
                  GW2x[cell] = mtrandom::gaussian();
                  GW2y[cell] = mtrandom::gaussian();
                  GW2z[cell] = mtrandom::gaussian();
               }

            //   std::cout << "HERE2" << std::endl;

               //iff FFt is enabled calculate the demag fields.
              // #ifdef FFT
              // if (sim::time %demag_update_rate == 0) env::calculate_demag_fields();
              // #endif
            if (sim::time %demag_update_rate == 0)   env::calculate_demag_fields();
               //std::vector < double > mm_env_exchange(num_cells,0.0);
               //std::vector < double > env_mm_exchange(mm::num_cells,0.0);

               //for (int i = 0; i < environment::num_interactions; i ++ ){
                  //int mm_cell = environment::list_of_mm_cells_with_neighbours[i];
                  //int env_cell = environment::list_of_env_cells_with_neighbours[i];
                  //double overlap = list_of_overlap_area[i];

               //}

               const double kB = 1.3806503e-23;

//std::cout << "HERE3" <<std::endl;
               for (int i = my_env_start_index; i < my_env_end_index; i++){
                  int cell = env::none_atomistic_cells[i];
                  m[0] = x_array[cell];
                  m[1] = y_array[cell];
                  m[2] = z_array[cell];
                  //std::cout << "HERE3" << "\t" << m[0] << '\t' << m[1] << "\t" << m[2] << std::endl;
                  spin_field = env::calculate_llb_fields(m, temperature, cell, x_array,y_array,z_array);
                //  std::cout << "field\t" << spin_field[0] << '\t' << spin_field[1] << '\t' << spin_field[2] <<std::endl;
                  //calcualte the noise terms
                  double sigma_para = sqrt(2*kB*temperature*env::alpha_para/(env::Ms[cell]*dt));
                  double sigma_perp = sqrt(2*kB*temperature*(env::alpha_perp-env::alpha_para)/(dt*env::Ms[cell]*env::alpha_perp*env::alpha_perp));
                  const double H[3] = {spin_field[0], spin_field[1], spin_field[2]};
                  //std::cout << spin_field[0] << std::endl;
                  //saves the noise terms to an array for easy access
                  const double GW2t[3] = {GW2x[cell],GW2y[cell],GW2z[cell]};
                  const double one_o_m_squared = 1.0/(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
                  const double SdotH = m[0]*H[0] + m[1]*H[1] + m[2]*H[2];
                  double xyz[3] = {0.0,0.0,0.0};
                  //calculates the LLB equation
                  xyz[0]=  - (m[1]*H[2]-m[2]*H[1])
                  + env::alpha_para*m[0]*SdotH*one_o_m_squared
                  - env::alpha_perp*(m[1]*(m[0]*H[1]-m[1]*H[0])-m[2]*(m[2]*H[0]-m[0]*H[2]))*one_o_m_squared
                  + GW1x[cell]*sigma_para
                  - env::alpha_perp*(m[1]*(m[0]*GW2t[1]-m[1]*GW2t[0])-m[2]*(m[2]*GW2t[0]-m[0]*GW2t[2]))*one_o_m_squared*sigma_perp;

                  xyz[1]=  - (m[2]*H[0]-m[0]*H[2])
                  + env::alpha_para*m[1]*SdotH*one_o_m_squared
                  - env::alpha_perp*(m[2]*(m[1]*H[2]-m[2]*H[1])-m[0]*(m[0]*H[1]-m[1]*H[0]))*one_o_m_squared
                  + GW1y[cell]*sigma_para
                  - env::alpha_perp*(m[2]*(m[1]*GW2t[2]-m[2]*GW2t[1])-m[0]*(m[0]*GW2t[1]-m[1]*GW2t[0]))*one_o_m_squared*sigma_perp;

                  xyz[2]=	 - (m[0]*H[1]-m[1]*H[0])
                  + env::alpha_para*m[2]*SdotH*one_o_m_squared
                  - env::alpha_perp*(m[0]*(m[2]*H[0]-m[0]*H[2])-m[1]*(m[1]*H[2]-m[2]*H[1]))*one_o_m_squared
                  + GW1z[cell]*sigma_para
                  - env::alpha_perp*(m[0]*(m[2]*GW2t[0]-m[0]*GW2t[2])-m[1]*(m[1]*GW2t[2]-m[2]*GW2t[1]))*one_o_m_squared*sigma_perp;

                  x_euler_array[cell] = xyz[0];
                  y_euler_array[cell] = xyz[1];
                  z_euler_array[cell] = xyz[2];
                //  std::cout << "xyz\t" << xyz[0] << '\t' << xyz[1] << '\t' << xyz[2] << "\t" << m[0] << '\t' << m[1] << '\t' << m[2] <<"\t" <<  SdotH << "\t" << one_o_m_squared << '\t' << env::alpha_para << '\t' << env::alpha_perp << std::endl;

               }
      //         std::cout << "HERE4" << std::endl;
        //    std::cin.get();
               //save the euler step to a spin storage array
               for (int i = my_env_start_index; i < my_env_end_index; i++){
                  int cell = env::none_atomistic_cells[i];
                  x_spin_storage_array[cell] = x_array[cell] + x_euler_array[cell]*dt;
                  y_spin_storage_array[cell] = y_array[cell] + y_euler_array[cell]*dt;
                  z_spin_storage_array[cell] = z_array[cell] + z_euler_array[cell]*dt;
          //        std::cout << "storage" << '\t' << x_spin_storage_array[cell] << '\t' << y_spin_storage_array[cell] << '\t' << z_spin_storage_array[cell] << std::endl;

               }

        //    std::cin.get();
//                    for (int i =0; i <num_env_cells; i++){
  //        if (vmpi::my_rank == 1)     std::cerr << i << "\t" << x_spin_storage_array[i] << "\t" << my_env_start_index << '\t' << my_env_end_index << "\t" << num_env_cells << std::endl;
    //         }
            //   std::cout << x_spin_storage_array.size() <<  "\t" << num_cells << std::endl;
                #ifdef MPICF
              MPI_Allreduce(MPI_IN_PLACE, &x_spin_storage_array[0],     num_env_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &y_spin_storage_array[0],     num_env_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &z_spin_storage_array[0],     num_env_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
                #endif

         //     std::cout << "HERE5" << std::endl;


               for (int i = my_env_start_index; i < my_env_end_index; i++){
                  int cell = env::none_atomistic_cells[i];

                  m[0] = x_spin_storage_array[cell];
                  m[1] = y_spin_storage_array[cell];
                  m[2] = z_spin_storage_array[cell];

            //      std::cout << "storage" << '\t' << x_spin_storage_array[cell] << '\t' << y_spin_storage_array[cell] << '\t' << z_spin_storage_array[cell] << std::endl;


                  //calcualte spin fields
                  spin_field = env::calculate_llb_fields(m, temperature, cell, x_spin_storage_array, y_spin_storage_array, z_spin_storage_array);
            //      std::cout << "field\t" << spin_field[0] << '\t' << spin_field[1] << '\t' << spin_field[2] <<std::endl;

                  double sigma_para = sqrt(2.0*kB*temperature*env::alpha_para/(env::Ms[cell]*dt)); //why 1e-27
                  double sigma_perp = sqrt(2.0*kB*temperature*(env::alpha_perp-env::alpha_para)/(dt*env::Ms[cell]*env::alpha_perp*env::alpha_perp));
                  const double H[3] = {spin_field[0], spin_field[1], spin_field[2]};

                  //saves the noise terms to an array
                  const double GW2t[3] = {GW2x[cell],GW2y[cell],GW2z[cell]};
                  const double one_o_m_squared = 1.0/(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
                  const double SdotH = m[0]*H[0] + m[1]*H[1] + m[2]*H[2];
                  double xyz[3] = {0.0,0.0,0.0};
                  //calculates the LLB equation
                  xyz[0]=  - (m[1]*H[2]-m[2]*H[1])
                  + env::alpha_para*m[0]*SdotH*one_o_m_squared
                  - env::alpha_perp*(m[1]*(m[0]*H[1]-m[1]*H[0])-m[2]*(m[2]*H[0]-m[0]*H[2]))*one_o_m_squared
                  + GW1x[cell]*sigma_para
                  - env::alpha_perp*(m[1]*(m[0]*GW2t[1]-m[1]*GW2t[0])-m[2]*(m[2]*GW2t[0]-m[0]*GW2t[2]))*one_o_m_squared*sigma_perp;

                  xyz[1]=  - (m[2]*H[0]-m[0]*H[2])
                  + env::alpha_para*m[1]*SdotH*one_o_m_squared
                  - env::alpha_perp*(m[2]*(m[1]*H[2]-m[2]*H[1])-m[0]*(m[0]*H[1]-m[1]*H[0]))*one_o_m_squared
                  + GW1y[cell]*sigma_para
                  - env::alpha_perp*(m[2]*(m[1]*GW2t[2]-m[2]*GW2t[1])-m[0]*(m[0]*GW2t[1]-m[1]*GW2t[0]))*one_o_m_squared*sigma_perp;

                  xyz[2]=	 - (m[0]*H[1]-m[1]*H[0])
                  + env::alpha_para*m[2]*SdotH*one_o_m_squared
                  - env::alpha_perp*(m[0]*(m[2]*H[0]-m[0]*H[2])-m[1]*(m[1]*H[2]-m[2]*H[1]))*one_o_m_squared
                  + GW1z[cell]*sigma_para
                  - env::alpha_perp*(m[0]*(m[2]*GW2t[0]-m[0]*GW2t[2])-m[1]*(m[1]*GW2t[2]-m[2]*GW2t[1]))*one_o_m_squared*sigma_perp;

                  //heun delta = xyz
                  x_heun_array[cell] = xyz[0];
                  y_heun_array[cell] = xyz[1];
                  z_heun_array[cell] = xyz[2];
        //          std::cout << "xyz\t" << xyz[0] << '\t' << xyz[1] << '\t' << xyz[2] << "\t" << m[0] << '\t' << m[1] << '\t' << m[2] <<"\t" <<  SdotH << "\t" << one_o_m_squared << '\t' << env::alpha_para << '\t' << env::alpha_perp << std::endl;

               }

          //  std::cin.get();

               for (int i = my_env_start_index; i < my_env_end_index; i++){
                  int cell = env::none_atomistic_cells[i];
                  env::x_mag_array[cell] = 0.0;
                  env::y_mag_array[cell] = 0.0;
                  env::z_mag_array[cell] = 0.0;
               }


               //calcualtes the new magnetisations from the heun and euler deltas
               for (int i = my_env_start_index; i < my_env_end_index; i++){
                  int cell = env::none_atomistic_cells[i];

                  x_array[cell] = x_initial_spin_array[cell] + 0.5*dt*(x_euler_array[cell] + x_heun_array[cell]);
                  y_array[cell] = y_initial_spin_array[cell] + 0.5*dt*(y_euler_array[cell] + y_heun_array[cell]);
                  z_array[cell] = z_initial_spin_array[cell] + 0.5*dt*(z_euler_array[cell] + z_heun_array[cell]);

                  env::x_mag_array[cell] = x_array[cell]*env::Ms[cell];
                  env::y_mag_array[cell] = y_array[cell]*env::Ms[cell];
                  env::z_mag_array[cell] = z_array[cell]*env::Ms[cell];
          //        std::cout << env::x_mag_array[cell] << '\t' <<env::y_mag_array[cell] << '\t' << env::z_mag_array[cell] <<std::endl;
               }

          //     std::cin.get();

//
               // #ifdef MPICF
               // MPI_Allreduce(MPI_IN_PLACE, &env::x_mag_array[0],     num_env_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
               // MPI_Allreduce(MPI_IN_PLACE, &env::y_mag_array[0],     num_env_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
               // MPI_Allreduce(MPI_IN_PLACE, &env::z_mag_array[0],     num_env_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
               // #endif

               //for (int i = my_env_start_index; i < my_env_end_index; i++){
                  //int cell = env::none_atomistic_cells[i];

              //    std::cout << env::x_mag_array[cell] << '\t' <<env::y_mag_array[cell] << '\t' << env::z_mag_array[cell] <<std::endl;
               //}
          //     std::cin.get();
//
// //               outputs magnetisation
           if (sim::time %(vout::output_rate*1000) == 0 && vmpi::my_rank == 0 ) env::output();
//

         return 0;

      }


   }
