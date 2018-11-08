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

#include "random.hpp"
#include "errors.hpp"
#include "atoms.hpp"
#include "cells.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>

namespace mm = micromagnetic::internal;

namespace micromagnetic_arrays_llg{

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

   bool LLG_set=false; ///< Flag to define state of LLG arrays (initialised/uninitialised)

}


namespace micromagnetic{

   int micromagnetic_init_llg(int num_cells){

      // check calling of routine if error checking is activated
      if(err::check==true) std::cout << "LLG_init has been called" << std::endl;

      using namespace micromagnetic_arrays_llg;

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

      LLG_set=true;

      return EXIT_SUCCESS;
   }


   int LLG( std::vector <int> local_cell_array,
      int num_steps,
      int num_cells,
      int num_local_cells,
      double temperature,
      std::vector<double>& x_mag_array,
      std::vector<double>& y_mag_array,
      std::vector<double>& z_mag_array,
      double Hx,
      double Hy,
      double Hz,
      double H,
      double dt,
      std::vector <double> volume_array
   ){

      // check calling of routine if error checking is activated
      if(err::check==true){std::cout << "micromagnetic::LLG_Heun has been called" << std::endl;}

      using namespace micromagnetic_arrays_llg;

      // Check for initialisation of LLG integration arrays
      if(LLG_set== false) micromagnetic::micromagnetic_init_llg(num_cells);
      // Local variables for system integration

      //calculte chi(T).
      mm::one_o_chi_para =  mm::calculate_chi_para(num_local_cells,local_cell_array, num_cells, temperature);
      mm::one_o_chi_perp =  mm::calculate_chi_perp(num_local_cells,local_cell_array, num_cells, temperature);

      //The external fields equal the length of the field times the applied field vector.
      //This is saved to an array.
      mm::ext_field[0] = H*Hx;
      mm::ext_field[1] = H*Hy;
      mm::ext_field[2] = H*Hz;

      // const double real_time = sim::time*mp::dt_SI;
   	// const double omega = sim::fmr_field_frequency*1.e9; // Hz
   	// const double Hfmrx = sim::fmr_field_unit_vector[0];
   	// const double Hfmry = sim::fmr_field_unit_vector[1];
   	// const double Hfmrz = sim::fmr_field_unit_vector[2];
   	// const double Hsinwt = sim::fmr_field_strength * sin(2.0 * M_PI * omega * real_time);
   	// mm::fmr_H[0] = Hfmrx * Hsinwt;
   	// mm::fmr_H[1] = Hfmry * Hsinwt;
   	// mm::fmr_H[2] = Hfmrz * Hsinwt;


      //save this new m as the initial value, so it can be saved and used in the final equation.
      //normalises the x,y,z magnetisations to have a length of 1
      for (int lc = 0; lc < num_local_cells; lc++){
         int cell = local_cell_array[lc];
         x_array[cell] = x_mag_array[cell];
         y_array[cell] = y_mag_array[cell];
         z_array[cell] = z_mag_array[cell];
         double m_squared = sqrt(x_array[cell]*x_array[cell] + y_array[cell]*y_array[cell] + z_array[cell]*z_array[cell]);
         x_array[cell] = x_array[cell]/m_squared;
         y_array[cell] = y_array[cell]/m_squared;
         z_array[cell] = z_array[cell]/m_squared;
         x_initial_spin_array[cell] = x_array[cell];
         y_initial_spin_array[cell] = y_array[cell];
         z_initial_spin_array[cell] = z_array[cell];
      }

      std::vector<double> spin_field(3,0.0);
      std::vector<double> m(3,0.0);

      //calcualtes the euler gradient
      for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){

         // get cell ID
         int cell = list_of_micromagnetic_cells[lc];

         // load magnetization into temporary variable
         m[0] = x_array[cell];
         m[1] = y_array[cell];
         m[2] = z_array[cell];

         // calculates spin fields
         spin_field = mm::calculate_llg_fields(m, temperature, num_cells, cell, x_array,y_array,z_array);

         //calcualtes 1/(1+a^2) and a/(1+a^2) for llg
         const double alpha = mm::alpha[cell];
         const double ialpha2 = 1.0/(1.0+alpha*alpha);
         const double one_oneplusalpha_sq = ialpha2; // material specific alpha and gamma
         const double alpha_oneplusalpha_sq = alpha*ialpha2;

         const double S[3] = {m[0],m[1],m[2]};
         const double H[3] = {-spin_field[0], -spin_field[1], -spin_field[2]};
         //calcautes the delta S stores in x,y,z
         double xyz[3] = {0.0,0.0,0.0};
         xyz[0]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
         xyz[1]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
         xyz[2]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

         // Store dS in euler array
         x_euler_array[cell] = xyz[0];
         y_euler_array[cell] = xyz[1];
         z_euler_array[cell] = xyz[2];

      }

      double S_new[3] = {0.0,0.0,0.0};
      double mod_S;

      //saves the new S array from euler step and normalises
      for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
         int cell = list_of_micromagnetic_cells[lc];

         S_new[0] = x_array[cell] + x_euler_array[cell]*dt;
         S_new[1] = y_array[cell] + y_euler_array[cell]*dt;
         S_new[2] = z_array[cell] + z_euler_array[cell]*dt;

         // Normalise Spin Length
         mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

         S_new[0]=S_new[0]*mod_S;
         S_new[1]=S_new[1]*mod_S;
         S_new[2]=S_new[2]*mod_S;

         x_spin_storage_array[cell] = S_new[0];
         y_spin_storage_array[cell] = S_new[1];
         z_spin_storage_array[cell] = S_new[2];

      }

      // heun step
      for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
         int cell = list_of_micromagnetic_cells[lc];

         m[0] = x_spin_storage_array[cell];
         m[1] = y_spin_storage_array[cell];
         m[2] = z_spin_storage_array[cell];

         spin_field = mm::calculate_llg_fields(m, temperature, num_cells, cell, x_spin_storage_array,y_spin_storage_array,z_spin_storage_array);

         //save magnetisation and field to arrays for easy access.
         const double S[3] = {m[0],m[1],m[2]};
         const double H[3] = {-spin_field[0], -spin_field[1], -spin_field[2]};

         //calcualtes 1/(1+a^2) and a/(1+a^2) for llg
         const double one_oneplusalpha_sq = 1.0/(1.0+mm::alpha[cell]*mm::alpha[cell]); // material specific alpha and gamma
         const double alpha_oneplusalpha_sq = mm::alpha[cell]/(1+mm::alpha[cell]*mm::alpha[cell]);

         //saves delta to xyz
         double xyz[3] = {0.0,0.0,0.0};
         xyz[0]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
         xyz[1]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
         xyz[2]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

         // Store dS in euler array
         x_heun_array[cell] = xyz[0];
         y_heun_array[cell] = xyz[1];
         z_heun_array[cell] = xyz[2];
      }

      //calculates new spin arrays from heun and euler steps
      for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
         int cell = list_of_micromagnetic_cells[lc];

         S_new[0]= x_initial_spin_array[cell]+(dt/2.0)*(x_euler_array[cell]+x_heun_array[cell]);
         S_new[1]= y_initial_spin_array[cell]+(dt/2)*(y_euler_array[cell]+y_heun_array[cell]);
         S_new[2]= z_initial_spin_array[cell]+(dt/2)*(z_euler_array[cell]+z_heun_array[cell]);

         // Normalise Spin Length
         mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

         x_array[cell]=S_new[0]*mod_S;
         y_array[cell]=S_new[1]*mod_S;
         z_array[cell]=S_new[2]*mod_S;

         cells::mag_array_x[cell] = x_array[cell]*mm::ms[cell];
         cells::mag_array_y[cell] = y_array[cell]*mm::ms[cell];
         cells::mag_array_z[cell] = z_array[cell]*mm::ms[cell];

      }

      //updates atom magnetisations
      if (discretisation_type  == 2 || sim::time%vout::output_rate -1){
         for(int atom_list=0;atom_list<number_of_none_atomistic_atoms;atom_list++){
            int atom = list_of_none_atomistic_atoms[atom_list];
            int cell = cells::atom_cell_id_array[atom];
            atoms::x_spin_array[atom] = x_array[cell]*mm::m_e[cell];
            atoms::y_spin_array[atom] = y_array[cell]*mm::m_e[cell];
            atoms::z_spin_array[atom] = z_array[cell]*mm::m_e[cell];
            atoms::m_spin_array[atom] = mm::m_e[cell];
         }
      }

      /*std::ofstream pfile;
      pfile.open("cell_config3");
       if(sim::time>10000){
      for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){

         int cell = list_of_micromagnetic_cells[lc];

         pfile << cell << '\t' << internal::cell_material_array[cell] << '\t'<< cells::pos_and_mom_array[4*cell+0] << '\t' <<
                                   cells::pos_and_mom_array[4*cell+1] << '\t' <<
                                   cells::pos_and_mom_array[4*cell+2] << '\t' <<
                                   cells::mag_array_x[cell] << '\t' <<
                                   cells::mag_array_y[cell] << '\t' <<
                                   cells::mag_array_z[cell] << '\t' << std::endl;

      }
      std::cin.get();
   }*/

      if (enable_resistance && mm::resistance_layer_2 != mm::resistance_layer_1)  micromagnetic::MR_resistance = mm::calculate_resistance();

      // if(sim::time>10000){
      // 	for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
      // 		int cell = list_of_micromagnetic_cells[lc];
      //
      // 		double S[3] = {x_array[cell],y_array[cell],z_array[cell]};
      // 		double mz = (S[2]+1)/2.0;
      // 		if (mz > 1.0) mz = 1.0;
      // 		if (mz < 0.0) mz = 0.0;
      // 		double mx = sqrt(S[0]*S[0]+S[1]*S[1]);
      // 		if (mx > 1.0) mx = 1.0;
      // 		double mag_m = mz;
      // 		int para = int(mz*100.0);
      // 		int perp = int(mx*100.0);
      // 		int para1D = int(mag_m*1000.0);
      // 		P[para][perp]+=1.0;
      // 		P1D[para1D]+=1.0;
      // 		mean_M+=mag_m;
      // 		counter+=1.0;
      // 	}
      // }
      return 0;

   }
}
