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

// file scope variables (for performance)
std::vector<double> spin_field(3,0.0);
std::vector<double> m(3,0.0);

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

   std::vector <double> x_total_external_field_array;
   std::vector <double> y_total_external_field_array;
   std::vector <double> z_total_external_field_array;

   std::vector <double> x_total_spin_field_array;
   std::vector <double> y_total_spin_field_array;
   std::vector <double> z_total_spin_field_array;

   bool LLG_set=false; ///< Flag to define state of LLG arrays (initialised/uninitialised)

}

namespace micromagnetic{

   int micromagnetic_init_llg(int num_cells){

      // check calling of routine if error checking is activated
      if(err::check==true) std::cout << "LLB_init has been called" << std::endl;

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

      x_total_external_field_array.resize(num_cells,0.0);
      y_total_external_field_array.resize(num_cells,0.0);
      z_total_external_field_array.resize(num_cells,0.0);

      x_total_spin_field_array.resize(num_cells,0.0);
      y_total_spin_field_array.resize(num_cells,0.0);
      z_total_spin_field_array.resize(num_cells,0.0);

      LLG_set=true;

      return EXIT_SUCCESS;
   }


   int LLG( std::vector <int> &local_cell_array,
            int num_steps,
            int num_cells,                     // number of cells (all processors)
            int num_local_cells,               // number of cells (local processor)
            double temperature,                // system temperature
            std::vector<double>& x_mag_array,  // un-normalised moment array (mu_B)
            std::vector<double>& y_mag_array,
            std::vector<double>& z_mag_array,
            double Hx, double Hy, double Hz,   // External applied direction
            double H,                          // External applied field (Tesla)
            double dt,                         // timestep (gamma*dt)
            std::vector <double> &volume_array  // volume of each magnetic_cell
   ){

      using namespace micromagnetic_arrays_llg;

      // Check for initialisation of LLG integration arrays
      if(LLG_set== false) micromagnetic::micromagnetic_init_llg(num_cells);

      // set cell temperatures
      for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
         int cell = list_of_micromagnetic_cells[lc];
         mm::T[cell] = temperature;
      }

      // calculate chi(T).
      mm::calculate_chi_perp(number_of_micromagnetic_cells, list_of_micromagnetic_cells, mm::one_o_chi_perp, mm::T, mm::Tc);

      // Set the current external field
      mm::ext_field[0] = H*Hx;
      mm::ext_field[1] = H*Hy;
      mm::ext_field[2] = H*Hz;

      // The LLG assumes unit vectors - normalises the x,y,z magnetisations to have a length of 1
      for (int cell = 0; cell < num_cells; cell++){

         double mx = x_mag_array[cell];
         double my = y_mag_array[cell];
         double mz = z_mag_array[cell];
         const double im_squared = 1.0/sqrt(mx*mx + my*my +mz*mz);
         mx = mx*im_squared;
         my = my*im_squared;
         mz = mz*im_squared;

         x_array[cell] = mx; // normalised magnetizations
         y_array[cell] = my;
         z_array[cell] = mz;

         x_initial_spin_array[cell] = mx; // normalised magnetizations
         y_initial_spin_array[cell] = my;
         z_initial_spin_array[cell] = mz;
        // std::cout << "start " << cell << "\t"<< mx << "\t" << my << "\t" << mz << std::endl;
      }

      // Calculate spin dependent and external fields
      mm::calculate_llg_spin_fields    (temperature, num_cells, x_array,y_array,z_array, x_total_spin_field_array, y_total_spin_field_array, z_total_spin_field_array);
      mm::calculate_llg_external_fields(temperature, num_cells, x_array,y_array,z_array, x_total_external_field_array, y_total_external_field_array, z_total_external_field_array);

      //---------------------------------------------------------------------------
      // calculates the euler step
      //---------------------------------------------------------------------------
      for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){

         int cell = list_of_micromagnetic_cells[lc];

         //calculates 1/(1+a^2) and a/(1+a^2) for llg
         const double one_oneplusalpha_sq = -1.0/(1.0+mm::alpha_perp[cell]*mm::alpha_perp[cell]); // material specific alpha and gamma
         const double alpha_oneplusalpha_sq = mm::alpha_perp[cell] * one_oneplusalpha_sq;

         const double S[3] = {x_array[cell], y_array[cell], z_array[cell]};
         const double H[3] = {x_total_spin_field_array[cell] + x_total_external_field_array[cell],
                              y_total_spin_field_array[cell] + y_total_external_field_array[cell],
                              z_total_spin_field_array[cell] + z_total_external_field_array[cell]};

         // calculates the delta S stores in euler array
         x_euler_array[cell] = (one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
         y_euler_array[cell] = (one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
         z_euler_array[cell] = (one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

      }

      // For parallel version set arrays to large negative number every time
      // to allow parallel reduction to work (MPI_MAX always picks in positive order)
      #ifdef MPICF
         for(int cell=0; cell< x_spin_storage_array.size(); cell++){
            x_spin_storage_array[cell] = 0.0;
            y_spin_storage_array[cell] = 0.0;
            z_spin_storage_array[cell] = 0.0;
         }
      #endif

      //saves the new S array from euler step and normalises
      for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){

         int cell = list_of_micromagnetic_cells[lc];

         const double mx_new = x_array[cell] + x_euler_array[cell]*dt;
         const double my_new = y_array[cell] + y_euler_array[cell]*dt;
         const double mz_new = z_array[cell] + z_euler_array[cell]*dt;

         // Normalise Spin Length
         const double imm = 1.0/sqrt(mx_new*mx_new + my_new*my_new + mz_new*mz_new);

         x_spin_storage_array[cell] = mx_new * imm;
         y_spin_storage_array[cell] = my_new * imm;
         z_spin_storage_array[cell] = mz_new * imm;

      }

      // Reduce cell magnetizations on all processors to enable correct exchange field calculations
      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, &x_spin_storage_array[0], num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &y_spin_storage_array[0], num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &z_spin_storage_array[0], num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif

      mm::calculate_llg_spin_fields(temperature, num_cells, x_spin_storage_array,     y_spin_storage_array,     z_spin_storage_array,
                                                            x_total_spin_field_array, y_total_spin_field_array, z_total_spin_field_array);

      //---------------------------------------------------------------------------
      // Heun step
      //---------------------------------------------------------------------------
      for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){

         int cell = list_of_micromagnetic_cells[lc];

         const double one_oneplusalpha_sq = -1.0/(1.0+mm::alpha_perp[cell]*mm::alpha_perp[cell]); // material specific alpha and gamma
         const double alpha_oneplusalpha_sq = mm::alpha_perp[cell] * one_oneplusalpha_sq;

         const double S[3] = {x_spin_storage_array[cell], y_spin_storage_array[cell], z_spin_storage_array[cell]};
         const double H[3] = {x_total_spin_field_array[cell] + x_total_external_field_array[cell],
                              y_total_spin_field_array[cell] + y_total_external_field_array[cell],
                              z_total_spin_field_array[cell] + z_total_external_field_array[cell]};

         //saves delta to xyz
         x_heun_array[cell] = (one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
         y_heun_array[cell] = (one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
         z_heun_array[cell] = (one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

      }

      // For parallel version set arrays to large negative number every time
      // to allow parallel reduction to work
      #ifdef MPICF
         for(int cell=0; cell< x_array.size(); cell++){
            x_array[cell] = 0.0;
            y_array[cell] = 0.0;
            z_array[cell] = 0.0;
            cells::mag_array_x[cell] = 0.0;
            cells::mag_array_y[cell] = 0.0;
            cells::mag_array_z[cell] = 0.0;
         }
      #endif

      const double hdt = 0.5*dt;

      //saves the new S array from euler step and normalises
      for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){

         int cell = list_of_micromagnetic_cells[lc];

         const double mx_new = x_array[cell] + x_euler_array[cell]*dt;
         const double my_new = y_array[cell] + y_euler_array[cell]*dt;
         const double mz_new = z_array[cell] + z_euler_array[cell]*dt;

         // Normalise Spin Length
         const double imm = 1.0/sqrt(mx_new*mx_new + my_new*my_new + mz_new*mz_new);

         x_spin_storage_array[cell] = mx_new * imm;
         y_spin_storage_array[cell] = my_new * imm;
         z_spin_storage_array[cell] = mz_new * imm;

      }

      //calculates new spin arrays from heun and euler steps
      for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
         int cell = list_of_micromagnetic_cells[lc];

         const double mx_new = x_initial_spin_array[cell]+(hdt)*(x_euler_array[cell]+x_heun_array[cell]);
         const double my_new = y_initial_spin_array[cell]+(hdt)*(y_euler_array[cell]+y_heun_array[cell]);
         const double mz_new = z_initial_spin_array[cell]+(hdt)*(z_euler_array[cell]+z_heun_array[cell]);

         // Normalise Spin Length
         const double ms = mm::ms[cell];
         const double imm = 1.0/sqrt(mx_new*mx_new + my_new*my_new + mz_new*mz_new);

         x_array[cell] = mx_new * imm;
         y_array[cell] = my_new * imm;
         z_array[cell] = mz_new * imm;

         cells::mag_array_x[cell] = x_array[cell]*ms;
         cells::mag_array_y[cell] = y_array[cell]*ms;
         cells::mag_array_z[cell] = z_array[cell]*ms;
     //    std::cout << "end " << cell << "\t"<< cells::mag_array_x[cell]/ms  << "\t" << cells::mag_array_y[cell]/ms  << "\t" << cells::mag_array_z[cell]/ms  << std::endl;

      }


      // Reduce unit vectors and moments to all processors
      #ifdef MPICF
      	MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_x[0], num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      	MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_y[0], num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      	MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_z[0], num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      	MPI_Allreduce(MPI_IN_PLACE, &x_array[0],            num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      	MPI_Allreduce(MPI_IN_PLACE, &y_array[0],            num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      	MPI_Allreduce(MPI_IN_PLACE, &z_array[0],            num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif

     // updates atom magnetisations
       if (discretisation_type  ==2){
         for(int atom_list=0;atom_list<number_of_none_atomistic_atoms;atom_list++){
            int atom = list_of_none_atomistic_atoms[atom_list];
            int cell = cells::atom_cell_id_array[atom];
            double me = mm::m_e[cell];
            atoms::x_spin_array[atom] = x_array[cell]*me;
            atoms::y_spin_array[atom] = y_array[cell]*me;
            atoms::z_spin_array[atom] = z_array[cell]*me;
            atoms::m_spin_array[atom] = me;
         }
      }

      // if (enable_resistance && mm::resistance_layer_2 != mm::resistance_layer_1)  {
      //    micromagnetic::MR_resistance = mm::calculate_resistance();
      //   // std::cout << micromagnetic::MR_resistance  << std::endl;
      // }
      // if (sim::time%sim::partial_time == 0) {

      //    //std::cout << sim::time%sim::partial_time << "\t" << sim::time << "\t" << sim::partial_time << std::endl;
      //    mm::outputs();
      // }

      return 0;

   }

} // end of micromagnetic namespace
