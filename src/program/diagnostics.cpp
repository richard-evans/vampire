//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains the standard benchmark program
///
/// @details Simulates a system for a number of timesteps
///
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    09/03/2011
/// @internal
///   Created:    05/02/2011
///   Revision:   09/03/2011
///=====================================================================================
///

// Standard Libraries
#include "constants.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>

// Vampire Header files
#include "anisotropy.hpp"
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"


namespace program{

   int timestep_scaling(){

      // check calling of routine if error checking is activated
      if(err::check==true){std::cout << "program::timestep_scaling has been called" << std::endl;}

      std::cout << " Diagnostic - Timestep Scaling " << std::endl;

      // loop over timesteps
      for(int powerv=18; powerv > 13; powerv--){
         for(int value=1;value<10;value++){
            mp::dt_SI=double(value)*pow(10.0,-1.0*powerv);

            uint64_t timesteps = 5.0e-12/mp::dt_SI;

            std::cout << timesteps << std::endl;

            // reset derived parameters
            mp::set_derived_parameters();

            double sx = 0.01;
            double sy = 0.0;
            double sz = 1.0;

            double modS = 1.0/sqrt(sx*sx + sy*sy + sz*sz);

            sx*=modS;
            sy*=modS;
            sz*=modS;

            for(int atom=0; atom<atoms::num_atoms; atom++){
               atoms::x_spin_array[atom] = sx;
               atoms::y_spin_array[atom] = sy;
               atoms::z_spin_array[atom] = sz;
            }

            sim::integrate(timesteps);
            stats::reset();
            uint64_t start_time = sim::time;
            // Simulate system
            while( sim::time < timesteps+start_time ){
               sim::integrate(1);

               // Calculate mag_m, mag after sim::partial_time steps
               stats::update();

            } // end of time loop
            zmag << mp::dt_SI << "\t";
            std::cout << mp::dt_SI << "\t";
            vout::data();
         }
      }

      return EXIT_SUCCESS;
   }

   void boltzmann_dist(){

      // check calling of routine if error checking is activated
      if(err::check==true) std::cout << "program::boltzmann_dist has been called" << std::endl;

      // array for binning spin angle
      std::vector<double> bin(181,0.0);

      // Equilibrate system
      sim::integrate(sim::equilibration_time);

      // Simulate system
      while(sim::time<sim::total_time+sim::equilibration_time){

         // Integrate system
         sim::integrate(sim::partial_time);

         // Calculate magnetisation statistics
         for(int atom=0; atom<atoms::num_atoms; atom++){
            double angle = acos(atoms::z_spin_array[atom])*180.0/M_PI;
            double id = vmath::iround(angle+0.5);
            bin[id]+=1.0;
         }
      }

      // Find max probability and max P
      double maxPK = 0.0;
      double maxPH = 0.0;
      double maxP  = 0.0;
      for(int b=0;b<181;b++){
         double energyK = anisotropy::get_anisotropy_constant(0)/(sim::temperature*constants::kB); // get anisotropy constant for material 0
         double energyH = sim::H_applied*mp::material[0].mu_s_SI/(sim::temperature*constants::kB); // get anisotropy constant for material 0
         double PK = sin(double (b)*M_PI/180.0)*exp(energyK*sin(double (b)*M_PI/180.0)*sin(double (b)*M_PI/180.0));
         double PH = sin(double (b)*M_PI/180.0)*exp(energyH*cos(double (b)*M_PI/180.0));
         if((bin[b])>maxP) maxP=bin[b];
         if(PK>maxPK) maxPK=PK;
         if(PH>maxPH) maxPH=PH;
      }

      std::ofstream ofile("boltzmann-distribution.txt");

      // Output data
      ofile << "# Anisotropy:          " << anisotropy::get_anisotropy_constant(0) << "\t" << anisotropy::get_anisotropy_constant(0)/(sim::temperature*constants::kB) << std::endl;
      ofile << "# Field:               " << sim::H_applied << "\t" << sim::H_applied*mp::material[0].mu_s_SI/(sim::temperature*constants::kB) << std::endl;
      ofile << "# Moment, Temperature: " << mp::material[0].mu_s_SI/constants::muB << "\t" << sim::temperature << std::endl;
      for(int b=0;b<181;b++){
         double energyK = anisotropy::get_anisotropy_constant(0)/(sim::temperature*constants::kB); // get anisotropy constant for material 0
         double energyH = sim::H_applied*mp::material[0].mu_s_SI/(sim::temperature*constants::kB); // get anisotropy constant for material 0
         double PK = sin(double (b)*M_PI/180.0)*exp(energyK*sin(double (b)*M_PI/180.0)*sin(double (b)*M_PI/180.0));
         double PH = sin(double (b)*M_PI/180.0)*exp(energyH*cos(double (b)*M_PI/180.0));
         ofile << b << "\t" << bin[b]/maxP << "\t" << (bin[b]+bin[180-b])/(2.0*maxP) << "\t" << PK/maxPK << "\t" << PH/maxPH << std::endl;
      }

      // close Boltzmann file
      ofile.close();

   }

   void boltzmann_dist_micromagnetic_llg(){

      // check calling of routine if error checking is activated
      if(err::check==true) std::cout << "program::boltzmann_dist has been called" << std::endl;

      // array for binning spin angle
      std::vector<double> bin(181,0.0);

      // Equilibrate system
      sim::integrate(sim::equilibration_time);

      // Simulate system
      while(sim::time<sim::total_time+sim::equilibration_time){

         // Integrate system
         sim::integrate(sim::partial_time);

         // Calculate magnetisation statistics
         for(int atom=0; atom<atoms::num_atoms; atom++){
            const double mx = atoms::x_spin_array[atom];
            const double my = atoms::y_spin_array[atom];
            const double mz = atoms::z_spin_array[atom];
            const double mm = sqrt(mx*mx + my*my + mz*mz);
            double angle = acos(mz/mm)*180.0/M_PI;
            double id = vmath::iround(angle+0.5);
            //std::cout << id << "\t" << angle << "\t" << mx << "\t" << my << "\t" << mz << "\t" << mm << std::endl;
            bin[id]+=1.0;
         }

         // Calculate magnetisation statistics
         stats::update();

         // Output data
         vout::data();

      }

      // Find max probability and max P
      double maxPK = 0.0;
      double maxPH = 0.0;
      double maxP  = 0.0;
      for(int b=0;b<181;b++){
         double energyK = anisotropy::get_anisotropy_constant(0)/(sim::temperature*constants::kB); // get anisotropy constant for material 0
         double energyH = sim::H_applied*mp::material[0].mu_s_SI/(sim::temperature*constants::kB); // get anisotropy constant for material 0
         double PK = sin(double (b)*M_PI/180.0)*exp(energyK*sin(double (b)*M_PI/180.0)*sin(double (b)*M_PI/180.0));
         double PH = sin(double (b)*M_PI/180.0)*exp(energyH*cos(double (b)*M_PI/180.0));
         if((bin[b])>maxP) maxP=bin[b];
         if(PK>maxPK) maxPK=PK;
         if(PH>maxPH) maxPH=PH;
      }

      std::ofstream ofile("boltzmann-distribution.txt");

      // Output data
      ofile << "# Anisotropy:          " << anisotropy::get_anisotropy_constant(0) << "\t" << anisotropy::get_anisotropy_constant(0)/(sim::temperature*constants::kB) << std::endl;
      ofile << "# Field:               " << sim::H_applied << "\t" << sim::H_applied*mp::material[0].mu_s_SI/(sim::temperature*constants::kB) << std::endl;
      ofile << "# Moment, Temperature: " << mp::material[0].mu_s_SI/constants::muB << "\t" << sim::temperature << std::endl;
      for(int b=0;b<181;b++){
         double energyK = anisotropy::get_anisotropy_constant(0)/(sim::temperature*constants::kB); // get anisotropy constant for material 0
         double energyH = sim::H_applied*mp::material[0].mu_s_SI/(sim::temperature*constants::kB); // get anisotropy constant for material 0
         double PK = sin(double (b)*M_PI/180.0)*exp(energyK*sin(double (b)*M_PI/180.0)*sin(double (b)*M_PI/180.0));
         double PH = sin(double (b)*M_PI/180.0)*exp(energyH*cos(double (b)*M_PI/180.0));
         ofile << b << "\t" << bin[b]/maxP << "\t" << (bin[b]+bin[180-b])/(2.0*maxP) << "\t" << PK/maxPK << "\t" << PH/maxPH << std::endl;
      }

      // close Boltzmann file
      ofile.close();

   }

}//end of namespace program
