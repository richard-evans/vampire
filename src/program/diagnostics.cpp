//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains the standard benchmark program
///
/// @details Simulates a system for a number of timesteps
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    09/03/2011
/// @internal
///	Created:		05/02/2011
///	Revision:	09/03/2011
///=====================================================================================
///

// Standard Libraries
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
         double energyK = anisotropy::get_anisotropy_constant(0)/(sim::temperature*1.3806503e-23); // get anisotropy constant for material 0
         double energyH = sim::H_applied*mp::material[0].mu_s_SI/(sim::temperature*1.3806503e-23); // get anisotropy constant for material 0
         double PK = sin(double (b)*M_PI/180.0)*exp(energyK*sin(double (b)*M_PI/180.0)*sin(double (b)*M_PI/180.0));
         double PH = sin(double (b)*M_PI/180.0)*exp(energyH*cos(double (b)*M_PI/180.0));
         if((bin[b])>maxP) maxP=bin[b];
         if(PK>maxPK) maxPK=PK;
         if(PH>maxPH) maxPH=PH;
      }

      std::ofstream ofile("boltzmann-distribution.txt");

      // Output data
      ofile << "# Anisotropy:          " << anisotropy::get_anisotropy_constant(0) << "\t" << anisotropy::get_anisotropy_constant(0)/(sim::temperature*1.3806503e-23) << std::endl;
      ofile << "# Field:               " << sim::H_applied << "\t" << sim::H_applied*mp::material[0].mu_s_SI/(sim::temperature*1.3806503e-23) << std::endl;
      ofile << "# Moment, Temperature: " << mp::material[0].mu_s_SI/9.274e-24 << "\t" << sim::temperature << std::endl;
      for(int b=0;b<181;b++){
         double energyK = anisotropy::get_anisotropy_constant(0)/(sim::temperature*1.3806503e-23); // get anisotropy constant for material 0
         double energyH = sim::H_applied*mp::material[0].mu_s_SI/(sim::temperature*1.3806503e-23); // get anisotropy constant for material 0
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
         double energyK = anisotropy::get_anisotropy_constant(0)/(sim::temperature*1.3806503e-23); // get anisotropy constant for material 0
         double energyH = sim::H_applied*mp::material[0].mu_s_SI/(sim::temperature*1.3806503e-23); // get anisotropy constant for material 0
         double PK = sin(double (b)*M_PI/180.0)*exp(energyK*sin(double (b)*M_PI/180.0)*sin(double (b)*M_PI/180.0));
         double PH = sin(double (b)*M_PI/180.0)*exp(energyH*cos(double (b)*M_PI/180.0));
         if((bin[b])>maxP) maxP=bin[b];
         if(PK>maxPK) maxPK=PK;
         if(PH>maxPH) maxPH=PH;
      }

      std::ofstream ofile("boltzmann-distribution.txt");

      // Output data
      ofile << "# Anisotropy:          " << anisotropy::get_anisotropy_constant(0) << "\t" << anisotropy::get_anisotropy_constant(0)/(sim::temperature*1.3806503e-23) << std::endl;
      ofile << "# Field:               " << sim::H_applied << "\t" << sim::H_applied*mp::material[0].mu_s_SI/(sim::temperature*1.3806503e-23) << std::endl;
      ofile << "# Moment, Temperature: " << mp::material[0].mu_s_SI/9.274e-24 << "\t" << sim::temperature << std::endl;
      for(int b=0;b<181;b++){
         double energyK = anisotropy::get_anisotropy_constant(0)/(sim::temperature*1.3806503e-23); // get anisotropy constant for material 0
         double energyH = sim::H_applied*mp::material[0].mu_s_SI/(sim::temperature*1.3806503e-23); // get anisotropy constant for material 0
         double PK = sin(double (b)*M_PI/180.0)*exp(energyK*sin(double (b)*M_PI/180.0)*sin(double (b)*M_PI/180.0));
         double PH = sin(double (b)*M_PI/180.0)*exp(energyH*cos(double (b)*M_PI/180.0));
         ofile << b << "\t" << bin[b]/maxP << "\t" << (bin[b]+bin[180-b])/(2.0*maxP) << "\t" << PK/maxPK << "\t" << PH/maxPH << std::endl;
      }

      // close Boltzmann file
      ofile.close();

   }

}//end of namespace program
