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

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "montecarlo.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

namespace program{

/// Function to update lagrange lambda
void update_lagrange_lambda(){

   // Save initial value of lambda
   const double lamda_old_x = program::lagrange_lambda_x;
   const double lamda_old_y = program::lagrange_lambda_y;
   const double lamda_old_z = program::lagrange_lambda_z;

   // Calculate magnetisation
   stats::update();

   //std::cout << "xx " << program::lagrange_m << std::endl;

   const std::vector<double> mm = stats::system_magnetization.get_magnetization();

   // Calculate unit vector of magnetisation
   const double mx = mm[0];
   const double my = mm[1];
   const double mz = mm[2];
   program::lagrange_m = mm[3];

   // Constraint vector
   const double nu_x=cos(montecarlo::constraint_theta*M_PI/180.0)*sin(montecarlo::constraint_phi*M_PI/180.0);
   const double nu_y=sin(montecarlo::constraint_theta*M_PI/180.0)*sin(montecarlo::constraint_phi*M_PI/180.0);
   const double nu_z=cos(montecarlo::constraint_phi*M_PI/180.0);

   double lx = -program::lagrange_N*(mx-nu_x);
   double ly = -program::lagrange_N*(my-nu_y);
   double lz = -program::lagrange_N*(mz-nu_z);

   program::lagrange_lambda_x = lamda_old_x + lx*mp::dt;
   program::lagrange_lambda_y = lamda_old_y + ly*mp::dt;
   program::lagrange_lambda_z = lamda_old_z + lz*mp::dt;

}

}

namespace program{

//-----------------------------------------------------------
///  Function to perform constrained energy minimisation
///  using the LaGrange multiplier method
//
///  (c) R F L Evans 2013
//
///  D.A.Garanin Phys. Rev. Lett. 90, 065504 (2003)
//
//-----------------------------------------------------------
void lagrange_multiplier(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "program::lagrange_multiplier has been called" << std::endl;}

   // Enable LMM fields
   program::lagrange_multiplier_enabled=true;

   // Set prefactor in LaGrange multiplier (Tesla)
   program::lagrange_N=10.0;

   // Initialise LaGrange m
   program::update_lagrange_lambda();

   // Set temperature to zero
   sim::temperature=0.0;
   sim::thermal_field_enabled = false;

   // set minimum rotational angle
   montecarlo::constraint_theta=montecarlo::constraint_theta_min;

   // perform rotational angle sweep
   while(montecarlo::constraint_theta<=montecarlo::constraint_theta_max){

      // set minimum azimuthal angle
      montecarlo::constraint_phi=montecarlo::constraint_phi_min;

      // perform azimuthal angle sweep
      while(montecarlo::constraint_phi<=montecarlo::constraint_phi_max){

         zlog << zTs() << "Constraint of theta =  " << montecarlo::constraint_theta << ", phi = " << montecarlo::constraint_phi << std::endl;

         // Reset start time
         int start_time=sim::time;

         // Simulate system
         while(sim::time<sim::loop_time+start_time){

            // Integrate system
            sim::integrate(sim::partial_time);

            // Check for torque criteria
            double torque=stats::max_torque();

            //std::cout << sim::time << "\t" << torque << std::endl;

            if((torque<1.0e-6) && (sim::time-start_time>100)){
               break;
            }

         }

         // Calculate magnetisation statistics
         stats::update();

         // Output data
         vout::data();

         // Increment azimuthal angle
         montecarlo::constraint_phi+=montecarlo::constraint_phi_delta;
         montecarlo::constraint_phi_changed=true;

      } // End of azimuthal angle sweep
      if(vout::gnuplot_array_format) zmag << std::endl;

      // Increment rotational angle
      montecarlo::constraint_theta+=montecarlo::constraint_theta_delta;
      montecarlo::constraint_theta_changed=true;

   } // End of rotational angle sweep

   return;
}

}//end of namespace program
