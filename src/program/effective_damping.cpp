//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2014. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>
#include <cmath>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// program module headers
#include "internal.hpp"

namespace program{

//------------------------------------------------------------------------------
// Function to rotate all spins by angle (degrees) around a unit vector axis
// using the Rodrigues rotation formula
//
//    S' = S cos(phi) + (u x S) sin(phi) + u (u . S) (1 - cos(phi))
//
//------------------------------------------------------------------------------
void rotate_spins_around_axis(const std::vector<double>& axis, const double angle){

   const double phi = angle * M_PI / 180.0;
   const double cos_phi = cos(phi);
   const double sin_phi = sin(phi);
   const double one_minus_cos_phi = 1.0 - cos_phi;

   // unit vector axis components (normalised on input)
   const double u = axis[0];
   const double v = axis[1];
   const double w = axis[2];

   // loop over all spins and rotate by phi around the axis
   for(int atom = 0; atom < atoms::num_atoms; atom++){

      // Load spin coordinates
      const double x = atoms::x_spin_array[atom];
      const double y = atoms::y_spin_array[atom];
      const double z = atoms::z_spin_array[atom];

      const double dot = x*u + y*v + z*w;

      // Set new spin positions
      atoms::x_spin_array[atom] = u * dot * one_minus_cos_phi + x * cos_phi + (-w * y + v * z) * sin_phi;
      atoms::y_spin_array[atom] = v * dot * one_minus_cos_phi + y * cos_phi + ( w * x - u * z) * sin_phi;
      atoms::z_spin_array[atom] = w * dot * one_minus_cos_phi + z * cos_phi + (-v * x + u * y) * sin_phi;

   }

   return;

}

//-----------------------------------------------------------------------------
//
//   Program to calculate effective gilbert damping from the atomistic model.
//
//   System is first equilibrated at sim:temperature along the z-axis. All
//   spins are then rotated around sim:effective-damping-rotation-axis (default x)
//   by sim:effective-damping-rotation-angle (default 30 degrees) and the magnetisation
//   then relaxes towards the z-direction under the application of a field
//   or magnetocrystalline anisotropy. The magnetization trace can then be
//   fitted to a macroscopic equation to extract the effective damping, alpha.
//
//   Ref. M O A Ellis et al, Physical Review B 86, 174418 (2012)
//
//-----------------------------------------------------------------------------
void effective_damping(){

   // check calling of routine if error checking is activated
   if(err::check==true) std::cout << "program::effective_damping has been called" << std::endl;

   // save system temperature
   double temp=sim::temperature;

   // Set equilibration temperature and field
   sim::temperature=sim::Teq;

   // Equilibrate system
   while(sim::time<sim::equilibration_time){

      sim::integrate(sim::partial_time);

      // Calculate magnetisation statistics
      stats::update();

      // Output data
      vout::data();
   }

   // Rotate all spins around the user defined axis by the user defined angle (default 30 degrees about the x-axis)
   rotate_spins_around_axis(program::internal::effective_damping_rotation_axis, program::internal::effective_damping_rotation_angle);

   // reset system temperature
   sim::temperature=temp;

   // Perform Time Series with relaxation
   while(sim::time<sim::equilibration_time+sim::total_time){

      // Integrate system
      sim::integrate(sim::partial_time);

      // Calculate magnetisation statistics
      stats::update();

      // Output data
      vout::data();

   }

}

}//end of namespace program
