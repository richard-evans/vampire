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
/// @brief Contains functions to calculate energy for a spin/system
///
/// @details None
///
/// @section notes Implementation Notes
/// This is a list of other notes, not related to functionality but rather to implementation.
/// Also include references, formulae and other notes here.
///
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
/// @internal
///   Created:    07/02/2011
///   Revision:     ---
///=====================================================================================
///

// Standard Libraries
#include <algorithm>
#include <cmath>
#include <iostream>

// Vampire Header files
#include "anisotropy.hpp"
#include "atoms.hpp"
#include "exchange.hpp"
#include "integrate.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "dipole.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "spintransport.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "constants.hpp"

// sim module header
#include "internal.hpp"

namespace sim{



/// @brief Calculates the applied field energy for a single spin.
///
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @param[in] atom atom number
/// @param[in] imaterial material of local atom
/// @param[in] Sx x-spin of local atom
/// @param[in] Sy y-spin of local atom
/// @param[in] Sz z-spin of local atom
/// @return applied field energy
///
/// @internal
///   Created:    07/02/2011
///   Revision:     ---
///=====================================================================================
///
double spin_applied_field_energy(const double Sx, const double Sy, const double Sz){

   return -sim::H_applied*(sim::H_vec[0]*Sx + sim::H_vec[1]*Sy + sim::H_vec[2]*Sz);

}

double spin_local_applied_field_energy(const int mat, const double sx, const double sy, const double sz){

   const double B = mp::material[mat].applied_field_strength;
   const double hx = B * mp::material[mat].applied_field_unit_vector[0];
   const double hy = B * mp::material[mat].applied_field_unit_vector[1];
   const double hz = B * mp::material[mat].applied_field_unit_vector[2];

   return -(hx*sx + hy*sy + hz*sz);

}

/// @brief Calculates the magnetostatic energy for a single spin.
///
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    08/06/2012
///
/// @param[in] atom atom number
/// @param[in] imaterial material of local atom
/// @param[in] Sx x-spin of local atom
/// @param[in] Sy y-spin of local atom
/// @param[in] Sz z-spin of local atom
/// @return magnetostatic energy
///
/// @internal
///   Created:    08/06/2012
///   Revision:     ---
///=====================================================================================
///
double spin_magnetostatic_energy(const int atom, const double Sx, const double Sy, const double Sz){
   //return -1.0*(dipole::atom_dipolar_field_array_x[atom]*Sx+dipole::atom_dipolar_field_array_y[atom]*Sy+dipole::atom_dipolar_field_array_z[atom]*Sz);
   return -1.0*(dipole::atom_mu0demag_field_array_x[atom]*Sx+dipole::atom_mu0demag_field_array_y[atom]*Sy+dipole::atom_mu0demag_field_array_z[atom]*Sz);
}

// Calculates LSF energy
double spin_longitudinal_energy(const int atom){

   // Standard Landau Hamiltonian
   const int imaterial = atoms::type_array[atom];
   const double Sx = atoms::x_spin_array[atom];
   const double Sy = atoms::y_spin_array[atom];
   const double Sz = atoms::z_spin_array[atom];
   double spinlength = sqrt(Sx * Sx + Sy * Sy + Sz * Sz);
   double mod_S2_i = spinlength * spinlength;
   double mod_S4_i = mod_S2_i * mod_S2_i;
   double mod_S6_i = mod_S4_i * mod_S2_i;

   double A, B, C;
   integrate::get_lsf_coefficients(imaterial, A, B, C);

   double landau_energy = (1 / mp::material[imaterial].mu_s_SI) * ((A * mod_S2_i) + (B * mod_S4_i) + (C * mod_S6_i));

   return landau_energy; // Tesla
}

/// @brief Calculates the total energy for a single spin.
///
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @param[in] atom atom number
/// @return total spin energy
///
/// @internal
///   Created:    07/02/2011
///   Revision:     ---
///=====================================================================================
///
double calculate_spin_energy(const int atom){

   // check calling of routine if error checking is activated
   if(err::check==true) std::cout << "calculate_spin_energy has been called" << std::endl;

   // Local spin value
   const double Sx=atoms::x_spin_array[atom];
   const double Sy=atoms::y_spin_array[atom];
   const double Sz=atoms::z_spin_array[atom];

   // Determine neighbour material
   const int imaterial=atoms::type_array[atom];

   // Initialise energy to zero
   double energy=0.0;

   // Calculate total spin energy
   energy += exchange::single_spin_energy(atom, Sx, Sy, Sz);
   energy += exchange::single_spin_biquadratic_energy(atom, Sx, Sy, Sz);
   energy += exchange::single_spin_four_spin_energy(atom, Sx, Sy, Sz);

   // calculate anisotropy energy for atom
   energy += anisotropy::single_spin_energy(atom, imaterial, Sx, Sy, Sz, sim::temperature);

   // Landau energy
   if(integrate::integrator==integrate::integrator_t::lsf_mc){
      energy += spin_longitudinal_energy(atom);
   }

   energy+=spin_applied_field_energy(Sx, Sy, Sz);
   energy+=spin_magnetostatic_energy(atom, Sx, Sy, Sz);

   // local applied fields
   if(sim::local_applied_field) energy += spin_local_applied_field_energy(imaterial, Sx, Sy, Sz);

   // vcma energy
   const double vcma = program::fractional_electric_field_strength * spin_transport::get_voltage() * sim::internal::vcmak[imaterial];
   energy -= vcma * Sz * Sz;

   return energy; // Tesla
}

} // end of namespace sim
