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
#ifndef PROGRAM_H_
#define PROGRAM_H_

//==========================================================
// Namespace program
//==========================================================
namespace program
{
	// program functions
	extern int bmark();
	extern void time_series();
	extern int hysteresis();
	extern int static_hysteresis();
   extern void partial_hysteresis_loop();
	extern int curie_temperature();
	extern void field_cool();
	extern void temperature_pulse();
	extern void hamr();
	extern void cmc_anisotropy();
	extern void hybrid_cmc();
   extern void reverse_hybrid_cmc();
   extern void lagrange_multiplier();
   extern void localised_temperature_pulse();
   extern void effective_damping();
   extern void fmr();

	// Sundry programs and diagnostics not under general release
	extern int LLB_Boltzmann();
	extern int timestep_scaling();
	extern void boltzmann_dist();
	
}

#endif /*PROGRAM_H_*/
