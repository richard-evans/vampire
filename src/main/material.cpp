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
#include "material.hpp"

namespace mp{
	// Constructor
materials_t::materials_t ():
	name("material#n"),
	element("Fe"),
	alpha(1.0),
	mu_s_SI(1.72*9.27400915e-24),
	magnetisation(0.0),
	gamma_rel(1.0),
	random_spins(false),
	geometry(0),
	core_shell_size(1.0),
	interface_roughness(0.0),
	density(1.0),
	cutoff(0.8),
	alloy_master(false),
	alloy_class(-1),
	continuous(false),
	moment_flag(true),
	one_oneplusalpha_sq(0.5),
	alpha_oneplusalpha_sq(0.5),
	H_th_sigma(0.0),
	constrained(false),
	temperature(0.0),
	maximum_temperature(0.0),
	minimum_temperature(0.0),
	couple_to_phonon_temperature(false),
	applied_field_strength(0.0),
	applied_field_unit_vector(3,0.0),
	fmr_field_strength(0.0),
	fmr_field_frequency(0.0),
	fmr_field_unit_vector(3,0.0),
   fill(false),
   temperature_rescaling_alpha(1.0),
	temperature_rescaling_Tc(0.0),
	non_magnetic(0)
	{

	//std::cout << "constructor " << anis_flag << "\t" << ianis_flag << std::endl;
	// derived parameters
	for(int i=0;i<100;i++){
		geometry_coords[i][0]=0.0;
		geometry_coords[i][1]=0.0;
	}

	SAF.resize(mp::max_materials);
	enable_SAF = false;
	// array variables

 	EF_MM.resize(mp::max_materials);
 	override_atomsitic.resize(mp::max_materials);

	// array variables
   for(int i=0; i<mp::max_materials; i++){
		intermixing[i]=0.0;
		alloy[i]=0.0;
	}
	initial_spin[0]=0.0;
	initial_spin[1]=0.0;
	initial_spin[2]=1.0;

	// Applied field direction default initialisation
	applied_field_unit_vector.at(0)=0.0;
	applied_field_unit_vector.at(1)=0.0;
	applied_field_unit_vector.at(2)=1.0;

	pinning_field_unit_vector.resize(3,0.0);
	
	// FMR field direction default initialisation
	fmr_field_unit_vector.at(0)=0.0;
	fmr_field_unit_vector.at(1)=0.0;
	fmr_field_unit_vector.at(2)=1.0;
}

int materials_t::print(){

	std::cout << "----------------------------------------------------------------" << std::endl;
	std::cout << " Material " << name << std::endl;
	std::cout << "----------------------------------------------------------------" << std::endl;
	std::cout << "alpha          = " << alpha << std::endl;
	std::cout << "mu_s_SI        = " << mu_s_SI << std::endl;
	std::cout << "gamma_rel      = " << gamma_rel << std::endl;

	return 0;

}

}
