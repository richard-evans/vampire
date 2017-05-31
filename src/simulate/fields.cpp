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
//====================================================================================================
//
//       				                    	Fields
//
//  			 		Subroutines to calculate fields for the hamiltonian
//
//									Version 1.0 R Evans 20/10/2008
//
//====================================================================================================
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
//#include "demag.hpp"
#include "dipole.hpp"
#include "ltmp.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "spintorque.hpp"
#include "stats.hpp"
#include "vmpi.hpp"

// sim module header
#include "internal.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

//========================
//function prototypes
//========================

int calculate_exchange_fields(const int,const int);
int calculate_anisotropy_fields(const int,const int);
void calculate_second_order_uniaxial_anisotropy_fields(const int,const int);
void calculate_sixth_order_uniaxial_anisotropy_fields(const int,const int);
void calculate_spherical_harmonic_fields(const int,const int);
void calculate_random_spherical_harmonic_fields(const int,const int);
void calculate_lattice_anisotropy_fields(const int, const int);
int calculate_cubic_anisotropy_fields(const int,const int);
int calculate_applied_fields(const int,const int);
int calculate_thermal_fields(const int,const int);
int calculate_dipolar_fields(const int,const int);
void calculate_hamr_fields(const int,const int);
void calculate_fmr_fields(const int,const int);
void calculate_surface_anisotropy_fields(const int,const int);
void calculate_lagrange_fields(const int,const int);
void calculate_full_spin_fields(const int start_index,const int end_index);


namespace sim{

void calculate_spin_fields(const int start_index,const int end_index){
	///======================================================
	/// 		Subroutine to calculate spin dependent fields
	///
	///			Version 1.0 R Evans 20/10/2008
	///======================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_spin_fields has been called" << std::endl;}

	// Initialise Total Spin Fields to zero
	fill (atoms::x_total_spin_field_array.begin()+start_index,atoms::x_total_spin_field_array.begin()+end_index,0.0);
	fill (atoms::y_total_spin_field_array.begin()+start_index,atoms::y_total_spin_field_array.begin()+end_index,0.0);
	fill (atoms::z_total_spin_field_array.begin()+start_index,atoms::z_total_spin_field_array.begin()+end_index,0.0);

	// Exchange Fields
	if(sim::hamiltonian_simulation_flags[0]==1) calculate_exchange_fields(start_index,end_index);

	// Anisotropy Fields
	if(sim::UniaxialScalarAnisotropy || sim::TensorAnisotropy) calculate_anisotropy_fields(start_index,end_index);
   if(sim::second_order_uniaxial_anisotropy) calculate_second_order_uniaxial_anisotropy_fields(start_index,end_index);
   if(sim::sixth_order_uniaxial_anisotropy) calculate_sixth_order_uniaxial_anisotropy_fields(start_index,end_index);
   if(sim::spherical_harmonics && sim::random_anisotropy==false) calculate_spherical_harmonic_fields(start_index,end_index);
   if(sim::random_anisotropy && sim::spherical_harmonics) calculate_random_spherical_harmonic_fields(start_index,end_index);
   if(sim::lattice_anisotropy_flag) calculate_lattice_anisotropy_fields(start_index,end_index);
   if(sim::CubicScalarAnisotropy) calculate_cubic_anisotropy_fields(start_index,end_index);
	//if(sim::hamiltonian_simulation_flags[1]==3) calculate_local_anis_fields();
	if(sim::surface_anisotropy==true) calculate_surface_anisotropy_fields(start_index,end_index);
	// Spin Dependent Extra Fields
	if(sim::lagrange_multiplier==true) calculate_lagrange_fields(start_index,end_index);

	calculate_full_spin_fields(start_index,end_index);

}

void calculate_external_fields(const int start_index,const int end_index){
	///======================================================
	/// 		Subroutine to calculate external fields
	///
	///			Version 1.0 R Evans 20/10/2008
	///======================================================
	//const int num_atoms = atoms::num_atoms;

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_external_fields has been called" << std::endl;}

	// Initialise Total External Fields to zero
	fill (atoms::x_total_external_field_array.begin()+start_index,atoms::x_total_external_field_array.begin()+end_index,0.0);
	fill (atoms::y_total_external_field_array.begin()+start_index,atoms::y_total_external_field_array.begin()+end_index,0.0);
	fill (atoms::z_total_external_field_array.begin()+start_index,atoms::z_total_external_field_array.begin()+end_index,0.0);

	if(sim::program==7) calculate_hamr_fields(start_index,end_index);
   else if(sim::program==13){

      // Local thermal Fields
      ltmp::get_localised_thermal_fields(atoms::x_total_external_field_array,atoms::y_total_external_field_array,
                                         atoms::z_total_external_field_array, start_index, end_index);

      // Applied Fields
      if(sim::hamiltonian_simulation_flags[2]==1) calculate_applied_fields(start_index,end_index);

   }
	else{

		// Thermal Fields
		if(sim::hamiltonian_simulation_flags[3]==1) calculate_thermal_fields(start_index,end_index);

		// Applied Fields
		if(sim::hamiltonian_simulation_flags[2]==1) calculate_applied_fields(start_index,end_index);

	}

   // Get updated spin torque fields
   st::get_spin_torque_fields(atoms::x_total_external_field_array, atoms::y_total_external_field_array, atoms::z_total_external_field_array, start_index, end_index);

	// FMR Fields only for fmr program
	if(sim::enable_fmr) calculate_fmr_fields(start_index,end_index);

	// Dipolar Fields
	calculate_dipolar_fields(start_index,end_index);

}
}

int calculate_exchange_fields(const int start_index,const int end_index){
	///======================================================
	/// 		Subroutine to calculate exchange fields
	///
	///			Version 2.0 Richard Evans 08/09/2011
	///======================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_exchange_fields has been called" << std::endl;}

	// Use appropriate function for exchange calculation
	switch(atoms::exchange_type){
		case 0: // isotropic
			for(int atom=start_index;atom<end_index;atom++){
				double Hx=0.0;
				double Hy=0.0;
				double Hz=0.0;
				const int start=atoms::neighbour_list_start_index[atom];
				const int end=atoms::neighbour_list_end_index[atom]+1;
				for(int nn=start;nn<end;nn++){
					const int natom = atoms::neighbour_list_array[nn];
					const double Jij=atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij;
					Hx -= Jij*atoms::x_spin_array[natom];
					Hy -= Jij*atoms::y_spin_array[natom];
					Hz -= Jij*atoms::z_spin_array[natom];
				}
				atoms::x_total_spin_field_array[atom] += Hx;
				atoms::y_total_spin_field_array[atom] += Hy;
				atoms::z_total_spin_field_array[atom] += Hz;
			}
			break;
		case 1: // vector
			for(int atom=start_index;atom<end_index;atom++){
				double Hx=0.0;
				double Hy=0.0;
				double Hz=0.0;
				const int start=atoms::neighbour_list_start_index[atom];
				const int end=atoms::neighbour_list_end_index[atom]+1;
				for(int nn=start;nn<end;nn++){
					const int natom = atoms::neighbour_list_array[nn];
					const int iid = atoms::neighbour_interaction_type_array[nn]; // interaction id
					const double Jij[3]={atoms::v_exchange_list[iid].Jij[0],
												atoms::v_exchange_list[iid].Jij[1],
												atoms::v_exchange_list[iid].Jij[2]};

					Hx -= Jij[0]*atoms::x_spin_array[natom];
					Hy -= Jij[1]*atoms::y_spin_array[natom];
					Hz -= Jij[2]*atoms::z_spin_array[natom];
				}
				atoms::x_total_spin_field_array[atom] += Hx;
				atoms::y_total_spin_field_array[atom] += Hy;
				atoms::z_total_spin_field_array[atom] += Hz;
			}
			break;
		case 2: // tensor
			for(int atom=start_index;atom<end_index;atom++){
				double Hx=0.0;
				double Hy=0.0;
				double Hz=0.0;
				const int start=atoms::neighbour_list_start_index[atom];
				const int end=atoms::neighbour_list_end_index[atom]+1;
				for(int nn=start;nn<end;nn++){
					const int natom = atoms::neighbour_list_array[nn];
					const int iid = atoms::neighbour_interaction_type_array[nn]; // interaction id
					const double Jij[3][3]={{atoms::t_exchange_list[iid].Jij[0][0],
													 atoms::t_exchange_list[iid].Jij[0][1],
													 atoms::t_exchange_list[iid].Jij[0][2]},

													{atoms::t_exchange_list[iid].Jij[1][0],
													 atoms::t_exchange_list[iid].Jij[1][1],
													 atoms::t_exchange_list[iid].Jij[1][2]},

													{atoms::t_exchange_list[iid].Jij[2][0],
													 atoms::t_exchange_list[iid].Jij[2][1],
													 atoms::t_exchange_list[iid].Jij[2][2]}};

					const double S[3]={atoms::x_spin_array[natom],atoms::y_spin_array[natom],atoms::z_spin_array[natom]};

					Hx -= (Jij[0][0]*S[0] + Jij[0][1]*S[1] +Jij[0][2]*S[2]);
					Hy -= (Jij[1][0]*S[0] + Jij[1][1]*S[1] +Jij[1][2]*S[2]);
					Hz -= (Jij[2][0]*S[0] + Jij[2][1]*S[1] +Jij[2][2]*S[2]);
				}
				atoms::x_total_spin_field_array[atom] += Hx;
				atoms::y_total_spin_field_array[atom] += Hy;
				atoms::z_total_spin_field_array[atom] += Hz;
			}
			break;
		}

		return EXIT_SUCCESS;
	}

int calculate_anisotropy_fields(const int start_index,const int end_index){
	///======================================================
	/// 	Subroutine to calculate uniaxial anisotropy fields
	///
	///			Version 1.0 R Evans 20/10/2008
	///======================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_anisotropy_fields has been called" << std::endl;}

		// Use appropriate function for anisotropy calculation
	switch(sim::AnisotropyType){
		case 0: // scalar
			for(int atom=start_index;atom<end_index;atom++){
				const int imaterial=atoms::type_array[atom];
				atoms::z_total_spin_field_array[atom] -= 2.0*mp::MaterialScalarAnisotropyArray[imaterial].K*atoms::z_spin_array[atom];
			}
			break;
		case 1: // tensor
			for(int atom=start_index;atom<end_index;atom++){
				const int imaterial=atoms::type_array[atom];

				const double K[3][3]={{2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[0][0],
											  2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[0][1],
											  2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[0][2]},

											 {2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[1][0],
											  2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[1][1],
											  2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[1][2]},

										    {2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[2][0],
											  2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[2][1],
											  2.0*mp::MaterialTensorAnisotropyArray[imaterial].K[2][2]}};

				const double S[3]={atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};

				atoms::x_total_spin_field_array[atom] -= (K[0][0]*S[0] + K[0][1]*S[1] +K[0][2]*S[2]);
				atoms::y_total_spin_field_array[atom] -= (K[1][0]*S[0] + K[1][1]*S[1] +K[1][2]*S[2]);
				atoms::z_total_spin_field_array[atom] -= (K[2][0]*S[0] + K[2][1]*S[1] +K[2][2]*S[2]);
			}
			break;
   }
	return EXIT_SUCCESS;
}

///------------------------------------------------------
///  Function to calculate second order uniaxial
///  anisotropy fields
///
///  (c) R F L Evans 2013
///
///  E = k4*(S . e)^4
///  Hx = -4*k4*(S . e)^3 e_x
///  Hy = -4*k4*(S . e)^3 e_y
///  Hz = -4*k4*(S . e)^3 e_z
///
///------------------------------------------------------
void calculate_second_order_uniaxial_anisotropy_fields(const int start_index,const int end_index){
   for(int atom=start_index;atom<end_index;atom++){
      const int imaterial=atoms::type_array[atom];
      const double ex = mp::material.at(imaterial).UniaxialAnisotropyUnitVector.at(0);
      const double ey = mp::material.at(imaterial).UniaxialAnisotropyUnitVector.at(1);
      const double ez = mp::material.at(imaterial).UniaxialAnisotropyUnitVector.at(2);
      const double Sx = atoms::x_spin_array[atom];
      const double Sy = atoms::y_spin_array[atom];
      const double Sz = atoms::z_spin_array[atom];
      const double Ku2 = 4.0*mp::material_second_order_anisotropy_constant_array[imaterial];
      const double Sdote = (Sx*ex + Sy*ey + Sz*ez);
      const double Sdote3 = Sdote*Sdote*Sdote;

      atoms::x_total_spin_field_array[atom] -= Ku2*ex*Sdote3;
      atoms::y_total_spin_field_array[atom] -= Ku2*ey*Sdote3;
      atoms::z_total_spin_field_array[atom] -= Ku2*ez*Sdote3;
   }
   return;
}

///------------------------------------------------------
///  Function to calculate sixth order uniaxial
///  anisotropy fields
///
///  (c) R F L Evans 2013
///
///  E = k6*(S . e)^6
///  Hx = -6*k6*(S . e)^5 e_x
///  Hy = -6*k6*(S . e)^5 e_y
///  Hz = -6*k6*(S . e)^5 e_z
///
///------------------------------------------------------
void calculate_sixth_order_uniaxial_anisotropy_fields(const int start_index,const int end_index){
   for(int atom=start_index;atom<end_index;atom++){
      const int imaterial=atoms::type_array[atom];
      const double ex = mp::material.at(imaterial).UniaxialAnisotropyUnitVector.at(0);
      const double ey = mp::material.at(imaterial).UniaxialAnisotropyUnitVector.at(1);
      const double ez = mp::material.at(imaterial).UniaxialAnisotropyUnitVector.at(2);
      const double Sx = atoms::x_spin_array[atom];
      const double Sy = atoms::y_spin_array[atom];
      const double Sz = atoms::z_spin_array[atom];
      const double Ku3 = 6.0*mp::material_sixth_order_anisotropy_constant_array[imaterial];
      const double Sdote = (Sx*ex + Sy*ey + Sz*ez);
      const double Sdote5 = Sdote*Sdote*Sdote*Sdote*Sdote;

      atoms::x_total_spin_field_array[atom] -= Ku3*ex*Sdote5;
      atoms::y_total_spin_field_array[atom] -= Ku3*ey*Sdote5;
      atoms::z_total_spin_field_array[atom] -= Ku3*ez*Sdote5;
   }
   return;
}

///--------------------------------------------------------------------------------------------------------------
///  Function to calculate spherical harmonic anisotropy fields
///
///  (c) R F L Evans 2015
///
///  Higher order anisotropies generally need to be described using spherical harmonics. The usual form (a
///  series in S leads to cross pollution of terms, giving strange temperature dependencies.
///
///  The harmonics are described with Legendre polynomials with even order, which for 2nd, 4th and 6th are:
///  ( http://en.wikipedia.org/wiki/Legendre_polynomials )
///
///  k_2(sz) = (1/2) *(3sz^2 - 1)
///  k_4(sz) = (1/8) *(35sz^4 - 30sz^2 + 3)
///  k_6(sz) = (1/16)*(231sz^6 - 315*sz^4 + 105sz^2 - 5)
///
///  The harmonics feature an arbritrary 2/3 factor compared with the usual form, and so in VAMPIRE these are
///  renormalised to maintain consistency for the 2nd order terms.
///
///  The field induced by the harmonics is given by the first derivative w.r.t. sz. This can be projected onto
///  any arbritrary direction ex,ey,ez allowing higher order anisotropy terms along any direction. This
///  direction is shared with the other uniaxial anisotropy coefficients since they should not be used
///  simultaneously.
///
///--------------------------------------------------------------------------------------------------------------
void calculate_spherical_harmonic_fields(const int start_index,const int end_index){

   // rescaling prefactor
   const double scale = 2.0/3.0; // Factor to rescale anisotropies to usual scale

   // constant factors
   const double oneo8 = 1.0/8.0;
   const double oneo16 = 1.0/16.0;

   // loop over all atoms
   for(int atom=start_index; atom<end_index; atom++){

      // Determine atom type
      const int imaterial=atoms::type_array[atom];

      // determine harmonic constants for material
      const double k2 = mp::material_spherical_harmonic_constants_array[3*imaterial + 0];
      const double k4 = mp::material_spherical_harmonic_constants_array[3*imaterial + 1];
      const double k6 = mp::material_spherical_harmonic_constants_array[3*imaterial + 2];

      // determine anisotropy direction and dot product
      const double ex = mp::material[imaterial].UniaxialAnisotropyUnitVector[0];
      const double ey = mp::material[imaterial].UniaxialAnisotropyUnitVector[1];
      const double ez = mp::material[imaterial].UniaxialAnisotropyUnitVector[2];
      const double sx = atoms::x_spin_array[atom];
      const double sy = atoms::y_spin_array[atom];
      const double sz = atoms::z_spin_array[atom];

      const double sdote = (sx*ex + sy*ey + sz*ez);
      const double sdote3 = sdote*sdote*sdote;
      const double sdote5 = sdote3*sdote*sdote;

      // calculate field (double negative from scale factor and negative derivative)
      atoms::x_total_spin_field_array[atom] += scale*ex*(k2*3.0*sdote + k4*oneo8*(140.0*sdote3 - 60.0*sdote) + k6*oneo16*(1386.0*sdote5 - 1260.0*sdote3 + 210.0*sdote));
      atoms::y_total_spin_field_array[atom] += scale*ey*(k2*3.0*sdote + k4*oneo8*(140.0*sdote3 - 60.0*sdote) + k6*oneo16*(1386.0*sdote5 - 1260.0*sdote3 + 210.0*sdote));
      atoms::z_total_spin_field_array[atom] += scale*ez*(k2*3.0*sdote + k4*oneo8*(140.0*sdote3 - 60.0*sdote) + k6*oneo16*(1386.0*sdote5 - 1260.0*sdote3 + 210.0*sdote));

   }

   return;

}

///--------------------------------------------------------------------------------------------------------------
///  Function to calculate random spherical harmonic anisotropy fields
///
///  (c) R F L Evans 2015
///
///  In this function uniaxial anisotropy is calculated using spherical harmonics,
///  except each atom is allowed a locally defined anisotropy axis. This comes with
///  a performance cost, and so this version is only caled if needed (defined by the
///  sim::random_anisotropy flag).
///
///--------------------------------------------------------------------------------------------------------------
void calculate_random_spherical_harmonic_fields(const int start_index,const int end_index){

  // rescaling prefactor
  const double scale = 2.0/3.0; // Factor to rescale anisotropies to usual scale

  // constant factors
  const double oneo8 = 1.0/8.0;
  const double oneo16 = 1.0/16.0;

  // loop over all atoms
  for(int atom=start_index; atom<end_index; atom++){

    // Determine atom type
    const int imaterial=atoms::type_array[atom];

    // determine harmonic constants for material
    const double k2 = mp::material_spherical_harmonic_constants_array[3*imaterial + 0];
    const double k4 = mp::material_spherical_harmonic_constants_array[3*imaterial + 1];
    const double k6 = mp::material_spherical_harmonic_constants_array[3*imaterial + 2];

    // determine anisotropy direction and dot product
    const double ex = atoms::uniaxial_anisotropy_vector_x[atom];
    const double ey = atoms::uniaxial_anisotropy_vector_y[atom];
    const double ez = atoms::uniaxial_anisotropy_vector_z[atom];
    const double sx = atoms::x_spin_array[atom];
    const double sy = atoms::y_spin_array[atom];
    const double sz = atoms::z_spin_array[atom];

    const double sdote = (sx*ex + sy*ey + sz*ez);
    const double sdote3 = sdote*sdote*sdote;
    const double sdote5 = sdote3*sdote*sdote;

    // calculate field (double negative from scale factor and negative derivative)
    atoms::x_total_spin_field_array[atom] += scale*ex*(k2*3.0*sdote + k4*oneo8*(140.0*sdote3 - 60.0*sdote) + k6*oneo16*(1386.0*sdote5 - 1260.0*sdote3 + 210.0*sdote));
    atoms::y_total_spin_field_array[atom] += scale*ey*(k2*3.0*sdote + k4*oneo8*(140.0*sdote3 - 60.0*sdote) + k6*oneo16*(1386.0*sdote5 - 1260.0*sdote3 + 210.0*sdote));
    atoms::z_total_spin_field_array[atom] += scale*ez*(k2*3.0*sdote + k4*oneo8*(140.0*sdote3 - 60.0*sdote) + k6*oneo16*(1386.0*sdote5 - 1260.0*sdote3 + 210.0*sdote));

  }

  return;

}

//------------------------------------------------------
///  Function to calculate lattice anisotropy fields
//
///  (c) R F L Evans 2013
//
//------------------------------------------------------
void calculate_lattice_anisotropy_fields(const int start_index,const int end_index){

   // Precalculate material lattice anisotropy constants
   std::vector<double> klatt_array(0);
   klatt_array.reserve(mp::num_materials);
   for(int imat=0; imat<mp::num_materials; imat++) klatt_array.push_back(2.0*mp::material[imat].Klatt*mp::material[imat].lattice_anisotropy.get_lattice_anisotropy_constant(sim::temperature));

   // Precalculate unit vectors
   std::vector<double> ex(0);
   std::vector<double> ey(0);
   std::vector<double> ez(0);

   ex.reserve(mp::num_materials);
   ey.reserve(mp::num_materials);
   ez.reserve(mp::num_materials);

   for(int imat=0; imat<mp::num_materials; imat++) ex.push_back(mp::material.at(imat).UniaxialAnisotropyUnitVector.at(0));
   for(int imat=0; imat<mp::num_materials; imat++) ey.push_back(mp::material.at(imat).UniaxialAnisotropyUnitVector.at(1));
   for(int imat=0; imat<mp::num_materials; imat++) ez.push_back(mp::material.at(imat).UniaxialAnisotropyUnitVector.at(2));

   // Now calculate fields
   for(int atom=start_index;atom<end_index;atom++){
      const int imaterial=atoms::type_array[atom];
      const double Sx = atoms::x_spin_array[atom];
      const double Sy = atoms::y_spin_array[atom];
      const double Sz = atoms::z_spin_array[atom];
      const double Sdote = (Sx*ex[imaterial] + Sy*ey[imaterial] + Sz*ez[imaterial]);

      atoms::x_total_spin_field_array[atom] -= klatt_array[imaterial]*ex[imaterial]*Sdote;
      atoms::y_total_spin_field_array[atom] -= klatt_array[imaterial]*ey[imaterial]*Sdote;
      atoms::z_total_spin_field_array[atom] -= klatt_array[imaterial]*ez[imaterial]*Sdote;

   }

   return;

}

int calculate_cubic_anisotropy_fields(const int start_index,const int end_index){
	///------------------------------------------------------
	/// 	Function to calculate cubic anisotropy fields
	///
	///			Version 1.0 R Evans 28/07/2012
	///
	///		E = -0.5 Kc (Sx^4 + Sy^4 + Sz^4)
	///		Hx = +2 Kc*(Sx^3)
	///		Hy = +2 Kc*(Sy^3)
	///		Hz = +2 Kc*(Sz^3)
	///
	///------------------------------------------------------
	//std::cout << "here" << std::endl;
	for(int atom=start_index;atom<end_index;atom++){
		const int imaterial=atoms::type_array[atom];
		const double Kc=2.0*mp::MaterialCubicAnisotropyArray[imaterial];

		const double Sx=atoms::x_spin_array[atom];
		atoms::x_total_spin_field_array[atom] -= Kc*Sx*Sx*Sx;

		const double Sy=atoms::y_spin_array[atom];
		atoms::y_total_spin_field_array[atom] -= Kc*Sy*Sy*Sy;

		const double Sz=atoms::z_spin_array[atom];
		atoms::z_total_spin_field_array[atom] -= Kc*Sz*Sz*Sz;

	}
	return EXIT_SUCCESS;
}

void calculate_surface_anisotropy_fields(const int start_index,const int end_index){
	///======================================================
	/// 		Subroutine to calculate surface anisotropy fields
	///
	///			Version 1.0 Richard Evans 13/09/2011
	///======================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_surface_anisotropy_fields has been called" << std::endl;}

	for(int atom=start_index;atom<end_index;atom++){
		// only calculate for surface atoms
		if(atoms::surface_array[atom]==true){
			const int imaterial=atoms::type_array[atom];
			const double Ks=0.5*2.0*mp::material[imaterial].Ks; // note factor two here from differentiation
			const double S[3]={atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};

			for(int nn=atoms::nearest_neighbour_list_si[atom];nn<atoms::nearest_neighbour_list_ei[atom];nn++){
				const double si_dot_eij=(S[0]*atoms::eijx[nn]+S[1]*atoms::eijy[nn]+S[2]*atoms::eijz[nn]);
				atoms::x_total_spin_field_array[atom]-=Ks*si_dot_eij*atoms::eijx[nn];
				atoms::y_total_spin_field_array[atom]-=Ks*si_dot_eij*atoms::eijy[nn];
				atoms::z_total_spin_field_array[atom]-=Ks*si_dot_eij*atoms::eijz[nn];
			}
		}
	}

	return;
}

int calculate_applied_fields(const int start_index,const int end_index){
	///==========================================================================
	///
	/// 	Function to calculate applied fields
	///
	///		Version 1.0 R Evans 20/10/2008
	///		Version 2.0 R F L Evans 18/11/2012
	///
	///==========================================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_applied_fields has been called" << std::endl;}

	// Declare constant temporaries for global field
	const double Hx=sim::H_vec[0]*sim::H_applied;
	const double Hy=sim::H_vec[1]*sim::H_applied;
	const double Hz=sim::H_vec[2]*sim::H_applied;

	// Declare array for local (material specific) applied field
	std::vector<double> Hlocal(0);

	// Check for local applied field
	if(sim::local_applied_field==true){
		Hlocal.reserve(3*mp::material.size());

		// Loop over all materials
		for(unsigned int mat=0;mat<mp::material.size();mat++){
			Hlocal.push_back(mp::material[mat].applied_field_strength*mp::material[mat].applied_field_unit_vector[0]);
			Hlocal.push_back(mp::material[mat].applied_field_strength*mp::material[mat].applied_field_unit_vector[1]);
			Hlocal.push_back(mp::material[mat].applied_field_strength*mp::material[mat].applied_field_unit_vector[2]);
		}

		// Add local field AND global field
		for(int atom=start_index;atom<end_index;atom++){
			const int imaterial=atoms::type_array[atom];
			atoms::x_total_external_field_array[atom] += Hx + Hlocal[3*imaterial + 0];
			atoms::y_total_external_field_array[atom] += Hy + Hlocal[3*imaterial + 1];
			atoms::z_total_external_field_array[atom] += Hz + Hlocal[3*imaterial + 2];
		}
	}
	else{
		// Calculate global field
		for(int atom=start_index;atom<end_index;atom++){
			atoms::x_total_external_field_array[atom] += Hx;
			atoms::y_total_external_field_array[atom] += Hy;
			atoms::z_total_external_field_array[atom] += Hz;
		}

	}

	// Add external field from thin film sample
	if(sim::ext_demag==true){

      const std::vector<double> m_l = stats::system_magnetization.get_magnetization();

		// calculate global demag field -mu_0 M D, M = m/V
		const double mu_0= -4.0*M_PI*1.0e-7/(cs::system_dimensions[0]*cs::system_dimensions[1]*cs::system_dimensions[2]*1.0e-30);
      const double HD[3]={	mu_0*sim::demag_factor[0]*m_l[0],
                           mu_0*sim::demag_factor[1]*m_l[1],
                           mu_0*sim::demag_factor[2]*m_l[2]};

		//std::cout << "mu_0" << "\t" << mu_0 << std::endl;
		//std::cout << "Magnetisation " << stats::total_mag_actual[0] << "\t" << stats::total_mag_actual[1] << "\t" << stats::total_mag_actual[2] << std::endl;
		//std::cout << "External Demag Field " << HD[0] << "\t" << HD[1] << "\t" << HD[2] << std::endl;
		for(int atom=start_index;atom<end_index;atom++){
			atoms::x_total_external_field_array[atom] += HD[0];
			atoms::y_total_external_field_array[atom] += HD[1];
			atoms::z_total_external_field_array[atom] += HD[2];
		}
	}

	return 0;
}

int calculate_thermal_fields(const int start_index,const int end_index){
	///======================================================
	/// 		Subroutine to calculate thermal fields
	///
   ///      Version 1.2 R Evans 12/08/2014
	///======================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_thermal_fields has been called" << std::endl;}

   // unroll sigma for speed
   std::vector<double> sigma_prefactor(0);
   sigma_prefactor.reserve(mp::material.size());

   // Calculate material temperature (with optional rescaling)
   for(unsigned int mat=0;mat<mp::material.size();mat++){
      double temperature = sim::temperature;
      // Check for localised temperature
      if(sim::local_temperature) temperature = mp::material[mat].temperature;
      // Calculate temperature rescaling
      double alpha = mp::material[mat].temperature_rescaling_alpha;
      double Tc = mp::material[mat].temperature_rescaling_Tc;
      // if T<Tc T/Tc = (T/Tc)^alpha else T = T
      double rescaled_temperature = temperature < Tc ? Tc*pow(temperature/Tc,alpha) : temperature;
      double sqrt_T=sqrt(rescaled_temperature);
      sigma_prefactor.push_back(sqrt_T*mp::material[mat].H_th_sigma);
   }

 	generate (atoms::x_total_external_field_array.begin()+start_index,atoms::x_total_external_field_array.begin()+end_index, mtrandom::gaussian);
	generate (atoms::y_total_external_field_array.begin()+start_index,atoms::y_total_external_field_array.begin()+end_index, mtrandom::gaussian);
	generate (atoms::z_total_external_field_array.begin()+start_index,atoms::z_total_external_field_array.begin()+end_index, mtrandom::gaussian);

	for(int atom=start_index;atom<end_index;atom++){

		const int imaterial=atoms::type_array[atom];
      const double H_th_sigma = sigma_prefactor[imaterial];

		atoms::x_total_external_field_array[atom] *= H_th_sigma;
		atoms::y_total_external_field_array[atom] *= H_th_sigma;
		atoms::z_total_external_field_array[atom] *= H_th_sigma;
	}

	return EXIT_SUCCESS;
}

int calculate_dipolar_fields(const int start_index,const int end_index){
	///======================================================
	/// 		Subroutine to calculate dipolar fields
	///
	///			Version 1.0 R Evans 02/11/2009
	///======================================================

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_dipolar_fields has been called" << std::endl;}

	// Add dipolar fields
	if(dipole::activated){
	   for(int atom=start_index;atom<end_index;atom++){
		   atoms::x_total_external_field_array[atom] += dipole::atom_dipolar_field_array_x[atom];
			atoms::y_total_external_field_array[atom] += dipole::atom_dipolar_field_array_y[atom];
			atoms::z_total_external_field_array[atom] += dipole::atom_dipolar_field_array_z[atom];
			/*std::cout << atoms::x_total_external_field_array[atom] << "\t" <<  dipole::atom_dipolar_field_array_x[atom] << "\t";
			std::cout << atoms::y_total_external_field_array[atom] << "\t" <<  dipole::atom_dipolar_field_array_y[atom] << "\t";
			std::cout << atoms::z_total_external_field_array[atom] << "\t" <<  dipole::atom_dipolar_field_array_z[atom] << std::endl;*/
      }
	}

	return 0;
}

void calculate_hamr_fields(const int start_index,const int end_index){

	if(err::check==true){std::cout << "calculate_hamr_fields has been called" << std::endl;}

	// Declare hamr variables
	const double fwhm=200.0; // A
	const double fwhm2=fwhm*fwhm;
	const double px = sim::head_position[0];
	const double py = sim::head_position[1];
	const double DeltaT=sim::Tmax-sim::Tmin;

	// declare head-field variables
	const double H_bounds_min[2]={-400.0,-250.0}; // A
	const double H_bounds_max[2]={-100.0,+250.0}; // A
	const double H_osc_freq=200.0; // A
	const double Hloc_min_x=sim::head_position[0]+H_bounds_min[0];
	const double Hloc_min_y=sim::head_position[1]+H_bounds_min[1];
	const double Hloc_max_x=sim::head_position[0]+H_bounds_max[0];
	const double Hloc_max_y=sim::head_position[1]+H_bounds_max[1];
	const double Hloc_parity_field=sim::H_applied*double(2*(int(sim::head_position[0]/H_osc_freq)%2)-1);
	const double Hvecx=sim::H_vec[0];
	const double Hvecy=sim::H_vec[1];
	const double Hvecz=sim::H_vec[2];

	// Add localised thermal field
	generate (atoms::x_total_external_field_array.begin()+start_index,atoms::x_total_external_field_array.begin()+end_index, mtrandom::gaussian);
	generate (atoms::y_total_external_field_array.begin()+start_index,atoms::y_total_external_field_array.begin()+end_index, mtrandom::gaussian);
	generate (atoms::z_total_external_field_array.begin()+start_index,atoms::z_total_external_field_array.begin()+end_index, mtrandom::gaussian);

	if(sim::head_laser_on){
		for(int atom=start_index;atom<end_index;atom++){
			const int imaterial=atoms::type_array[atom];
			const double cx = atoms::x_coord_array[atom];
			const double cy = atoms::y_coord_array[atom];
			const double r2 = (cx-px)*(cx-px)+(cy-py)*(cy-py);
			const double sqrt_T = sqrt(sim::Tmin+DeltaT*exp(-r2/fwhm2));
			const double H_th_sigma = sqrt_T*mp::material[imaterial].H_th_sigma;
			atoms::x_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
			atoms::y_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
			atoms::z_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
		}

		// Add localised applied field
		for(int atom=start_index;atom<end_index;atom++){
			const double cx = atoms::x_coord_array[atom];
			const double cy = atoms::y_coord_array[atom];
			double Hx=0.0;
			double Hy=0.0;
			double Hz=0.0;
			if((cx >= Hloc_min_x) && (cx <= Hloc_max_x) && (cy >= Hloc_min_y) && (cy <= Hloc_max_y)){
				Hx=Hvecx*Hloc_parity_field;
				Hy=Hvecy*Hloc_parity_field;
				Hz=Hvecz*Hloc_parity_field;
			}
			atoms::x_total_external_field_array[atom] += Hx;
			atoms::y_total_external_field_array[atom] += Hy;
			atoms::z_total_external_field_array[atom] += Hz;
		}
	}
	else{
		// Otherwise just use global temperature
		double sqrt_T=sqrt(sim::temperature);
		for(int atom=start_index;atom<end_index;atom++){
			const int imaterial=atoms::type_array[atom];
			const double H_th_sigma = sqrt_T*material_parameters::material[imaterial].H_th_sigma;
			atoms::x_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
			atoms::y_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
			atoms::z_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
		}
	}
}

void calculate_fmr_fields(const int start_index,const int end_index){

	if(err::check==true){std::cout << "calculate_fmr_fields has been called" << std::endl;}

	// Calculate fmr constants
	const double real_time = sim::time*mp::dt_SI;
	const double omega = sim::fmr_field_frequency*1.e9; // Hz
	const double Hfmrx = sim::fmr_field_unit_vector[0];
	const double Hfmry = sim::fmr_field_unit_vector[1];
	const double Hfmrz = sim::fmr_field_unit_vector[2];
	const double Hsinwt = sim::fmr_field_strength * sin(2.0 * M_PI * omega * real_time);
	const double Hx = Hfmrx * Hsinwt;
	const double Hy = Hfmry * Hsinwt;
	const double Hz = Hfmrz * Hsinwt;

	// Save fmr field strength for possible output
	sim::fmr_field = Hsinwt;

	if(sim::local_fmr_field==true){

		std::vector<double> H_fmr_local;
		H_fmr_local.reserve(3*mp::material.size());

		// Loop over all materials
		for(unsigned int mat=0;mat<mp::material.size();mat++){
			const double Hsinwt_local=mp::material[mat].fmr_field_strength*sin(2.0*M_PI*real_time*mp::material[mat].fmr_field_frequency);

			H_fmr_local.push_back(Hsinwt_local*mp::material[mat].fmr_field_unit_vector[0]);
			H_fmr_local.push_back(Hsinwt_local*mp::material[mat].fmr_field_unit_vector[1]);
			H_fmr_local.push_back(Hsinwt_local*mp::material[mat].fmr_field_unit_vector[2]);
		}

		// Add local field AND global field
		for(int atom=start_index;atom<end_index;atom++){
			const int imaterial=atoms::type_array[atom];
			atoms::x_total_external_field_array[atom] += Hx + H_fmr_local[3*imaterial + 0];
			atoms::y_total_external_field_array[atom] += Hy + H_fmr_local[3*imaterial + 1];
			atoms::z_total_external_field_array[atom] += Hz + H_fmr_local[3*imaterial + 2];
		}
	}
	else{
		// Add fmr field
		for(int atom=start_index;atom<end_index;atom++){
			atoms::x_total_external_field_array[atom] += Hx;
			atoms::y_total_external_field_array[atom] += Hy;
			atoms::z_total_external_field_array[atom] += Hz;
		}
	}

	return;
}

///------------------------------------------------------
///  Function to calculate LaGrange multiplier fields for
///  constrained minimization
///
///  (c) R F L Evans 2013
///
///------------------------------------------------------
void calculate_lagrange_fields(const int start_index,const int end_index){

   // LaGrange Multiplier
   const double lx=sim::lagrange_lambda_x;
   const double ly=sim::lagrange_lambda_y;
   const double lz=sim::lagrange_lambda_z;

   // Constraint vector
   const double nu_x=cos(sim::constraint_theta*M_PI/180.0)*sin(sim::constraint_phi*M_PI/180.0);
   const double nu_y=sin(sim::constraint_theta*M_PI/180.0)*sin(sim::constraint_phi*M_PI/180.0);
   const double nu_z=cos(sim::constraint_phi*M_PI/180.0);

   // Magnetisation
   const double imm=1.0/sim::lagrange_m;
   const double imm3=1.0/(sim::lagrange_m*sim::lagrange_m*sim::lagrange_m);

   const double N=sim::lagrange_N;

   // Calculate LaGrange fields
   for(int atom=start_index;atom<end_index;atom++){
      const double sx=atoms::x_spin_array[atom];
      const double sy=atoms::y_spin_array[atom];
      const double sz=atoms::z_spin_array[atom];

      //std::cout << "S " << sx << "\t" << sy << "\t" << sz << std::endl;
      //std::cout << "L " << lx << "\t" << ly << "\t" << lz << std::endl;
      //std::cout << imm << "\t" << imm3 << std::endl;

      const double lambda_dot_s = lx*sx + ly*sy + lz*sz;

      atoms::x_total_spin_field_array[atom]+=N*(lx*imm - lambda_dot_s*sx*imm3 - nu_x);
      atoms::y_total_spin_field_array[atom]+=N*(ly*imm - lambda_dot_s*sy*imm3 - nu_y);
      atoms::z_total_spin_field_array[atom]+=N*(lz*imm - lambda_dot_s*sz*imm3 - nu_z);

      //std::cout << "\t" << N*(lx*imm - lambda_dot_s*sx*imm3 - nu_x) << std::endl;
      //std::cout << "\t" << N*(ly*imm - lambda_dot_s*sy*imm3 - nu_y) << std::endl;
      //std::cout << "\t" << N*(lz*imm - lambda_dot_s*sz*imm3 - nu_z) << std::endl;
      //std::cin.get();
   }
   return;

}

//------------------------------------------------------------------------------
// Master function to calculate fields in large loop
//------------------------------------------------------------------------------
void calculate_full_spin_fields(const int start_index,const int end_index){

	using namespace sim::internal;

   for(int atom=start_index;atom<end_index;atom++){

		// temporary variables for field components
		double hx = 0.0;
		double hy = 0.0;
		double hz = 0.0;

		// temporary constant for spin components
		const double sx = atoms::x_spin_array[atom];
		const double sy = atoms::x_spin_array[atom];
		const double sz = atoms::x_spin_array[atom];

		// get material parameter
		const int material=atoms::type_array[atom];

		//----------------------------------------------------------------------------------
		// Slonczewski spin torque field
		//----------------------------------------------------------------------------------

		// save polarization to temporary constant
		const double stpx = slonczewski_spin_polarization_unit_vector[0];
		const double stpy = slonczewski_spin_polarization_unit_vector[1];
		const double stpz = slonczewski_spin_polarization_unit_vector[2];

		const double staj = slonczewski_aj[material];
		const double stbj = slonczewski_bj[material];

		// calculate field
		hx += staj*(sy*stpz - sz*stpy) + stbj*stpx;
		hy += staj*(sz*stpx - sx*stpz) + stbj*stpy;
		hz += staj*(sx*stpy - sy*stpx) + stbj*stpz;

		// save field to spin field array
		atoms::x_total_spin_field_array[atom]+=hx;
		atoms::y_total_spin_field_array[atom]+=hy;
		atoms::z_total_spin_field_array[atom]+=hz;

	}

	return;

}
