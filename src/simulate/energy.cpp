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
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///

// Standard Libraries
#include <algorithm>
#include <cmath>
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "demag.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

namespace sim{

/// @brief Calculates the exchange energy for a single spin (isotropic).
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @param[in] int atom number 
/// @param[in] imaterial material of local atom
/// @param[in] Sx x-spin of local atom  
/// @param[in] Sy y-spin of local atom 
/// @param[in] Sz z-spin of local atom 
/// @return exchange energy
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
inline double spin_exchange_energy_isotropic(const int atom, const int imaterial, const double Sx, const double Sy, const double Sz){
	
	// energy
	double energy=0.0;
	
	// Loop over neighbouring spins to calculate exchange
	for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){
			
		const int natom = atoms::neighbour_list_array[nn];
		const int jmaterial=atoms::type_array[natom]; 
		const double Jij = material_parameters::material[imaterial].Jij_matrix[jmaterial];

		energy+=Jij*(atoms::x_spin_array[natom]*Sx + atoms::y_spin_array[natom]*Sy + atoms::z_spin_array[natom]*Sz);
	}
		
	return energy;
	
}

/// @brief Calculates the exchange energy for a single spin (anisotropic).
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
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
/// @return exchange energy
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
inline double spin_exchange_energy_anisotropic(const int atom, const int imaterial, const double Sx, const double Sy, const double Sz){
	
	// energy
	double energy=0.0;
	
	// Loop over neighbouring spins to calculate exchange
	for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){
			
		const int natom = atoms::neighbour_list_array[nn];
		const int jmaterial=atoms::type_array[natom]; 
		const double Jij_x = material_parameters::material[imaterial].Jij_matrix[jmaterial];
		const double Jij_y = material_parameters::material[imaterial].Jij_matrix[jmaterial];
		const double Jij_z = material_parameters::material[imaterial].Jij_matrix[jmaterial];

		energy+=(Jij_x*atoms::x_spin_array[natom]*Sx + Jij_y*atoms::y_spin_array[natom]*Sy + Jij_z*atoms::z_spin_array[natom]*Sz);
	}
		
	return energy;
	
}

/// @brief Calculates the exchange energy for a single spin (matrix).
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
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
/// @return exchange energy
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
inline double spin_exchange_energy_matrix(const int atom, const int imaterial, const double Sx, const double Sy, const double Sz){
	
	// energy
	double energy=0.0;
	
	// Loop over neighbouring spins to calculate exchange
	for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){
			
		const int natom = atoms::neighbour_list_array[nn];
		const int jmaterial=atoms::type_array[natom]; 
		const double Jij_xx = material_parameters::material[imaterial].Jij_matrix[jmaterial];
		const double Jij_xy = 0.0; //material_parameters::material[imaterial].Jij_matrix[jmaterial];
		const double Jij_xz = 0.0; //material_parameters::material[imaterial].Jij_matrix[jmaterial];
		const double Jij_yx = 0.0; //material_parameters::material[imaterial].Jij_matrix[jmaterial];
		const double Jij_yy = material_parameters::material[imaterial].Jij_matrix[jmaterial];
		const double Jij_yz = 0.0; //material_parameters::material[imaterial].Jij_matrix[jmaterial];
		const double Jij_zx = 0.0; //material_parameters::material[imaterial].Jij_matrix[jmaterial];
		const double Jij_zy = 0.0; //material_parameters::material[imaterial].Jij_matrix[jmaterial];
		const double Jij_zz = material_parameters::material[imaterial].Jij_matrix[jmaterial];

		energy+=(Jij_xx*atoms::x_spin_array[natom]*Sx + Jij_xy*atoms::y_spin_array[natom]*Sx + Jij_xz*atoms::z_spin_array[natom]*Sx + 
					Jij_yx*atoms::x_spin_array[natom]*Sy + Jij_yy*atoms::y_spin_array[natom]*Sy + Jij_yz*atoms::z_spin_array[natom]*Sy + 
					Jij_zx*atoms::x_spin_array[natom]*Sz + Jij_zy*atoms::y_spin_array[natom]*Sz + Jij_zz*atoms::z_spin_array[natom]*Sz);
	}
		
	return energy;
	
}

/// @brief Calculates the uniaxial anisotropy energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
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
/// @return exchange energy
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
inline double spin_uniaxial_energy(const int atom, const int imaterial, const double Sx, const double Sy, const double Sz){
	
	return mp::material[imaterial].Ku*Sz*Sz;
	
}

/// @brief Calculates the applied field energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
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
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
inline double spin_applied_field_energy(const int atom, const int imaterial, const double Sx, const double Sy, const double Sz){;

	return -sim::H_applied*(sim::H_vec[0]*Sx + sim::H_vec[1]*Sy + sim::H_vec[2]*Sz);

}

/// @brief Calculates the surface anisotropy energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
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
/// @return exchange energy
///
/// @internal
///	Created:		13/09/2011
///	Revision:	  ---
///=====================================================================================
///
inline double spin_surface_anisotropy_energy(const int atom, const int imaterial, const double Sx, const double Sy, const double Sz){
	
	double energy=0.0;

	if(atoms::surface_array[atom]==true && sim::surface_anisotropy==true){
		const double Ks=mp::material[imaterial].Ks;
		for(int nn=atoms::nearest_neighbour_list_si[atom];nn<atoms::nearest_neighbour_list_ei[atom];nn++){
			const double si_dot_eij=(Sx*atoms::eijx[nn]+Sy*atoms::eijy[nn]+Sz*atoms::eijz[nn]);
			energy+=Ks*si_dot_eij*si_dot_eij;
		}
	}
	
	return energy;
}

/// @brief Calculates the magnetostatic energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2012. All Rights Reserved.
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
///	Created:		08/06/2012
///	Revision:	  ---
///=====================================================================================
///
inline double spin_magnetostatic_energy(const int atom, const int imaterial, const double Sx, const double Sy, const double Sz){
	
	return -1.0*(atoms::x_dipolar_field_array[atom]*Sx+atoms::y_dipolar_field_array[atom]*Sy+atoms::z_dipolar_field_array[atom]*Sz);
}

/// @brief Calculates the total energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
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
///	Created:		07/02/2011
///	Revision:	  ---
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
	energy+=spin_exchange_energy_isotropic(atom, imaterial, Sx, Sy, Sz);
	energy+=spin_applied_field_energy(atom, imaterial, Sx, Sy, Sz);
	energy+=spin_uniaxial_energy(atom, imaterial, Sx, Sy, Sz);
	energy+=spin_surface_anisotropy_energy(atom, imaterial, Sx, Sy, Sz);
	energy+=spin_magnetostatic_energy(atom, imaterial, Sx, Sy, Sz);
	
	return energy; // Tesla
}

} // end of namespace sim

