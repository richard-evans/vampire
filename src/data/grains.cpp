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
#include "atoms.hpp"  
#include "grains.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

#include <cmath>
#include <iostream>

namespace grains{

	int num_grains=1; // always assume 1 grain 

	std::vector <int> grain_size_array(0);

	std::vector <double> x_coord_array(0);
	std::vector <double> y_coord_array(0);
	std::vector <double> z_coord_array(0);

	std::vector <double> x_mag_array(0);
	std::vector <double> y_mag_array(0);
	std::vector <double> z_mag_array(0);
	std::vector <double> mag_m_array(0);
	std::vector <double> sat_mag_array(0);

	std::vector <double> x_mat_mag_array(0);
	std::vector <double> y_mat_mag_array(0);
	std::vector <double> z_mat_mag_array(0);
	std::vector <double> mat_mag_m_array(0);
	std::vector <double> mat_sat_mag_array(0);


int set_properties(){
	//========================================================================================================
	///		 				Function to calculate number of atoms in each grain
	//
	///														Version 1.0
	//
	///												R F Evans 15/07/2009
	//
	//========================================================================================================

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "grains::set_properties has been called" << std::endl;}	
	#ifdef MPICF
		const unsigned int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
	#else
		const unsigned int num_local_atoms = atoms::num_atoms;
	#endif

	// resize and initialise grain arrays
	if(grains::num_grains > 0){
		grains::grain_size_array.resize(grains::num_grains,0);
		grains::x_coord_array.resize(grains::num_grains,0.0);
		grains::y_coord_array.resize(grains::num_grains,0.0);
		grains::z_coord_array.resize(grains::num_grains,0.0);
		grains::x_mag_array.resize(grains::num_grains,0.0);
		grains::y_mag_array.resize(grains::num_grains,0.0);
		grains::z_mag_array.resize(grains::num_grains,0.0);
		grains::mag_m_array.resize(grains::num_grains,0.0);
		grains::sat_mag_array.resize(grains::num_grains,0.0);
		
		if(mp::num_materials>1){
			grains::x_mat_mag_array.resize(grains::num_grains*mp::num_materials,0.0);
			grains::y_mat_mag_array.resize(grains::num_grains*mp::num_materials,0.0);
			grains::z_mat_mag_array.resize(grains::num_grains*mp::num_materials,0.0);
			grains::mat_mag_m_array.resize(grains::num_grains*mp::num_materials,0.0);
			grains::mat_sat_mag_array.resize(grains::num_grains*mp::num_materials,0.0);
		}
	}
	else{
		terminaltextcolor(RED);
		std::cerr << "Warning - no grains detected!" << std::endl;
		terminaltextcolor(WHITE);
	}

	// loop over atoms to determine grain properties
	for(unsigned int atom=0;atom< num_local_atoms;atom++){
		const int grain = atoms::grain_array[atom];
		const int mat = atoms::type_array[atom];

		// check grain is within allowable bounds
		if((grain>=0) && (grain<grains::num_grains)){
			grains::grain_size_array[grain]+=1;
			grains::x_coord_array[grain]+=atoms::x_coord_array[atom];
			grains::y_coord_array[grain]+=atoms::y_coord_array[atom];
			grains::z_coord_array[grain]+=atoms::z_coord_array[atom];
			grains::sat_mag_array[grain]+=mp::material[mat].mu_s_SI;
			if(mp::num_materials>1) grains::mat_sat_mag_array[grain*mp::num_materials+mat]+=mp::material[mat].mu_s_SI;
		}
		else{
			terminaltextcolor(RED);
			std::cerr << "Error - atom " << atom << " belongs to grain " << grain << " which is greater than maximum number of grains ";
			std::cerr << grains::num_grains << ". Exiting" << std::endl;
			terminaltextcolor(WHITE);
			err::vexit();
		}
	}

	// Reduce grain properties on all CPUs
	#ifdef MPICF
		//MPI::COMM_WORLD.Allreduce(&grains::grain_size_array[0], &grains::grain_size_array[0],grains::num_grains, MPI_INT,MPI_SUM);
		//MPI::COMM_WORLD.Allreduce(&grains::x_coord_array[0], &grains::x_coord_array[0],grains::num_grains, MPI_INT,MPI_SUM);
		//MPI::COMM_WORLD.Allreduce(&grains::y_coord_array[0], &grains::y_coord_array[0],grains::num_grains, MPI_INT,MPI_SUM);
		//MPI::COMM_WORLD.Allreduce(&grains::z_coord_array[0], &grains::z_coord_array[0],grains::num_grains, MPI_INT,MPI_SUM);
		//MPI::COMM_WORLD.Allreduce(&grains::sat_mag_array[0], &grains::sat_mag_array[0],grains::num_grains, MPI_INT,MPI_SUM);
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &grains::grain_size_array[0],grains::num_grains, MPI_INT,MPI_SUM);
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &grains::x_coord_array[0],grains::num_grains, MPI_DOUBLE,MPI_SUM);
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &grains::y_coord_array[0],grains::num_grains, MPI_DOUBLE,MPI_SUM);
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &grains::z_coord_array[0],grains::num_grains, MPI_DOUBLE,MPI_SUM);
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &grains::sat_mag_array[0],grains::num_grains, MPI_DOUBLE,MPI_SUM);
		if(mp::num_materials>1) MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &grains::mat_sat_mag_array[0],grains::num_grains*mp::num_materials, MPI_DOUBLE,MPI_SUM);
	#endif

	//vinfo << "-------------------------------------------------------------------------------------------------------------------" << std::endl;
	//vinfo << "# Grain number\tnum atoms\tx\ty\tz\t" << std::endl;
	// Calculate mean grains coordinates 
	for(int grain=0;grain<grains::num_grains;grain++){
		// check for grains with zero atoms
		if(grains::grain_size_array[grain]==0){
			//std::cerr << "Warning Grain " << grain << " has no constituent atoms!" << std::endl;
		}
		else{
			grains::x_coord_array[grain]/=double(grains::grain_size_array[grain]);
			grains::y_coord_array[grain]/=double(grains::grain_size_array[grain]);
			grains::z_coord_array[grain]/=double(grains::grain_size_array[grain]);
			// output grain properties to vinfo
			//vinfo << grain << "\t" << grains::grain_size_array[grain] << "\t" << grains::x_coord_array[grain];
			//vinfo << "\t" << grains::y_coord_array[grain] << "\t" << grains::z_coord_array[grain] << std::endl;
		}
	}

	return EXIT_SUCCESS;

}

int mag(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "grains::mag has been called" << std::endl;}

	#ifdef MPICF
		const unsigned int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
	#else
		const unsigned int num_local_atoms = atoms::num_atoms;
	#endif

	//initialise magnetisations
	for(int grain=0;grain<grains::num_grains;grain++){
		grains::x_mag_array[grain]=0.0;
		grains::y_mag_array[grain]=0.0;
		grains::z_mag_array[grain]=0.0;
		grains::mag_m_array[grain]=0.0;
		if(mp::num_materials>1){
			// loop over materials
			for(int mat=0;mat<mp::num_materials;mat++){
				grains::x_mat_mag_array[grain*mp::num_materials+mat]=0.0;
				grains::y_mat_mag_array[grain*mp::num_materials+mat]=0.0;
				grains::z_mat_mag_array[grain*mp::num_materials+mat]=0.0;
				grains::mat_mag_m_array[grain*mp::num_materials+mat]=0.0;
			}
		}
	}

	// function to calculate grain magnetisations
	for(unsigned int atom=0;atom< num_local_atoms;atom++){

		const int grain = atoms::grain_array[atom];
		const int mat = atoms::type_array[atom];

		// check grain is within allowable bounds
		if((grain>=0) && (grain<grains::num_grains)){
			grains::x_mag_array[grain]+=(atoms::x_spin_array[atom]*mp::material[mat].mu_s_SI);
			grains::y_mag_array[grain]+=(atoms::y_spin_array[atom]*mp::material[mat].mu_s_SI);
			grains::z_mag_array[grain]+=(atoms::z_spin_array[atom]*mp::material[mat].mu_s_SI);
			if(mp::num_materials>1){
				grains::x_mat_mag_array[grain*mp::num_materials+mat]+=(atoms::x_spin_array[atom]*mp::material[mat].mu_s_SI);
				grains::y_mat_mag_array[grain*mp::num_materials+mat]+=(atoms::y_spin_array[atom]*mp::material[mat].mu_s_SI);
				grains::z_mat_mag_array[grain*mp::num_materials+mat]+=(atoms::z_spin_array[atom]*mp::material[mat].mu_s_SI);
			}
			
		}
		else{
			terminaltextcolor(RED);
			std::cerr << "Error - atom " << atom << " belongs to grain " << grain << " which is greater than maximum number of grains ";
			std::cerr << grains::num_grains << ". Exiting" << std::endl;
			terminaltextcolor(WHITE);
			err::vexit();
		}
	}

	// Reduce grain properties on all CPUs
	#ifdef MPICF
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &grains::x_mag_array[0],grains::num_grains, MPI_DOUBLE,MPI_SUM);
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &grains::y_mag_array[0],grains::num_grains, MPI_DOUBLE,MPI_SUM);
		MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &grains::z_mag_array[0],grains::num_grains, MPI_DOUBLE,MPI_SUM);
		if(mp::num_materials>1) MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &grains::x_mat_mag_array[0],grains::num_grains*mp::num_materials, MPI_DOUBLE,MPI_SUM);
		if(mp::num_materials>1) MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &grains::y_mat_mag_array[0],grains::num_grains*mp::num_materials, MPI_DOUBLE,MPI_SUM);
		if(mp::num_materials>1) MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &grains::z_mat_mag_array[0],grains::num_grains*mp::num_materials, MPI_DOUBLE,MPI_SUM);
	#endif

	// calculate mag_m of each grain and normalised direction
	for(int grain=0;grain<grains::num_grains;grain++){
		// check for grains with zero atoms
		if(grains::grain_size_array[grain]==0){
			//std::cerr << "Warning Grain " << grain << " has no constituent atoms!" << std::endl;
		}
		else{
			double norm = 1.0/grains::sat_mag_array[grain];
			double mx = grains::x_mag_array[grain]*norm;
			double my = grains::y_mag_array[grain]*norm;
			double mz = grains::z_mag_array[grain]*norm;
			grains::mag_m_array[grain] = sqrt(mx*mx+my*my+mz*mz);
			// Normalise to unit magnetisation
			grains::x_mag_array[grain]*=norm;
			grains::y_mag_array[grain]*=norm;
			grains::z_mag_array[grain]*=norm;
			// loop over all materials and normalise
			if(mp::num_materials>1){
				for(int mat=0;mat<mp::num_materials;mat++){
					const unsigned int idx=grain*mp::num_materials+mat;
					const double immm=1.0/grains::mat_sat_mag_array[idx];
					double mx = grains::x_mat_mag_array[idx]*immm;
					double my = grains::y_mat_mag_array[idx]*immm;
					double mz = grains::z_mat_mag_array[idx]*immm;
					grains::mat_mag_m_array[idx] = sqrt(mx*mx+my*my+mz*mz);
					// Normalise to unit magnetisation
					grains::x_mat_mag_array[idx]*=immm;
					grains::y_mat_mag_array[idx]*=immm;
					grains::z_mat_mag_array[idx]*=immm;
				}
			}
		}
	}
	
	return EXIT_SUCCESS;

}

int output_mag(std::ofstream& ofile){
	// check calling of routine if error checking is activated
        if(err::check==true){std::cout << "grains::output_mag has been called" << std::endl;}

	// calculate grain magnetisations
	grains::mag();

		if(vmpi::my_rank==0){
		// check file stream is open
		if(!ofile.is_open()){
			terminaltextcolor(RED);
			std::cerr << "Error - file stream is unopened for grains::output_mag(), exiting!" << std::endl;
			terminaltextcolor(WHITE);
		err::vexit();
		}

		unsigned int id=0; // grain id (excluding grains with zero atoms)

		// calculate mag_m of each grain and normalised direction
		for(int grain=0;grain<grains::num_grains;grain++){
			// check for grains with zero atoms
			if(grains::grain_size_array[grain]==0){
				//std::cerr << "Warning Grain " << grain << " has no constituent atoms!" << std::endl;
			}
			else{
				ofile << grain << "\t" << id << "\t" << grains::x_mag_array[grain] << "\t" << grains::y_mag_array[grain];
				ofile << "\t" << grains::z_mag_array[grain] << "\t" << grains::mag_m_array[grain] << std::endl;
				id++;
			}
		}
	}

	return EXIT_SUCCESS;
}

void output_mat_mag(std::ostream& stream){
	// loop over all grains
	for(int grain=0;grain<grains::num_grains;grain++){
		// check for grains with zero atoms
		if(grains::grain_size_array[grain]!=0){
			if(mp::num_materials>1){
				// loop over all materials
				for(int mat=0;mat<mp::num_materials;mat++){
					const unsigned int idx=grain*mp::num_materials+mat;
					const double imagm = 1.0/grains::mat_mag_m_array[idx];
					//stream << grains::grain_size_array[grain] << "\t";
					//stream << grains::mat_sat_mag_array[idx] << "\t";
					stream << grains::x_mat_mag_array[idx]*imagm << "\t";
					stream << grains::y_mat_mag_array[idx]*imagm << "\t";
					stream << grains::z_mat_mag_array[idx]*imagm << "\t";
					stream << grains::mat_mag_m_array[idx] << "\t";
				}
			}
			else{
            const double imagm = 1.0/grains::mag_m_array[grain];
            stream << grains::x_mag_array[grain]*imagm << "\t";
            stream << grains::y_mag_array[grain]*imagm << "\t";
            stream << grains::z_mag_array[grain]*imagm << "\t";
            stream << grains::mag_m_array[grain] << "\t";
			}
		}
	}
}

} // End of namespace grains
