//====================================================================================================
//
//       				                    		Demag
//
//  			 					Subroutines to setup and calculate demag fields
//	 
//												Version 1.0 R Evans 18/09/2009
//
//==================================================================================================== 
#include "atoms.hpp"
#include "cells.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "demag.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

#include <cmath>
#include <iostream>

namespace demag{

	int num_demag_cells=0;
	int demag_resolution=5;
	int update_rate=1;
	int update_counter=0;

	const double prefactor=1.0e+23;

  std::vector<int> atom_demag_array(0);

  std::vector<double> demag_spin_x_array(0);
  std::vector<double> demag_spin_y_array(0);
  std::vector<double> demag_spin_z_array(0);

  std::vector<double> demag_coord_x_array(0);
  std::vector<double> demag_coord_y_array(0);
  std::vector<double> demag_coord_z_array(0);

  std::vector<double> demag_field_x_array(0);
  std::vector<double> demag_field_y_array(0);
  std::vector<double> demag_field_z_array(0);

	double** rij_matrix=NULL;

	bool demag_set=false;

}

int set_demag(){
	//--------------------------------------------------------------------
	//
	//											set_demag
	//
	//				Function to initialise demag arrays and variables
	//
	//							R F Evans 21/09/09
	//	
	//--------------------------------------------------------------------

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){
		std::cout << "set_demag has been called " << vmpi::my_rank << std::endl;
	}

	//----------------------------------------------------------
	// If already initialised print warning and do nothing
	//----------------------------------------------------------
	if(demag::demag_set==true){
		std::cerr << "Warning, demag already initialised, ignoring re-initialisation" << std::endl; 
	}

	//----------------------------------------------------------
	// If not initialised then set up arrays
	//----------------------------------------------------------
	if(demag::demag_set==false){
		vinfo << "=========================================================================" << std::endl;
		vinfo << " Setting up Demag Fields using Resolution = " << demag::demag_resolution << std::endl;
		vinfo << "=========================================================================" << std::endl;
		//----------------------------------------------------------
		// Calculate number of demag cells
		//----------------------------------------------------------
		//std::cout << material_parameters::int_system_dimensions[0] << "\t" << material_parameters::int_system_dimensions[1] << "\t" << material_parameters::int_system_dimensions[2] << std::endl;
		if(demag::demag_resolution>0){
			int nx_cell = ceil(double(material_parameters::int_system_dimensions[0])/double(2*demag::demag_resolution));
			int ny_cell = ceil(double(material_parameters::int_system_dimensions[1])/double(6*demag::demag_resolution));
			int nz_cell = ceil(double(material_parameters::int_system_dimensions[2])/double(2*demag::demag_resolution));
			demag::num_demag_cells=nx_cell*ny_cell*nz_cell;
			vinfo << "Number of demag cells (nx,ny,nz,total)\t" << nx_cell << "\t" << ny_cell << "\t"; 
			vinfo << nz_cell << "\t" << demag::num_demag_cells << std::endl;
		}
		else{
			demag::num_demag_cells=atoms::num_atoms; //check correct for mpi?
			vinfo << "Number of demag cells (= num_atoms)\t" << demag::num_demag_cells << std::endl;
		}
		//----------------------------------------------------------
		// Resize arrays for demag calculation
		//----------------------------------------------------------
		demag::atom_demag_array.resize(atoms::num_atoms);

		demag::demag_spin_x_array.resize(demag::num_demag_cells,0.0);
		demag::demag_spin_y_array.resize(demag::num_demag_cells,0.0);
		demag::demag_spin_z_array.resize(demag::num_demag_cells,0.0);

		demag::demag_coord_x_array.resize(demag::num_demag_cells,0.0);
		demag::demag_coord_y_array.resize(demag::num_demag_cells,0.0);
		demag::demag_coord_z_array.resize(demag::num_demag_cells,0.0);

		demag::demag_field_x_array.resize(demag::num_demag_cells,0.0);
		demag::demag_field_y_array.resize(demag::num_demag_cells,0.0);
		demag::demag_field_z_array.resize(demag::num_demag_cells,0.0);

		atoms::x_dipolar_field_array.resize(atoms::num_atoms,0.0);
		atoms::y_dipolar_field_array.resize(atoms::num_atoms,0.0);
		atoms::z_dipolar_field_array.resize(atoms::num_atoms,0.0);

		demag::update_counter=0;	// Always recalculate first time

		//-----------------------------------------------------------------------------------------
		// Find which atoms are in which demag cell and calculate mean coordinates for each cell
		//-----------------------------------------------------------------------------------------
		if(demag::demag_resolution>0){
			const int nx_cell = ceil(double(material_parameters::int_system_dimensions[0])/double(2*demag::demag_resolution));
			const int ny_cell = ceil(double(material_parameters::int_system_dimensions[1])/double(6*demag::demag_resolution));
			const int nz_cell = ceil(double(material_parameters::int_system_dimensions[2])/double(2*demag::demag_resolution));
			// Allocate 3D demag array
			int*** three_D_coord_array;
			try{three_D_coord_array=new int**[nx_cell];
				for(int i=0; i<nx_cell ; i++){
					three_D_coord_array[i]=new int*[ny_cell];
					for(int j=0; j<ny_cell ; j++){
						three_D_coord_array[i][j]=new int[nz_cell];
					}
				}
			}
			catch(...){std::cout << "error allocating 3D coord array in demag setup" << std::endl;exit(1);}
	
			// Populate cell numbers
			int cell=0;
			for(int i=0; i<nx_cell; i++){
				for(int j=0; j<ny_cell ; j++){
					for(int k=0; k<nz_cell ; k++){ 
						three_D_coord_array[i][j][k] = cell;
						cell++;
					}
				}
			}
			if(cell!=demag::num_demag_cells){
				std::cerr << "error numbering demag cells, exiting" << std::endl;
				exit(1);
			} 
			for(int atom=0;atom<atoms::num_atoms;atom++){
				int cx=int(double(atoms::x_coord_array[atom])/double(2*demag::demag_resolution));
				int cy=int(double(atoms::y_coord_array[atom])/double(6*demag::demag_resolution));
				int cz=int(double(atoms::z_coord_array[atom])/double(2*demag::demag_resolution));
				demag::atom_demag_array[atom]=three_D_coord_array[cx][cy][cz];
			}
			// Deallocate space for 3D atom coord array
			try{for(int i=0; i<nx_cell ; i++){
					for(int j=0; j<ny_cell ;j++){
						delete [] three_D_coord_array[i][j];
					}
					delete [] three_D_coord_array[i];
				}
				delete [] three_D_coord_array;
				three_D_coord_array=NULL;
			}
			catch(...){std::cerr << "error deallocating 3D coord array" << std::endl; exit(1);}
			// Calculate number of spins per cell (initialised to zero)
			std::vector<int> num_spins_per_cell(0,demag::num_demag_cells);
			for(int atom=0;atom<atoms::num_atoms;atom++){
				const int cell = demag::atom_demag_array[atom];
				num_spins_per_cell[cell]++;
				// Calculate sum of coordinates for each cell
				demag::demag_coord_x_array[cell]+=double(atoms::x_coord_array[atom]);
				demag::demag_coord_y_array[cell]+=double(atoms::y_coord_array[atom]);
				demag::demag_coord_z_array[cell]+=double(atoms::z_coord_array[atom]);
			}
			// Calculate mean coordinates and output cell data to vinfo
			vinfo << "=========================================================================" << std::endl;
			vinfo << "Number of atoms/cell: cell number, num atoms, coord" << std::endl;
			vinfo << "=========================================================================" << std::endl;
			for(int cell=0;cell<demag::num_demag_cells;cell++){
				demag::demag_coord_x_array[cell]*=material_parameters::lattice_space_conversion[0];
				demag::demag_coord_y_array[cell]*=material_parameters::lattice_space_conversion[1];
				demag::demag_coord_z_array[cell]*=material_parameters::lattice_space_conversion[2];
				demag::demag_coord_x_array[cell]/=double(num_spins_per_cell[cell]);
				demag::demag_coord_y_array[cell]/=double(num_spins_per_cell[cell]);
				demag::demag_coord_z_array[cell]/=double(num_spins_per_cell[cell]);
				vinfo << cell << "\t" << num_spins_per_cell[cell] << "\t";
				vinfo << demag::demag_coord_x_array[cell] << "\t" << demag::demag_coord_y_array[cell];
				vinfo << "\t" << demag::demag_coord_z_array[cell] << "\t" << std::endl;
			}
		}
		else{
			for(int atom=0;atom<atoms::num_atoms;atom++){
				demag::atom_demag_array[atom]=atom; //check correct for mpi?
				demag::demag_coord_x_array[atom]=double(atoms::x_coord_array[atom])*material_parameters::lattice_space_conversion[0];
				demag::demag_coord_y_array[atom]=double(atoms::y_coord_array[atom])*material_parameters::lattice_space_conversion[1];
				demag::demag_coord_z_array[atom]=double(atoms::z_coord_array[atom])*material_parameters::lattice_space_conversion[2];
			}
		}
		//         mu_o       3.(m.r_hat)r_hat - m                  1.0 e-7
		// H = ---------- . ------------------------,   prefactor = ------- = 1.0 e+23
		//      4*pi a^3              |r|^3                        1.0 e-30
		//
		// Requires m = SUM(S.mu_s)
		//------------------------------------------------------------------
		// Attempt allocation of rij_matrix 
		//------------------------------------------------------------------
		vinfo << "=========================================================================" << std::endl;
		vinfo << "  Estimated memory requirements for rij_matrix: " << double(sizeof(double))*double(demag::num_demag_cells)*double(demag::num_demag_cells)/1.0e6 << " MB" << std::endl;
		bool super_allocation_success=true;
		bool sub_allocation_success=true;
		// Attempt allocation of array of pointers
		try{demag::rij_matrix=new double*[demag::num_demag_cells];}
		catch(...){
			super_allocation_success=false;
		}
		if(super_allocation_success==true){
			// Attempt allocation of internal array of pointers
			for(int i=0;i<demag::num_demag_cells;i++){
				try{demag::rij_matrix[i]=new double[demag::num_demag_cells];}
				catch(...){
					sub_allocation_success=false;
					break;
				}
			}
		}
		else{
			std::cerr << "Warning - insufficient memory for rij_matrix for demag fields" << std::endl;
			vinfo << "  Insufficient memory for rij matrix. This will severely impact performance of the code." << std::endl;
			vinfo << "  Try increasing the demag supercell size from " << demag::demag_resolution << " to reduce memory requirements." << std::endl;
			try{delete[] demag::rij_matrix; demag::rij_matrix=NULL;}
			catch(...){std::cerr << "Error deallocating rij_matrix, exiting" << std::endl; exit(1);}
		}
		
		// If allocation failed print error message and clean up memory
		if(sub_allocation_success==false){
			std::cerr << "Warning - insufficient memory for rij_matrix for demag fields" << std::endl;
			vinfo << "Insufficient memory for rij matrix. This will severely impact performance of the code." << std::endl;
			vinfo << "Try increasing the demag supercell size from " << demag::demag_resolution << " to reduce memory requirements." << std::endl;
			for(int i=0;i<demag::num_demag_cells;i++){
				try{if(demag::rij_matrix[i]!=NULL) delete[] demag::rij_matrix[i]; demag::rij_matrix[i]=NULL;}
				catch(...){std::cerr << "Error deallocating rij_matrix, exiting" << std::endl; exit(1);}
			}
			try{delete[] demag::rij_matrix; demag::rij_matrix=NULL;}
			catch(...){std::cerr << "Error deallocating rij_matrix, exiting" << std::endl; exit(1);}
		}
		else{vinfo << "  rij_matrix allocation successful" << std::endl;
			//------------------------------------------------------------------
			// Calculate rij matrix (demag cell cordinates are in angstroms)
			//------------------------------------------------------------------
			std::cout << "Calculating rij_matrix for Dipolar fields" << std::flush;
			for(int i=0;i<demag::num_demag_cells;i++){
				if(i%(demag::num_demag_cells/20)==0){
					std::cout << "."<< std::flush;
				}
				for(int j=0;j<demag::num_demag_cells;j++){
					if(i!=j){
						double dx = demag::demag_coord_x_array[j]-demag::demag_coord_x_array[i];
						double dy = demag::demag_coord_y_array[j]-demag::demag_coord_y_array[i];
						double dz = demag::demag_coord_z_array[j]-demag::demag_coord_z_array[i];
						double rij = sqrt(dx*dx+dy*dy+dz*dz);
						demag::rij_matrix[i][j]=1.0/rij;
					}
				}
			}
			std::cout << "Done!"<< std::endl;
		}

		vinfo << "=========================================================================" << std::endl;
		vinfo << "  Demag set-up completed" << std::endl;
		vinfo << "=========================================================================" << std::endl;

	//----------------------------------------------------------
	// End of if
	//----------------------------------------------------------
		demag::demag_set=true;
	}

	return 0;

}

int demag_field_update(){
	//--------------------------------------------------------------------
	//
	//								demag_field_update
	//
	//						Function to update demag fields
	//
	//							R F Evans 02/11/09
	//	
	//--------------------------------------------------------------------

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){
		std::cout << "demag_field_update has been called " << vmpi::my_rank << std::endl;
	}

	//----------------------------------------------------------
	// If not initialised print warning and do so
	//----------------------------------------------------------
	if(demag::demag_set==false){
		std::cerr << "Warning, demag not initialised - initialising" << std::endl;
		set_demag();
	}

	//----------------------------------------------------------
	// Update Dipolar Field Array
	//----------------------------------------------------------
  // Field has units ((T^2 m^2) / (N m^3)) * (J/T) = T
	for(int i=0;i<cells::num_cells;i++){
    // zero field arrays
    cells::x_field_array[i]=0.0;
    cells::y_field_array[i]=0.0;
    cells::z_field_array[i]=0.0;
		for(int j=0;j<cells::num_cells;j++){
			if(i!=j){
				
        const double mx = cells::x_mag_array[j];
				const double my = cells::y_mag_array[j];
				const double mz = cells::z_mag_array[j];

				const double dx = cells::x_coord_array[j]-cells::x_coord_array[i];
				const double dy = cells::y_coord_array[j]-cells::y_coord_array[i];
				const double dz = cells::z_coord_array[j]-cells::z_coord_array[i];

				const double drij = 1.0/sqrt(dx*dx+dy*dy+dz*dz);
        const double drij3 = 3.0*drij;

				const double ex = dx*drij;
				const double ey = dy*drij;
				const double ez = dz*drij;

        const double s_dot_e = (mx * ex + my * ey + mz * ez);

				cells::x_field_array[i]+=(3.0 * s_dot_e * ex - mx)*drij3;
				cells::y_field_array[i]+=(3.0 * s_dot_e * ey - my)*drij3;
				cells::z_field_array[i]+=(3.0 * s_dot_e * ez - mz)*drij3;
			}
		}
		cells::x_field_array[i]*=demag::prefactor;
		cells::y_field_array[i]*=demag::prefactor;
		cells::z_field_array[i]*=demag::prefactor;
		
    //vdp << "\t" << demag::demag_coord_x_array[i] << "\t" << demag::demag_coord_y_array[i] << "\t" << demag::demag_coord_z_array[i] << "\t";
		//vdp << demag::demag_field_x_array[i] << "\t" << demag::demag_field_y_array[i] << "\t" << demag::demag_field_z_array[i] << std::endl;
	}
	//----------------------------------------------------------
	// Update Atomistic Dipolar Field Array
	//----------------------------------------------------------
	for(int atom=0;atom<atoms::num_atoms;atom++){
		const int cell = atoms::cell_array[atom];

		// Copy field from macrocell to atomistic spin
		atoms::x_dipolar_field_array[atom]=cells::x_field_array[cell];
		atoms::y_dipolar_field_array[atom]=cells::y_field_array[cell];
		atoms::z_dipolar_field_array[atom]=cells::z_field_array[cell];
	}
	return 0;

}
