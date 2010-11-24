#include "atoms.hpp"  
#include "cells.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace cells{
	
	int num_cells=0;
	int size;
	int update_rate;
	int update_counter;

	std::vector <int> num_atoms_in_cell;

	std::vector <double> x_coord_array;
	std::vector <double> y_coord_array;
	std::vector <double> z_coord_array;

	std::vector <double> x_mag_array;
	std::vector <double> y_mag_array;
	std::vector <double> z_mag_array;
	
	std::vector <double> x_field_array;
	std::vector <double> y_field_array;
	std::vector <double> z_field_array;

	// Function to initialise cells
	int initialise(){
		// check calling of routine if error checking is activated
		if(err::check==true) std::cout << "cells::initialise has been called" << std::endl;
		
		// set initial cell variables
		cells::size=5; // In units of a
		cells::update_rate=10; // In timesteps
		cells::num_cells=0;
		cells::update_counter=0;
		
		// determine number of cells in each direction
		unsigned int ncellx = ceil((mp::system_dimensions[0]/mp::lattice_constant[0])/double(cells::size));
		unsigned int ncelly = ceil((mp::system_dimensions[1]/mp::lattice_constant[1])/double(cells::size));
		unsigned int ncellz = ceil((mp::system_dimensions[2]/mp::lattice_constant[2])/double(cells::size));
		
		//update total number of cells
		cells::num_cells=ncellx*ncelly*ncellz;
		
		std::cout << "Cells in x,y,z: " << ncellx << "\t" << ncelly << "\t" << ncellz << std::endl;
		std::cout << "Total number of cells: " << cells::num_cells << std::endl;
		
		// Determine number of cells in x,y,z
		const int d[3]={ncellx,ncelly,ncellz};
	
		// Set cell counter
		int cell=0;
		
		// Declare array for create space for 3D supercell array
		int*** supercell_array;
		std::cout << "Memory required for cell list calculation:" << 8.0*double(d[0])*double(d[1])*double(d[2])/1.0e6 << " MB" << std::endl;
		try{supercell_array=new int**[d[0]];
			for(int i=0; i<d[0] ; i++){
				supercell_array[i]=new int*[d[1]];
				for(int j=0; j<d[1] ; j++){
					supercell_array[i][j]=new int[d[2]];
					for(int k=0; k<d[2] ; k++){
						supercell_array[i][j][k]=cell;
						cell++;
					}
				}
			}
		}
		catch(...){std::cerr << "Error allocating supercell_array for cell list calculation" << std::endl;exit(EXIT_FAILURE);}
		
		// offset cells to prevent rounding error
		double atom_offset[3]={0.25*mp::lattice_constant[0],0.25*mp::lattice_constant[1],0.25*mp::lattice_constant[2]};

		// For MPI version, only add local atoms
		#ifdef MPICF
			int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
		#else
			int num_local_atoms = atoms::num_atoms;
		#endif

		// Assign atoms to cells                                                                                                                  
		for(int atom=0;atom<num_local_atoms;atom++){
			double c[3]={atoms::x_coord_array[atom]+atom_offset[0],atoms::y_coord_array[atom]+atom_offset[1],atoms::z_coord_array[atom]+atom_offset[2]};
			int scc[3]={0,0,0}; // super cell coordinates
			for(int i=0;i<3;i++){
				scc[i]=int(c[i]/(mp::lattice_constant[i]*double(cells::size))); // Always round down for supercell coordinates
				// Always check cell in range
				if(scc[i]<0 || scc[i]>= d[i]){
					std::cerr << "Error - atom out of supercell range in neighbourlist calculation!" << std::endl;
					#ifdef MPICF
					std::cerr << "\tCPU Rank: " << vmpi::my_rank << std::endl;
					#endif 
					std::cerr << "\tAtom number:      " << atom << std::endl;
					std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
					std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
					std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
					#ifdef MPICF
						MPI::COMM_WORLD.Abort(EXIT_FAILURE);
						exit(EXIT_FAILURE);
					#else
						exit(EXIT_FAILURE);
					#endif
				}
			}
			// If no error for range then assign atom to cell.
			atoms::cell_array[atom]=supercell_array[scc[0]][scc[1]][scc[2]];
		}
		
		// Deallocate supercell array
		try{
			for(int i=0; i<d[0] ; i++){
				for(int j=0; j<d[1] ;j++){
					delete [] supercell_array[i][j];
				}
				delete [] supercell_array[i];
			}
		delete [] supercell_array;
		supercell_array=NULL;
		}
		catch(...){std::cout << "error deallocating supercell_array" << std::endl; exit(1);}
		
		// Resize new cell arrays
		cells::x_coord_array.resize(cells::num_cells,0.0);
		cells::y_coord_array.resize(cells::num_cells,0.0);
		cells::z_coord_array.resize(cells::num_cells,0.0);
		
		cells::x_mag_array.resize(cells::num_cells,0.0);
		cells::y_mag_array.resize(cells::num_cells,0.0);
		cells::z_mag_array.resize(cells::num_cells,0.0);
		
		cells::x_field_array.resize(cells::num_cells,0.0);
		cells::y_field_array.resize(cells::num_cells,0.0);
		cells::z_field_array.resize(cells::num_cells,0.0);
		
		cells::num_atoms_in_cell.resize(cells::num_cells,0);
		
		// Now add atoms to each cell
		for(int atom=0;atom<num_local_atoms;atom++){
			int local_cell=atoms::cell_array[atom];
			cells::x_coord_array[local_cell]+=atoms::x_coord_array[atom];
			cells::y_coord_array[local_cell]+=atoms::y_coord_array[atom];
			cells::z_coord_array[local_cell]+=atoms::z_coord_array[atom];
			cells::num_atoms_in_cell[local_cell]++;
		}

		// For MPI sum coordinates from all CPUs
		#ifdef MPICF
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::num_atoms_in_cell[0],cells::num_cells,MPI_INT,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::x_coord_array[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::y_coord_array[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::z_coord_array[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
		#endif
		
		vinfo << "=========================================================================" << std::endl;
		vinfo << "Number of atoms/cell: cell number, num atoms, coord" << std::endl;
		vinfo << "=========================================================================" << std::endl;

		// Now find mean coordinates
		for(int local_cell=0;local_cell<cells::num_cells;local_cell++){
			cells::x_coord_array[local_cell]/=double(cells::num_atoms_in_cell[local_cell]);
			cells::y_coord_array[local_cell]/=double(cells::num_atoms_in_cell[local_cell]);
			cells::z_coord_array[local_cell]/=double(cells::num_atoms_in_cell[local_cell]);
			vinfo << local_cell << "\t" << cells::num_atoms_in_cell[local_cell] << "\t";
			vinfo << cells::x_coord_array[local_cell] << "\t" << cells::y_coord_array[local_cell];
			vinfo << "\t" << cells::z_coord_array[local_cell] << "\t" << std::endl;
		}
		
		// Now re-update num_atoms in cell for local atoms only
		for(int atom=0;atom<num_local_atoms;atom++){
			int local_cell=atoms::cell_array[atom];
			cells::num_atoms_in_cell[local_cell]++;
		}
		return EXIT_SUCCESS;
	};

	
	
	
	int mag();
	int output_mag(std::ofstream&);
	
} // End of namespace cells
