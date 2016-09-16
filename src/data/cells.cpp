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
///
/// @file
/// @brief Contains cells namespace and asssociated functions
///
/// @details Cells discretise the system into small blocks suitable for outputting
///          magnetisation data for large systems, and also for calculating demag fields
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    28/03/2011
/// @internal
///	Created:		11/11/2010
///	Revision:	  ---
///=====================================================================================
///

// Vampire Header files
#include "atoms.hpp"
#include "cells.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

// System header files
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace cells{

	int num_cells=0;
	int num_local_cells=0;
   int num_atoms_in_unit_cell=0;
	double size=7.0; // Angstroms

	bool initialised=false;

	std::vector <int> num_atoms_in_cell;
	std::vector <int> local_cell_array;

	std::vector <double> x_coord_array;
	std::vector <double> y_coord_array;
	std::vector <double> z_coord_array;

	std::vector <double> x_mag_array;
	std::vector <double> y_mag_array;
	std::vector <double> z_mag_array;

	std::vector <double> x_field_array;
	std::vector <double> y_field_array;
	std::vector <double> z_field_array;

   std::vector <double> volume_array;

/// @brief Cell initialiser function
///
/// @details Determines number of cells, assigns atoms to cells, and sets up cell arrays
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    28/03/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///=====================================================================================
///
	int initialise(){

		// check calling of routine if error checking is activated
		if(err::check==true) std::cout << "cells::initialise has been called" << std::endl;

		cells::num_cells=0;
		cells::num_local_cells=0;

		zlog << zTs() << "Macrocell size = " << cells::size << " Angstroms" << std::endl;

		// determine number of cells in each direction (with small shift to prevent the fence post problem)
		unsigned int ncellx = static_cast<unsigned int>(ceil((cs::system_dimensions[0]+0.01)/cells::size));
		unsigned int ncelly = static_cast<unsigned int>(ceil((cs::system_dimensions[1]+0.01)/cells::size));
		unsigned int ncellz = static_cast<unsigned int>(ceil((cs::system_dimensions[2]+0.01)/cells::size));

		//update total number of cells
		cells::num_cells=ncellx*ncelly*ncellz;

		zlog << zTs() << "Macrocells in x,y,z: " << ncellx << "\t" << ncelly << "\t" << ncellz << std::endl;
		zlog << zTs() << "Total number of macrocells: " << cells::num_cells << std::endl;
		zlog << zTs() << "Memory required for macrocell arrays: " << 80.0*double(cells::num_cells)/1.0e6 << " MB" << std::endl;

		// Determine number of cells in x,y,z
		const unsigned int d[3]={ncellx,ncelly,ncellz};

		// Set cell counter
		int cell=0;

		// Declare array for create space for 3D supercell array
		int*** supercell_array;
		//std::cout << "Memory required for cell list calculation:" << 8.0*double(d[0])*double(d[1])*double(d[2])/1.0e6 << " MB" << std::endl;
		try{supercell_array=new int**[d[0]];
			for(unsigned int i=0; i<d[0] ; i++){
				supercell_array[i]=new int*[d[1]];
				for(unsigned int j=0; j<d[1] ; j++){
					supercell_array[i][j]=new int[d[2]];
					for(unsigned int k=0; k<d[2] ; k++){
						supercell_array[i][j][k]=cell;
						cell++;
					}
				}
			}
		}
		catch(...){
			terminaltextcolor(RED);
			std::cerr << "Error allocating supercell_array for macrocell list calculation" << std::endl;
			terminaltextcolor(WHITE);
         zlog << zTs() << "Error allocating supercell_array for macrocell list calculation" << std::endl;
         err::vexit();
      }

		// slightly offset atomic coordinates to prevent fence post problem
      double atom_offset[3]={0.01,0.01,0.01};

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
				scc[i]=int(c[i]/cells::size); // Always round down for supercell coordinates
				// Always check cell in range
				if(scc[i]<0 || static_cast<unsigned int>(scc[i]) >= d[i]){
					terminaltextcolor(RED);
					std::cerr << "Error - atom out of supercell range in neighbourlist calculation!" << std::endl;
					terminaltextcolor(WHITE);
               zlog << zTs() << "Error - atom out of supercell range in neighbourlist calculation!" << std::endl;
					#ifdef MPICF
					terminaltextcolor(RED);
					std::cerr << "\tCPU Rank: " << vmpi::my_rank << std::endl;
					terminaltextcolor(WHITE);
					#endif
					terminaltextcolor(RED);
					std::cerr << "\tAtom number:      " << atom << std::endl;
					std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
					std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
					std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
					terminaltextcolor(WHITE);
					err::vexit();
				}
			}
			// If no error for range then assign atom to cell.
			atoms::cell_array[atom]=supercell_array[scc[0]][scc[1]][scc[2]];
		}

		// Deallocate supercell array
		try{
			for(unsigned int i=0; i<d[0] ; i++){
				for(unsigned int j=0; j<d[1] ;j++){
					delete [] supercell_array[i][j];
				}
				delete [] supercell_array[i];
			}
		delete [] supercell_array;
		supercell_array=NULL;
		}
		catch(...){
         zlog << zTs() << "error deallocating supercell_array" << std::endl;
         err::vexit();
      }

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
      cells::volume_array.resize(cells::num_cells,0.0);

      std::vector<double> total_moment_array(cells::num_cells,0.0);

		// Now add atoms to each cell as magnetic 'centre of mass'
		for(int atom=0;atom<num_local_atoms;atom++){
			int local_cell=atoms::cell_array[atom];
         int type = atoms::type_array[atom];
         const double mus = mp::material[type].mu_s_SI;
			cells::x_coord_array[local_cell]+=atoms::x_coord_array[atom]*mus;
			cells::y_coord_array[local_cell]+=atoms::y_coord_array[atom]*mus;
			cells::z_coord_array[local_cell]+=atoms::z_coord_array[atom]*mus;
         total_moment_array[local_cell]+=mus;
         cells::num_atoms_in_cell[local_cell]++;
		}

		// For MPI sum coordinates from all CPUs
		#ifdef MPICF
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::num_atoms_in_cell[0],cells::num_cells,MPI_INT,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::x_coord_array[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::y_coord_array[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::z_coord_array[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
         MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&total_moment_array[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
      #endif

		//if(vmpi::my_rank==0){
			//vinfo << "=========================================================================" << std::endl;
			//vinfo << "Number of atoms/cell: cell number, num atoms, coord" << std::endl;
			//vinfo << "=========================================================================" << std::endl;
		//}

      // Used to calculate magnetisation in each cell. Poor approximation when unit cell size ~ system size.
      const double atomic_volume = cs::unit_cell.dimensions[0]*cs::unit_cell.dimensions[1]*cs::unit_cell.dimensions[2]/cells::num_atoms_in_unit_cell;

		// Now find mean coordinates via magnetic 'centre of mass'
		for(int local_cell=0;local_cell<cells::num_cells;local_cell++){
			if(cells::num_atoms_in_cell[local_cell]>0){
            cells::x_coord_array[local_cell] = cells::x_coord_array[local_cell]/(total_moment_array[local_cell]);
            cells::y_coord_array[local_cell] = cells::y_coord_array[local_cell]/(total_moment_array[local_cell]);
            cells::z_coord_array[local_cell] = cells::z_coord_array[local_cell]/(total_moment_array[local_cell]);
            cells::volume_array[local_cell] = double(cells::num_atoms_in_cell[local_cell])*atomic_volume;
         }
			//if(vmpi::my_rank==0){
			//vinfo << local_cell << "\t" << cells::num_atoms_in_cell[local_cell] << "\t";
			//vinfo << cells::x_coord_array[local_cell] << "\t" << cells::y_coord_array[local_cell];
			//vinfo << "\t" << cells::z_coord_array[local_cell] << "\t" << std::endl;
			//}
		}

		//Set number of atoms in cell to zero
		for(int cell=0;cell<cells::num_cells;cell++){
		  cells::num_atoms_in_cell[cell]=0;
		}

		// Now re-update num_atoms in cell for local atoms only
		for(int atom=0;atom<num_local_atoms;atom++){
			int local_cell=atoms::cell_array[atom];
			cells::num_atoms_in_cell[local_cell]++;
		}

		// Calculate number of local cells
		for(int cell=0;cell<cells::num_cells;cell++){
			if(cells::num_atoms_in_cell[cell]!=0){
				cells::local_cell_array.push_back(cell);
				cells::num_local_cells++;
			}
		}

		zlog << zTs() << "Number of local macrocells on rank " << vmpi::my_rank << ": " << cells::num_local_cells << std::endl;

		cells::initialised=true;

      // Precalculate cell magnetisation
      cells::mag();

		return EXIT_SUCCESS;
	}

/// @brief Cell magnetisation function
///
/// @details Determines moment in each cell (J/T)
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    28/03/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///=====================================================================================
///
int mag() {

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "cells::mag has been called" << std::endl;

  // Check for initialised arrays
  //if(cells::initialised!=true) cells::initialise();

  for(int i=0; i<cells::num_cells; ++i) {
    cells::x_mag_array[i] = 0.0;
    cells::y_mag_array[i] = 0.0;
    cells::z_mag_array[i] = 0.0;
  }

#ifdef MPICF
    int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
#else
    int num_local_atoms = atoms::num_atoms;
#endif

  // calulate total moment in each cell
  for(int i=0;i<num_local_atoms;++i) {
    int cell = atoms::cell_array[i];
    int type = atoms::type_array[i];
    const double mus = mp::material[type].mu_s_SI;

    cells::x_mag_array[cell] += atoms::x_spin_array[i]*mus;
    cells::y_mag_array[cell] += atoms::y_spin_array[i]*mus;
    cells::z_mag_array[cell] += atoms::z_spin_array[i]*mus;
  }

#ifdef MPICF
  // Reduce magnetisation on all nodes
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::x_mag_array[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::y_mag_array[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,&cells::z_mag_array[0],cells::num_cells,MPI_DOUBLE,MPI_SUM);
#endif

  return EXIT_SUCCESS;
}

}  // End of namespace cells
