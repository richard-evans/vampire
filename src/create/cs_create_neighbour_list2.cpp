///
/// @file 
/// @brief Neighbourlist generation routine
///
/// Generate list of neighbours for each atom based on unit cell 
/// interaction template
///
/// @author Richard Evans, richard.evans@york.ac.uk
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans 2009-2011. All Rights Reserved.
///
///  @internal
///    Created  08/06/2009
///   Revision  3.0
///  Copyright  Copyright (c) 2011, Richard Evans
///
///=====================================================================================
///

// Vampire Header files
#include "create.hpp"
#include "errors.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

// Standard Libraries
#include <cmath>
#include <iostream>
#include <limits>

namespace cs{

///
/// @brief Generate atomic neighbourlist
///
/// Assigns atoms to unit cells and then calculates all interactions between cells
///
/// Partial cells can exist so ensure enough cells are generated
///
///    4    5    6    7    8
///    | ...|....|....|.   | 
///
///  In this example offset=4, and max_cell = 8. Therefore 4 cells are needed.
///
///
int create_neighbourlist(std::vector<cs::catom_t> & catom_array, std::vector<std::vector <int> > & cneighbourlist){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::create_neighbourlist has been called" << std::endl;}	
	
	// put number of atoms into temporary variable
	const int num_atoms = catom_array.size();

	// Reserve space for num_atoms
	cneighbourlist.reserve(num_atoms);

	// Reserve space for each atom in neighbour list according to material type
	for(int atom=0;atom<num_atoms;atom++){
		int max_nn = 4;
		cneighbourlist.push_back(std::vector<int>());
		cneighbourlist[atom].reserve(max_nn);
	}

	// Calculate system dimensions and number of supercells
	unsigned int max_val=std::numeric_limits<unsigned int>::max();
	unsigned int min[3]={max_val,max_val,max_val}; // lowest cell id
	unsigned int max[3]={0,0,0}; // highest cell id

	for(int atom=0;atom<num_atoms;atom++){
		unsigned int c[3]={catom_array[atom].scx,catom_array[atom].scy,catom_array[atom].scz};
		for(int i=0;i<3;i++){
			if(c[i]<min[i]){
				min[i]=c[i];
			}
			if(c[i]>max[i]){
				max[i]=c[i];
			}
		}
	}
	
	// calculate offset and cell maximum in whole unit cells
	const unsigned int offset[3] = {min[0], min[1], min[2]};
	const unsigned int max_cell[3] = {max[0],max[1],max[2]};
	
	// calculate number of cells needed = max-min+1 ( if max_cell = 25, then 0-25 = 26
	const unsigned int d[3]={max_cell[0]-offset[0]+1,max_cell[1]-offset[1]+1,max_cell[2]-offset[2]+1};

	// Declare array for create space for 3D supercell array
	std::vector<std::vector<std::vector<std::vector<int> > > > supercell_array;

	std::cout << "Memory required for neighbourlist calculation:" << 8.0*double(d[0])*double(d[1])*double(d[2])*double(unit_cell.atom.size())/1.0e6 << " MB" << std::endl;
	supercell_array.resize(d[0]);
	for(unsigned int i=0; i<d[0] ; i++){
		supercell_array[i].resize(d[1]);
		for(unsigned int j=0; j<d[1] ; j++){
			supercell_array[i][j].resize(d[2]);
			for(unsigned int k=0; k<d[2] ; k++){
				supercell_array[i][j][k].resize(unit_cell.atom.size(),-1);
			}
		}
	}

	// declare cell array to loop over
	const unsigned int num_cells=d[0]*d[1]*d[2];
	std::vector< std::vector <int> > cell_coord_array;
	cell_coord_array.reserve(num_cells);
	for(unsigned int i=0;i<num_cells;i++){
		cell_coord_array.push_back(std::vector<int>());
		cell_coord_array[i].resize(3);
	}

	// Initialise cell_array
	unsigned int cell=0;
	for(unsigned int x=0;x<d[0];x++){
		for(unsigned int y=0;y<d[1];y++){
			for(unsigned int z=0;z<d[2];z++){
				//save cell coordinates
				cell_coord_array[cell][0]=x;
				cell_coord_array[cell][1]=y;
				cell_coord_array[cell][2]=z;
				cell++;
			}
		}
	}

	// Populate supercell array with atom numbers
	for(int atom=0;atom<num_atoms;atom++){
		unsigned int scc[3]={catom_array[atom].scx-offset[0],catom_array[atom].scy-offset[1],catom_array[atom].scz-offset[2]};
		
		double c[3]={catom_array[atom].x,catom_array[atom].y,catom_array[atom].z};
		//std::cout << atom << "\t" << c[0] << "\t" << c[1] <<"\t" << c[2] << std::endl;
		for(int i=0;i<3;i++){
			scc[i]=int(c[i]/cs::unit_cell_size[i])-offset[i]; // Always round down for supercell coordinates
			// Always check cell in range
			if(scc[i]<0 || scc[i]>= d[i]){
				//std::cerr << "Error - atom out of supercell range in neighbourlist calculation!" << std::endl;
				#ifdef MPICF
				std::cerr << "\tCPU Rank: " << vmpi::my_rank << std::endl;
				#endif 
				std::cerr << "\tAtom number:      " << atom << std::endl;
				std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
				std::cerr << "\tmin coordinates:  " << min[0] << "\t" << min[1] << "\t" << min[2] << "\t" << std::endl;
				std::cerr << "\tmax coordinates:  " << max[0] << "\t" << max[1] << "\t" << max[2] << "\t" << std::endl;
				std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
				std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
				std::cerr << "\tCell offset:      " << offset[0] << "\t" << offset[1] << "\t" << offset[2] << std::endl;
				std::cerr << "\tCell offest (dp): " << offset[0]*cs::unit_cell_size[0] << "\t" << offset[1]*cs::unit_cell_size[1] << "\t" << offset[2]*cs::unit_cell_size[2] << std::endl;
				err::vexit();
			}
		}

		// Check for atoms greater than max_atoms_per_supercell
		if(catom_array[atom].uc_category<unit_cell.atom.size()){
			// Add atom to supercell
			supercell_array[scc[0]][scc[1]][scc[2]][catom_array[atom].uc_category]=atom;
		}
		else{
			std::cerr << "Error, number of atoms per supercell exceeded" << std::endl;
			std::cerr << "\tAtom number:      " << atom << std::endl;
			std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
			std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
			std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
			std::cerr << "\tCell offset:      " << offset[0] << "\t" << offset[1] << "\t" << offset[2] << std::endl;
			std::cerr << "\tAtoms in Current Cell:" << std::endl;
			for(unsigned int ix=0;ix<supercell_array[scc[0]][scc[1]][scc[2]].size();ix++){
				const int ixatom=supercell_array[scc[0]][scc[1]][scc[2]][ix];
				std::cerr << "\t\t [id x y z] "<< ix << "\t" << ixatom << "\t" << catom_array[ixatom].x << "\t" << catom_array[ixatom].y << "\t" << catom_array[ixatom].z << std::endl;
			}
			err::vexit();
		}
	}

	// Generate neighbour list
	std::cout <<"Generating Neighbour list"<< std::flush; 

	// Loop over all cells
	for(unsigned int cell=0;cell<num_cells;cell++){
		if(cell%(num_cells/10+1)==0){
			std::cout << "." << std::flush;
		}
		//std::cout << "cell: " << cell << ":>" << std::flush; // << std::endl;
		int scc[3]={cell_coord_array[cell][0],cell_coord_array[cell][1],cell_coord_array[cell][2]};
		// Loop over all interactions
		for(unsigned int i=0;i<cs::unit_cell.interaction.size();i++){
			const int atom=cs::unit_cell.interaction[i].i;
			const int natom=cs::unit_cell.interaction[i].j;
			const int nx=cs::unit_cell.interaction[i].dx+scc[0];
			const int ny=cs::unit_cell.interaction[i].dy+scc[1];
			const int nz=cs::unit_cell.interaction[i].dz+scc[2];
			// check for out-of-bounds access
			if((nx>=0 && nx<d[0]) && (ny>=0 && ny<d[1]) && (nz>=0 && nz<d[2])){
				if((supercell_array[nx][ny][nz][atom]!=-1) || (supercell_array[nx][ny][nz][natom]!=-1)){
					cneighbourlist[supercell_array[scc[0]][scc[1]][scc[2]][atom]].push_back(supercell_array[nx][ny][nz][natom]);
					//std::cout << supercell_array[nx][ny][nz][atom] << "\t" << supercell_array[nx][ny][nz][natom] << std::endl;
					//std::cin.get();
				}
			}
		}
	}
	
	std::cout << "done!" << std::endl;

	// Deallocate supercell array
	for(unsigned int i=0; i<d[0] ; i++){
		for(unsigned int j=0; j<d[1] ;j++){
			for(unsigned int k=0; k<d[2] ;k++){
				supercell_array[i][j][k].resize(0);
				}
				supercell_array[i][j].resize(0);
			}
			supercell_array[i].resize(0);
		}
	supercell_array.resize(0);

	// Print neighbour list
	//for(int atom=0;atom<catom_array.size();atom++){
	//	std::cout << atom << "\t";
	//	for(int nn=0;nn<cneighbourlist[atom].size();nn++){
	//		std::cout << cneighbourlist[atom][nn] << "\t";
	//	}
	//	std::cout << std::endl;
	//}
		
	// Mark surface atoms
	//for(int atom=0;atom<num_atoms;atom++){
	//	if(int(cneighbourlist[atom].size())!=mp::material[catom_array[atom].material].hamiltonian_num_neighbours){
	//		catom_array[atom].material=2;
	//	}
	//}

	return EXIT_SUCCESS;
}

} // End of namespace cs
