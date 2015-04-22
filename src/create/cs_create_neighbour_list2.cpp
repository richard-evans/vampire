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
#include "vio.hpp"
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
int create_neighbourlist(std::vector<cs::catom_t> & catom_array, std::vector<std::vector <neighbour_t> > & cneighbourlist){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::create_neighbourlist has been called" << std::endl;}	
	
	// put number of atoms into temporary variable
	const int num_atoms = catom_array.size();

	// Reserve space for num_atoms
	cneighbourlist.reserve(num_atoms);

	// estimate number of interactions per atom
	const int max_nn=int(1.1*(double(unit_cell.interaction.size())/double(unit_cell.atom.size())));
	
	// Reserve space for each atom in neighbour list according to material type
	for(int atom=0;atom<num_atoms;atom++){
		cneighbourlist.push_back(std::vector<neighbour_t>());
		cneighbourlist[atom].reserve(max_nn);
	}

   // Calculate system dimensions and number of supercells
   int max_val=std::numeric_limits<int>::max();
   int min[3]={max_val,max_val,max_val}; // lowest cell id
   int max[3]={0,0,0}; // highest cell id

	for(int atom=0;atom<num_atoms;atom++){
		int c[3]={catom_array[atom].scx,catom_array[atom].scy,catom_array[atom].scz};
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
	const int offset[3] = {min[0], min[1], min[2]};
	const int max_cell[3] = {max[0],max[1],max[2]};
	
	// calculate number of cells needed = max-min+1 ( if max_cell = 25, then 0-25 = 26
	const unsigned int d[3]={max_cell[0]-offset[0]+1,max_cell[1]-offset[1]+1,max_cell[2]-offset[2]+1};

	// Declare array for create space for 3D supercell array
	std::vector<std::vector<std::vector<std::vector<int> > > > supercell_array;

	zlog << zTs() << "Memory required for neighbourlist calculation:" << 8.0*double(d[0])*double(d[1])*double(d[2])*double(unit_cell.atom.size())/1.0e6 << " MB" << std::endl;
   zlog << zTs() << "Allocating memory for supercell array in neighbourlist calculation..."<< std::endl;
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
   zlog << zTs() << "\tDone"<< std::endl;

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
   zlog << zTs() << "Populating supercell array for neighbourlist calculation..."<< std::endl;
	// Populate supercell array with atom numbers
	for(int atom=0;atom<num_atoms;atom++){
		unsigned int scc[3]={catom_array[atom].scx-offset[0],catom_array[atom].scy-offset[1],catom_array[atom].scz-offset[2]};
		
		double c[3]={catom_array[atom].x,catom_array[atom].y,catom_array[atom].z};
		//std::cout << atom << "\t" << c[0] << "\t" << c[1] <<"\t" << c[2] << std::endl;
		for(int i=0;i<3;i++){
			//scc[i]=int(c[i]/cs::unit_cell_size[i])-offset[i]; // Always round down for supercell coordinates
			// Always check cell in range
         if(scc[i]>= d[i]){
			//if(scc[i]<0 || scc[i]>= d[i]){ // Chexk for scc < 0 not required since d and scc are unsigned
				//std::cerr << "Error - atom out of supercell range in neighbourlist calculation!" << std::endl;
				#ifdef MPICF
				terminaltextcolor(RED);
				std::cerr << "\tCPU Rank: " << vmpi::my_rank << std::endl;
				terminaltextcolor(WHITE);
				#endif 
				terminaltextcolor(RED);
				std::cerr << "\tAtom number:      " << atom << std::endl;
				std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
				std::cerr << "\tmin coordinates:  " << min[0] << "\t" << min[1] << "\t" << min[2] << "\t" << std::endl;
				std::cerr << "\tmax coordinates:  " << max[0] << "\t" << max[1] << "\t" << max[2] << "\t" << std::endl;
				std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
				std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
				std::cerr << "\tCell offset:      " << offset[0] << "\t" << offset[1] << "\t" << offset[2] << std::endl;
				std::cerr << "\tCell offest (dp): " << offset[0]*cs::unit_cell_size[0] << "\t" << offset[1]*cs::unit_cell_size[1] << "\t" << offset[2]*cs::unit_cell_size[2] << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		// Check for atoms greater than max_atoms_per_supercell
		if(catom_array[atom].uc_id<unit_cell.atom.size()){
			// Add atom to supercell
			supercell_array[scc[0]][scc[1]][scc[2]][catom_array[atom].uc_id]=atom;
		}
		else{
			terminaltextcolor(RED);
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
			terminaltextcolor(WHITE);
			err::vexit();
		}
	}

   zlog << zTs() << "\tDone"<< std::endl;

   // Get unit cell size
   const double ucdx=cs::unit_cell.dimensions[0];
   const double ucdy=cs::unit_cell.dimensions[1];
   const double ucdz=cs::unit_cell.dimensions[2];

	// Generate neighbour list
	std::cout <<"Generating neighbour list"<< std::flush;
   zlog << zTs() << "Memory required for neighbour list:" << 8.0*double(num_cells)*double(cs::unit_cell.interaction.size())/1.0e6 << " MB" << std::endl;
   zlog << zTs() << "Generating neighbour list..."<< std::endl;
	neighbour_t tmp_nt;
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

			int nx=cs::unit_cell.interaction[i].dx+scc[0];
			int ny=cs::unit_cell.interaction[i].dy+scc[1];
			int nz=cs::unit_cell.interaction[i].dz+scc[2];

         // vector from i->j
         double vx=0.0;
         double vy=0.0;
         double vz=0.0;

         #ifdef MPICF
           // Parallel periodic boundaries are handled explicitly elsewhere
         #else
         // Wrap around for periodic boundaries
         // Consider virtual atom position for position vector
         if(cs::pbc[0]==true){
            if(nx>=int(d[0])){
               nx=nx-d[0];
               vx=vx+d[0]*ucdx;
            }
            else if(nx<0){
               nx=nx+d[0];
               vx=vx-d[0]*ucdx;
            }
         }
         if(cs::pbc[1]==true){
            if(ny>=int(d[1])){
               ny=ny-d[1];
               vy=vy+d[1]*ucdy;
            }
            else if(ny<0){
               ny=ny+d[1];
               vy=vy-d[1]*ucdy;
            }
         }
         if(cs::pbc[2]==true){
            if(nz>=int(d[2])){
               nz=nz-d[2];
               vz=vz+d[2]*ucdz;
            }
            else if(nz<0){
               nz=nz+d[2];
               vz=vz-d[2]*ucdz;
            }
         }
         #endif
         // check for out-of-bounds access
         if((nx>=0 && nx<d[0]) && (ny>=0 && ny<d[1]) && (nz>=0 && nz<d[2])){
            // check for missing atoms
            if((supercell_array[scc[0]][scc[1]][scc[2]][atom]!=-1) && (supercell_array[nx][ny][nz][natom]!=-1)){

               // need actual atom numbers...
               int atomi = supercell_array[scc[0]][scc[1]][scc[2]][atom];
               int atomj = supercell_array[nx][ny][nz][natom];

               //std::cout << "int_id: " << i << "\tatom i: " << atomi << "\tatom j: " << atomj << "\tuc_i: " << atom << "\tuc_j: " << natom << std::endl;  

               double ix=catom_array[atomi].x; // Already in A
               double iy=catom_array[atomi].y;
               double iz=catom_array[atomi].z;
               double jx=catom_array[atomj].x;
               double jy=catom_array[atomj].y;
               double jz=catom_array[atomj].z;

               //std::cout << "\tpi:    " << ix << "\t" << iy << "\t" << iz << std::endl;
               //std::cout << "\tpj:    " << jx << "\t" << jy << "\t" << jz << std::endl;
               //std::cout << "\tv_uc:  " << vx << "\t" << vy << "\t" << vz << std::endl;

               vx+=jx-ix;
               vy+=jy-iy;
               vz+=jz-iz;

               //std::cout << "\tv:     " << jx-ix << "\t" << jy-iy << "\t" << jz-iz << std::endl;
               //std::cout << "\tv_eff: " << vx << "\t" << vy << "\t" << vz << std::endl;

               // get current index
               int index=cneighbourlist[supercell_array[scc[0]][scc[1]][scc[2]][atom]].size();

               // push back array of class
               cneighbourlist[supercell_array[scc[0]][scc[1]][scc[2]][atom]].push_back(tmp_nt);

               // now save atom id and interaction type
               cneighbourlist[supercell_array[scc[0]][scc[1]][scc[2]][atom]][index].nn=supercell_array[nx][ny][nz][natom];
               cneighbourlist[supercell_array[scc[0]][scc[1]][scc[2]][atom]][index].i=i;

               // Add position vector from i-> j
               cneighbourlist[supercell_array[scc[0]][scc[1]][scc[2]][atom]][index].vx=vx;
               cneighbourlist[supercell_array[scc[0]][scc[1]][scc[2]][atom]][index].vy=vy;
               cneighbourlist[supercell_array[scc[0]][scc[1]][scc[2]][atom]][index].vz=vz;

            }
			}
		}
	}
	terminaltextcolor(GREEN);
	std::cout << "done!" << std::endl;
	terminaltextcolor(WHITE);
   zlog << zTs() << "\tDone"<< std::endl;

	// Deallocate supercell array
   zlog << zTs() << "Deallocating supercell array for neighbour list calculation" << std::endl;
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

   zlog << zTs() << "\tDone" << std::endl;

	// Print neighbour list
	//for(int atom=0;atom<catom_array.size();atom++){
	//	std::cout << atom << "\t";
	//	for(int nn=0;nn<cneighbourlist[atom].size();nn++){
	//		std::cout << cneighbourlist[atom][nn].nn << "\t";
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
