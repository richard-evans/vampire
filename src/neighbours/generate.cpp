//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <iostream>

// Vampire headers
#include "create_atoms_class.hpp" // class definition for atoms in create module
#include "errors.hpp"
#include "neighbours.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

//-----------------------------------
// Fix for horrible windows compiler
//-----------------------------------
/*#ifdef WIN_COMPILE
   #define NOMINMAX
   #undef max
   #undef min
#endif*/

namespace neighbours{

//----------------------------------------------------------------------------------
// @brief Generate atomic neighbourlist for a generalised exchange template
//
// Assigns atoms to unit cells and then calculates all interactions between cells
//
// Partial cells can exist so ensure enough cells are generated
//
//    4    5    6    7    8
//    | ...|....|....|.   |
//
//  In this example offset=4, and max_cell = 8. Therefore 4 cells are needed.
//
//----------------------------------------------------------------------------------
void list_t::generate( std::vector<cs::catom_t>& atom_array,    // array of atoms (as reference for speed)
               unitcell::exchange_template_t& exchange, // exchange template to calculate neighbour list
               const unsigned int num_atoms_in_unit_cell,        // number of atoms in each cell to estimate interaction numbers
               const double ucdx,                       // unit cell size
               const double ucdy,
               const double ucdz
             ){

	// put number of atoms into temporary variable
	const int num_atoms = atom_array.size();

	// Reserve space for num_atoms
	list.reserve(num_atoms);

	// estimate number of interactions per atom
	const int64_t max_nn = int64_t( 1.1*( double(exchange.interaction.size()) / double(num_atoms_in_unit_cell) ) );

	// Reserve space for each atom in neighbour list according to material type
	for(int atom=0; atom < num_atoms; atom++){
		list.push_back(std::vector<neighbour_t>());
		list[atom].reserve(max_nn);
	}

   // Calculate system dimensions and number of supercells
   const int64_t max_val=1000000000000;
   int64_t min[3] = {max_val,max_val,max_val}; // lowest cell id
   int64_t max[3] = {0,0,0}; // highest cell id

   // find supercell range of atoms on this CPU
	for(int atom = 0; atom < num_atoms; atom++){

		int64_t c[3] = { atom_array[atom].scx,
                       atom_array[atom].scy,
                       atom_array[atom].scz};

      // loop over i,j,k
		for(int i = 0; i < 3; i++){
			if( c[i] < min[i] ){
				min[i] = c[i];
			}
			if( c[i] > max[i] ){
				max[i] = c[i];
			}
		}

	}

   // check for out of range value
   // loop over i,j,k
   for(int i = 0; i < 3; i++){
      if(min[i] > max_val){
         std::cerr << "Programmer error! too many supercells in atom list" << std::endl;
      }
   }

	// calculate offset and cell maximum in whole unit cells
	const int64_t offset[3]   = {min[0], min[1], min[2]};
	const int64_t max_cell[3] = {max[0], max[1], max[2]};

	// calculate number of cells needed = max-min+1
   // ( if max_cell = 25, then 0 to 25 = 26 cells)
	const int64_t d[3] = { ( max_cell[0] - offset[0] + 1 ),
                          ( max_cell[1] - offset[1] + 1 ),
                          ( max_cell[2] - offset[2] + 1 )};

	// Declare temporary array for 3D supercell array
	std::vector<std::vector<std::vector<std::vector<int> > > > supercell_array;

   // calculate total number of neighbours and inform user of memory needed
   double num_neighbours = double(d[0]) * double(d[1]) * double(d[2]) * double(num_atoms_in_unit_cell);
   #ifdef MPICF
      // calculate total interactions for entire system
      double total_neighbours = 0.0;
      MPI_Allreduce(&total_neighbours, &num_neighbours, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if(vmpi::master){
         zlog << zTs() << "Memory required for neighbourlist calculation (each cpu):" <<
         8.0*total_neighbours/(vmpi::num_processors * 1.0e6) << " MB" << std::endl;
         zlog << zTs() << "Memory required for neighbourlist calculation (all cpus):" <<
         8.0*total_neighbours/1.0e6 << " MB" << std::endl;
      }
   #else
      zlog << zTs() << "Memory required for neighbourlist calculation:" <<
      8.0*num_neighbours/1.0e6 << " MB" << std::endl;
   #endif

   // Inform user that neighbour list calculation is beginning
   zlog << zTs() << "Allocating memory for supercell array in neighbourlist calculation" << std::endl;

   // Allocate supercell array to list all atoms
	supercell_array.resize(d[0]);
	for(unsigned int i=0; i < d[0] ; i++){
		supercell_array[i].resize(d[1]);
		for(unsigned int j=0; j<d[1] ; j++){
			supercell_array[i][j].resize(d[2]);
			for(unsigned int k=0; k<d[2] ; k++){
				supercell_array[i][j][k].resize(num_atoms_in_unit_cell,-1);
			}
		}
	}
   zlog << zTs() << "\tAllocating memory done"<< std::endl;

	// declare 1D cell array to loop over
	const uint64_t num_cells = d[0]*d[1]*d[2];
	std::vector< std::vector <int> > cell_coord_array;
	cell_coord_array.reserve(num_cells);
	for(uint64_t i=0; i < num_cells; i++){
		cell_coord_array.push_back(std::vector<int>());
		cell_coord_array[i].resize(3);
	}

	// Initialise cell_array
	uint64_t cell=0; // cell counter
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

   // Inform user of time intensive process
   zlog << zTs() << "Populating supercell array for neighbourlist calculation..."<< std::endl;

	// Populate supercell array with atom numbers
	for(int atom=0; atom < num_atoms; atom++){

      // get supercell coordinates
      int64_t scc[3]={ atom_array[atom].scx - offset[0],
                       atom_array[atom].scy - offset[1],
                       atom_array[atom].scz - offset[2] };

      // get atom real position coordinates
		double c[3]={atom_array[atom].x, atom_array[atom].y, atom_array[atom].z};

      // check that atom is within valid range of supercell coordinates
		for(int i=0;i<3;i++){
			// Always check cell in range
         if( scc[i] >= d[i] ){
            terminaltextcolor(RED);
				std::cerr << "Error - atom out of supercell range in neighbourlist calculation!" << std::endl;
            terminaltextcolor(WHITE);
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
				std::cerr << "\tCell offest (dp): " << offset[0]*ucdx << "\t" << offset[1]*ucdy << "\t" << offset[2]*ucdz << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		// Check for atoms greater than max_atoms_per_supercell
		if(atom_array[atom].uc_id < num_atoms_in_unit_cell){
			// Add atom to supercell
			supercell_array[scc[0]][scc[1]][scc[2]][atom_array[atom].uc_id]=atom;
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
				std::cerr << "\t\t [id x y z] "<< ix << "\t" << ixatom << "\t" << atom_array[ixatom].x << "\t" << atom_array[ixatom].y << "\t" << atom_array[ixatom].z << std::endl;
			}
			terminaltextcolor(WHITE);
			err::vexit();
		}
	}

   // Inform user of progress
   zlog << zTs() << "\tPopulating supercell array completed"<< std::endl;

   // calculate total number of neighbours and inform user of memory needed
   num_neighbours = double(num_cells)*double(exchange.interaction.size());
   #ifdef MPICF
      // calculate total interactions for entire system
      total_neighbours = 0.0;
      MPI_Allreduce(&total_neighbours, &num_neighbours, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if(vmpi::master){
         zlog << zTs() << "Memory required for neighbour list (each cpu):" <<
         8.0*total_neighbours/(vmpi::num_processors * 1.0e6) << " MB" << std::endl;
         zlog << zTs() << "Memory required for neighbour list (all cpus):" <<
         8.0*total_neighbours/1.0e6 << " MB" << std::endl;
      }
   #else
      zlog << zTs() << "Memory required for neighbour list:" <<
      8.0*num_neighbours/1.0e6 << " MB" << std::endl;
   #endif

	// Generate neighbour list and inform user
	std::cout <<"Generating neighbour list"<< std::flush;
   zlog << zTs() << "Generating neighbour list..."<< std::endl;

   // temporary neighour to simplify memory management
	neighbour_t tmp_nt;

   // copy number of interactions to temporary constant
   const unsigned int num_interactions = exchange.interaction.size();

	// Loop over all cells
	for(uint64_t cell = 0; cell < num_cells; cell++){

      // Print out progress indicator to user
		if( cell % ( num_cells / 10 + 1 ) == 0 ){
			std::cout << "." << std::flush;
		}

      // get supercell coordinates of cell
		int scc[3]={ cell_coord_array[cell][0],
                   cell_coord_array[cell][1],
                   cell_coord_array[cell][2]};

		// Loop over all interactions in exchange template
		for(unsigned int i = 0; i < num_interactions; i++){

			const int atom=exchange.interaction[i].i;
			const int natom=exchange.interaction[i].j;

			int nx = exchange.interaction[i].dx + scc[0];
			int ny = exchange.interaction[i].dy + scc[1];
			int nz = exchange.interaction[i].dz + scc[2];

         // vector from i->j
         double vx=0.0;
         double vy=0.0;
         double vz=0.0;

         #ifdef MPICF
           // Parallel periodic boundaries are handled explicitly during the
           // halo region setup
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
         if( (nx >= 0 && static_cast<int64_t>(nx) < d[0] ) &&
             (ny >= 0 && static_cast<int64_t>(ny) < d[1] ) &&
             (nz >= 0 && static_cast<int64_t>(nz) < d[2] ) ){
            // check for missing atoms
            if((supercell_array[scc[0]][scc[1]][scc[2]][atom]!=-1) && (supercell_array[nx][ny][nz][natom]!=-1)){

               // need actual atom numbers...
               int atomi = supercell_array[scc[0]][scc[1]][scc[2]][atom];
               int atomj = supercell_array[nx][ny][nz][natom];

               //std::cout << "int_id: " << i << "\tatom i: " << atomi << "\tatom j: " << atomj << "\tuc_i: " << atom << "\tuc_j: " << natom << std::endl;

               // Load atom positions
               double ix = atom_array[atomi].x; // Already in A
               double iy = atom_array[atomi].y;
               double iz = atom_array[atomi].z;
               double jx = atom_array[atomj].x;
               double jy = atom_array[atomj].y;
               double jz = atom_array[atomj].z;

               //std::cout << "\tpi:    " << ix << "\t" << iy << "\t" << iz << std::endl;
               //std::cout << "\tpj:    " << jx << "\t" << jy << "\t" << jz << std::endl;
               //std::cout << "\tv_uc:  " << vx << "\t" << vy << "\t" << vz << std::endl;

               vx += jx-ix;
               vy += jy-iy;
               vz += jz-iz;

               //std::cout << "\tv:     " << jx-ix << "\t" << jy-iy << "\t" << jz-iz << std::endl;
               //std::cout << "\tv_eff: " << vx << "\t" << vy << "\t" << vz << std::endl;

               // get current index
               // int index = list[supercell_array[scc[0]][scc[1]][scc[2]][atom]].size(); // unused variable

               // set neighbour data
               tmp_nt.nn = supercell_array[nx][ny][nz][natom]; // atom ID of neighbour
               tmp_nt.i = i;                                   // interaction type
               tmp_nt.vx = vx;                                 // position vector i->j
               tmp_nt.vy = vy;
               tmp_nt.vz = vz;

               // push back array of class
               list[supercell_array[scc[0]][scc[1]][scc[2]][atom]].push_back(tmp_nt);

            }
			}
		}
	}

   // Inform user neighbour list calculation is complete
	if(vmpi::my_rank == 0){
		terminaltextcolor(GREEN);
		std::cout << "done!" << std::endl;
		terminaltextcolor(WHITE);
	}
   zlog << zTs() << "\tNeighbour list calculation complete"<< std::endl;

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

   zlog << zTs() << "\tSupercell array deallocated" << std::endl;

	return;
}

} // End of namespace neighbours
