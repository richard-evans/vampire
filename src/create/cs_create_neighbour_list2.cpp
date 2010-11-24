//======================================================================
//                         create_neighbourlist
//                Subroutine to create neighbour list
//
//======================================================================
#include "material.hpp"
#include "errors.hpp"
#include "create.hpp"
#include "vmpi.hpp"
#include <cmath>
#include <iostream>

namespace cs{
int create_neighbourlist(std::vector<cs::catom_t> & catom_array, std::vector<std::vector <int> > & cneighbourlist){
	//====================================================================================
	//
	//								cs_create_neighbourlist
	//
	//						Subroutine to calculate neighbour list
	//
	//							Version 1.0 R Evans 22/09/2008
	//
	//====================================================================================
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::create_neighbourlist has been called" << std::endl;}	
	
	const int num_atoms = catom_array.size();
	const int max_atoms_per_supercell=4;
	int range=1;	// Range of exchange interaction (unit cells)
	//-----------------------------------------------------
	// Setup neighbourlist
	//-----------------------------------------------------

	// Reserve space for num_atoms
	cneighbourlist.reserve(num_atoms);

	// Reserve space for each atom in neighbour list according to material type
	for(int atom=0;atom<num_atoms;atom++){
		int max_nn = mp::material[catom_array[atom].material].hamiltonian_num_neighbours;
		cneighbourlist.push_back(std::vector<int>());
		cneighbourlist[atom].reserve(max_nn);
	}

	// Calculate system dimensions and number of supercells
	double min[3]={1.0e20,1.0e20,1.0e20};
	double max[3]={-1.0e20,-1.0e20,-1.0e20};
	
	for(int atom=0;atom<num_atoms;atom++){
		double c[3]={catom_array[atom].x,catom_array[atom].y,catom_array[atom].z};
		for(int i=0;i<3;i++){
			if(c[i]<min[i]){
				min[i]=c[i];
			}
			if(c[i]>max[i]){
				max[i]=c[i];
			}
		}
	}
	// Added small correction to avoid rounding errors - RF 8/6/2010
	const int d[3]={1+round((max[0]-min[0])/mp::lattice_constant[0]+0.001),1+round((max[1]-min[1])/mp::lattice_constant[1]+0.001),1+round((max[2]-min[2])/mp::lattice_constant[2]+0.001)};
	
	// offset in whole unit cells
	const int offset[3] = {int((min[0])/mp::lattice_constant[0]), int((min[1])/mp::lattice_constant[1]), int((min[2])/mp::lattice_constant[2])};
	
		// Declare array for create space for 3D supercell array
		//int supercell_array[d[0]][d[1]][d[2]][max_atoms_per_supercell+2];
		int**** supercell_array;
		std::cout << "Memory required for neighbourlist calculation:" << 8.0*double(d[0])*double(d[1])*double(d[2])*double(max_atoms_per_supercell+2)/1.0e6 << " MB" << std::endl;
		try{supercell_array=new int***[d[0]];
			for(int i=0; i<d[0] ; i++){
				supercell_array[i]=new int**[d[1]];
				for(int j=0; j<d[1] ; j++){
					supercell_array[i][j]=new int*[d[2]];
					for(int k=0; k<d[2] ; k++){
						supercell_array[i][j][k]=new int[max_atoms_per_supercell+2];
					}
				}
			}
		}
		catch(...){std::cout << "error allocating supercell_array" << std::endl;exit(EXIT_FAILURE);}
		
		const int num_cells=d[0]*d[1]*d[2];
		std::vector< std::vector <int> > cell_coord_array;
		cell_coord_array.reserve(num_cells);
		for(int i=0;i<num_cells;i++){
			cell_coord_array.push_back(std::vector<int>());
			cell_coord_array[i].resize(3);
		}
		int cell=0;

		// Initialise array
		for(int z=0;z<d[2];z++){
			for(int y=0;y<d[1];y++){
				for(int x=0;x<d[0];x++){
					for(int n=0;n<max_atoms_per_supercell;n++){
						supercell_array[x][y][z][n]=0;
					}
					supercell_array[x][y][z][max_atoms_per_supercell]=0;
					supercell_array[x][y][z][max_atoms_per_supercell+1]=cell;
					//save cell coordinates
					cell_coord_array[cell][0]=x;
					cell_coord_array[cell][1]=y;
					cell_coord_array[cell][2]=z;
					cell++;
				}
			}
		}

		double atom_offset[3]={0.25*mp::lattice_constant[0],0.25*mp::lattice_constant[1],0.25*mp::lattice_constant[2]};

                // Populate supercell array with atom numbers                                                                                                                     
                for(int atom=0;atom<num_atoms;atom++){
		  double c[3]={catom_array[atom].x+atom_offset[0],catom_array[atom].y+atom_offset[1],catom_array[atom].z+atom_offset[2]};
			int scc[3]={0,0,0}; // super cell coordinates
			for(int i=0;i<3;i++){
				scc[i]=int(c[i]/mp::lattice_constant[i])-offset[i]; // Always round down for supercell coordinates
				// Always check cell in range
				if(scc[i]<0 || scc[i]>= d[i]){
				        std::cerr << "Error - atom out of supercell range in neighbourlist calculation!" << std::endl;
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
					std::cerr << "\tCell offest (dp): " << offset[0]*mp::lattice_constant[0] << "\t" << offset[1]*mp::lattice_constant[1] << "\t" << offset[2]*mp::lattice_constant[2] << std::endl;
					#ifdef MPICF
					MPI::COMM_WORLD.Abort(EXIT_FAILURE);
					exit(EXIT_FAILURE);
                                        #else
					exit(EXIT_FAILURE);
                                        #endif
				}
			}
			
			// Add atom to supercell
			int current=supercell_array[scc[0]][scc[1]][scc[2]][max_atoms_per_supercell];
			// Check for atoms greater than max_atoms_per_supercell
			if(current<max_atoms_per_supercell){
				supercell_array[scc[0]][scc[1]][scc[2]][current]=atom;
				supercell_array[scc[0]][scc[1]][scc[2]][max_atoms_per_supercell]++;
			}
			else{
				std::cerr << "Error, number of atoms per supercell exceeded, reduce system lattice constant" << std::endl;
				std::cerr << "\tAtom number:      " << atom << std::endl;
				std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
				std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
				std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
				std::cerr << "\tCell offset:      " << offset[0] << "\t" << offset[1] << "\t" << offset[2] << std::endl;
				#ifdef MPICF
				MPI::COMM_WORLD.Abort(EXIT_FAILURE);
				exit(EXIT_FAILURE);
                                #else
				exit(EXIT_FAILURE);
                                #endif

			}
		}
		
		// Generate neighbour list
		std::cout <<"Generating Neighbour list"; 
		//bool periodicity[3]={false,false,false};
		//int periodic_index[3];
		double periodic_offset[3];
	
		// Loop over all cells
		for(int cell=0;cell<num_cells;cell++){
			if(cell%(num_cells/10+1)==0){
				std::cout << "." << std::flush;
			}
			//std::cout << "cell: " << cell << ":>" << std::flush; // << std::endl;
			int scc[3]={cell_coord_array[cell][0],cell_coord_array[cell][1],cell_coord_array[cell][2]};
			const int num_atoms_in_cell=supercell_array[scc[0]][scc[1]][scc[2]][max_atoms_per_supercell];
			//loop over atoms in local cell
			for(int i=0;i<num_atoms_in_cell;i++){
				int atom=supercell_array[scc[0]][scc[1]][scc[2]][i];
				//std::cout << atom << "\t" << catom_array[atom].material << std::endl;
				//std::cout << "----------------------------------------------------------" << std::endl;
				//loop over neighbouring unit cells
				for(int z=scc[2]-range;z<=scc[2]+range;z++){
					if(z>=0 && z<d[2]){
						for(int y=scc[1]-range;y<=scc[1]+range;y++){
							//check for out of bounds cells
							if(y>=0 && y<d[1]){
								for(int x=scc[0]-range;x<=scc[0]+range;x++){
									//check for out of bounds cells
									if(x>=0 && x<d[0]){
										// loop over atoms in my own and neighboring cells
										const int num_atoms_in_neighbour_cell=supercell_array[x][y][z][max_atoms_per_supercell];
										for(int j=0;j<num_atoms_in_neighbour_cell;j++){
											int natom=supercell_array[x][y][z][j];
											// exclude local atom for distance calculation
											if(atom!=natom){
												double dx=catom_array[natom].x-catom_array[atom].x;
												double dy=catom_array[natom].y-catom_array[atom].y;
												double dz=catom_array[natom].z-catom_array[atom].z+periodic_offset[2];
												double drange=mp::material[catom_array[atom].material].cutoff*mp::lattice_constant[0];
												if(dx*dx+dy*dy+dz*dz<=drange*drange){
													if(cneighbourlist[atom].size()<=cneighbourlist[atom].capacity()){
														cneighbourlist[atom].push_back(natom);
													}
													else{
														std::cerr << "Error - number of neighbour atoms exceeded, increase hamiltonian num neighbours" << std::endl;
														std::cerr << "Size: " << cneighbourlist[atom].size() << "Capacity: " << cneighbourlist[atom].capacity() << std::endl;
														exit(EXIT_FAILURE);
													}
													//std::cout << "\t" << natom << "\t" << dx << "\t" << dy << "\t" << dz << "\t" <<dx*dx+dy*dy+dz*dz << "\t" << drange*drange << std::endl;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		std::cout << "done!" << std::endl;

		// Deallocate supercell array
		try{
			for(int i=0; i<d[0] ; i++){
				for(int j=0; j<d[1] ;j++){
					for(int k=0; k<d[2] ;k++){
						delete [] supercell_array[i][j][k];
					}
					delete [] supercell_array[i][j];
				}
				delete [] supercell_array[i];
			}
		delete [] supercell_array;
		supercell_array=NULL;
		}
	catch(...){std::cout << "error deallocating supercell_array" << std::endl; exit(1);}

	// Mark surface atoms
	//for(int atom=0;atom<num_atoms;atom++){
	//	if(int(cneighbourlist[atom].size())!=mp::material[catom_array[atom].material].hamiltonian_num_neighbours){
	//		catom_array[atom].material=2;
	//	}
	//}
	
	return EXIT_SUCCESS;
}

} // End of namespcae cs
