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
// Vampire Header files
#include "create.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

namespace cs{
	
int create_crystal_structure(std::vector<cs::catom_t> & catom_array){
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "cs::create_crystal_structure has been called" << std::endl;}	

	// Calculate number of global and local unit cells required (rounding up)
	cs::total_num_unit_cells[0]=int(vmath::iceil(cs::system_dimensions[0]/unit_cell.dimensions[0]));
	cs::total_num_unit_cells[1]=int(vmath::iceil(cs::system_dimensions[1]/unit_cell.dimensions[1]));
	cs::total_num_unit_cells[2]=int(vmath::iceil(cs::system_dimensions[2]/unit_cell.dimensions[2]));

	int min_bounds[3];
	int max_bounds[3];
	
	#ifdef MPICF
	if(vmpi::mpi_mode==0){
		min_bounds[0] = int(vmpi::min_dimensions[0]/unit_cell.dimensions[0]);
		min_bounds[1] = int(vmpi::min_dimensions[1]/unit_cell.dimensions[1]);
		min_bounds[2] = int(vmpi::min_dimensions[2]/unit_cell.dimensions[2]);
		max_bounds[0] = vmath::iceil(vmpi::max_dimensions[0]/unit_cell.dimensions[0]);
		max_bounds[1] = vmath::iceil(vmpi::max_dimensions[1]/unit_cell.dimensions[1]);
		max_bounds[2] = vmath::iceil(vmpi::max_dimensions[2]/unit_cell.dimensions[2]);
	}
	else{
		min_bounds[0]=0;
		min_bounds[1]=0;
		min_bounds[2]=0;
		max_bounds[0]=cs::total_num_unit_cells[0];
		max_bounds[1]=cs::total_num_unit_cells[1];
		max_bounds[2]=cs::total_num_unit_cells[2];
	}
	#else
		min_bounds[0]=0;
		min_bounds[1]=0;
		min_bounds[2]=0;
		max_bounds[0]=cs::total_num_unit_cells[0];
		max_bounds[1]=cs::total_num_unit_cells[1];
		max_bounds[2]=cs::total_num_unit_cells[2];
	#endif
		
	cs::local_num_unit_cells[0]=max_bounds[0]-min_bounds[0];
	cs::local_num_unit_cells[1]=max_bounds[1]-min_bounds[1];
	cs::local_num_unit_cells[2]=max_bounds[2]-min_bounds[2];
	
	int num_atoms=cs::local_num_unit_cells[0]*cs::local_num_unit_cells[1]*cs::local_num_unit_cells[2]*unit_cell.atom.size();

	// set catom_array size
	catom_array.reserve(num_atoms);
	
	// Initialise atoms number
	int atom=0;

	unsigned int maxlh=2;

	// Duplicate unit cell
	for(int z=min_bounds[2];z<max_bounds[2];z++){
		for(int y=min_bounds[1];y<max_bounds[1];y++){
			for(int x=min_bounds[0];x<max_bounds[0];x++){
				
				// need to change this to accept non-orthogonal lattices
				// Loop over atoms in unit cell
				for(int uca=0;uca<unit_cell.atom.size();uca++){
					double cx = (double(x)+unit_cell.atom[uca].x)*unit_cell.dimensions[0];
					double cy = (double(y)+unit_cell.atom[uca].y)*unit_cell.dimensions[1];
					double cz = (double(z)+unit_cell.atom[uca].z)*unit_cell.dimensions[2];

					#ifdef MPICF
						if(vmpi::mpi_mode==0){
							// only generate atoms within allowed dimensions
							if(	(cx>=vmpi::min_dimensions[0] && cx<vmpi::max_dimensions[0]) &&
									(cy>=vmpi::min_dimensions[1] && cy<vmpi::max_dimensions[1]) &&
									(cz>=vmpi::min_dimensions[2] && cz<vmpi::max_dimensions[2])){
						#endif
							if((cx<cs::system_dimensions[0]) && (cy<cs::system_dimensions[1]) && (cz<cs::system_dimensions[2])){
							catom_array.push_back(cs::catom_t());
							catom_array[atom].x=cx;
							catom_array[atom].y=cy;
							catom_array[atom].z=cz;
							//std::cout << atom << "\t" << cx << "\t" << cy <<"\t" << cz << std::endl;
							catom_array[atom].material=unit_cell.atom[uca].mat;
							catom_array[atom].uc_id=uca;
							catom_array[atom].lh_category=unit_cell.atom[uca].hc+z*maxlh;
							catom_array[atom].uc_category=unit_cell.atom[uca].lc;
							catom_array[atom].scx=x;
							catom_array[atom].scy=y;
							catom_array[atom].scz=z;
							atom++;
							}
						#ifdef MPICF
							}
						}
						else{
							if((cx<cs::system_dimensions[0]) && (cy<cs::system_dimensions[1]) && (cz<cs::system_dimensions[2])){
							catom_array.push_back(cs::catom_t());
							catom_array[atom].x=cx;
							catom_array[atom].y=cy;
							catom_array[atom].z=cz;
							catom_array[atom].material=unit_cell.atom[uca].mat;
							catom_array[atom].uc_id=uca;
							catom_array[atom].lh_category=unit_cell.atom[uca].hc+z*maxlh;
							catom_array[atom].uc_category=unit_cell.atom[uca].lc;
							catom_array[atom].scx=x;
							catom_array[atom].scy=y;
							catom_array[atom].scz=z;
							catom_array[atom].include=false; // assume no atoms until classification complete
							atom++;
							}
						}
						#endif
					}
				}
			}
		}

	// Check to see if actual and expected number of atoms agree, if not trim the excess
	if(atom!=num_atoms){
		std::vector<cs::catom_t> tmp_catom_array(num_atoms);
		tmp_catom_array=catom_array;
		catom_array.resize(atom);
		for(int a=0;a<atom;a++){
			catom_array[a]=tmp_catom_array[a];
		}
		tmp_catom_array.resize(0);
	}
	

	// If z-height material selection is enabled then do so
	if(cs::SelectMaterialByZHeight==true){
		
		// determine z-bounds for materials
		std::vector<double> mat_min(mp::num_materials);
		std::vector<double> mat_max(mp::num_materials);

		for(int mat=0;mat<mp::num_materials;mat++){
			mat_min[mat]=mp::material[mat].min*cs::system_dimensions[2];
			mat_max[mat]=mp::material[mat].max*cs::system_dimensions[2];
			// alloys generally are not defined by height, and so have max = 0.0
			if(mat_max[mat]<0.0000001) mat_max[mat]=-0.1;
		}
		
		// Assign materials to generated atoms
		for(unsigned int atom=0;atom<catom_array.size();atom++){
			for(int mat=0;mat<mp::num_materials;mat++){
				const double cz=catom_array[atom].z;
				if((cz>=mat_min[mat]) && (cz<mat_max[mat])){
					catom_array[atom].material=mat;
					catom_array[atom].include=true;
				}
			}
		}

		// Delete unneeded atoms
		clear_atoms(catom_array);
		
	}
	
	// Check to see if any atoms have been generated
	if(atom==0){
		std::cout << "Error - no atoms have been generated, increase system dimensions!" << std::endl;
		zlog << zTs() << "Error: No atoms have been generated. Increase system dimensions." << std::endl;
		err::vexit();
	}

	// Now unselect all atoms by default for particle shape cutting
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		catom_array[atom].include=false;
	}

	return EXIT_SUCCESS;
}

void read_unit_cell(unit_cell_t & unit_cell, std::string filename){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::read_unit_cell has been called" << std::endl;}	

	std::cout << "Reading in unit cell data..." << std::flush;

	// ifstream declaration
	std::ifstream inputfile;
	
	// Open file read only
	inputfile.open(filename.c_str());
	
	// Check for opening
	if(!inputfile.is_open()){
		std::cerr << "Error! - cannot open unit cell input file: " << filename.c_str() << " Exiting" << std::endl;
		err::vexit();
	}

	// keep record of current line
	unsigned int line_counter=0;
	unsigned int line_id=0;
	// Loop over all lines
	while (! inputfile.eof() ){
		line_counter++;
		// read in whole line
		std::string line;
		getline(inputfile,line);
		//std::cout << line.c_str() << std::endl;

		// ignore blank lines
		std::string empty="";
		if(line==empty) continue;
		
		// set character triggers
		const char* hash="#";	// Comment identifier
		
		bool has_hash=false;
		// Determine if line is a comment line
		for(int i=0;i<line.length();i++){
			char c=line.at(i);

			if(c== *hash){
					has_hash=true;
					break;
			}
		}
		// if hash character found then read next line
		if(has_hash==true) continue;
		
		// convert line to string stream
		std::istringstream iss(line,std::istringstream::in);
		
		// defaults for interaction list
		int exc_type=-1; // assume isotropic
		int num_interactions=0; // assume no interactions
		int interaction_range=1; // assume +-1 unit cell as default
				
		// non-comment line found - check for line number
		switch(line_id){
			case 0:
				iss >> unit_cell.dimensions[0] >> unit_cell.dimensions[1] >> unit_cell.dimensions[2];
				break;
			case 1:
				iss >> unit_cell.shape[0][0] >> unit_cell.shape[0][1] >> unit_cell.shape[0][2];
				break;
			case 2:
				iss >> unit_cell.shape[1][0] >> unit_cell.shape[1][1] >> unit_cell.shape[1][2];
				break;
			case 3:
				iss >> unit_cell.shape[2][0] >> unit_cell.shape[2][1] >> unit_cell.shape[2][2];
				break;
			case 4:
				int num_uc_atoms;
				iss >> num_uc_atoms;
				//std::cout << "Reading in " << num_uc_atoms << " atoms" << std::endl;
				// resize unit_cell.atom array if within allowable bounds
				if( (num_uc_atoms >0) && (num_uc_atoms <= 1000000)) unit_cell.atom.resize(num_uc_atoms);
				else { std::cerr << "Error! Requested number of atoms " << num_uc_atoms << " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is outside of valid range 1-1,000,000. Exiting" << std::endl; err::vexit();}  
				// loop over all atoms and read into class
				for (int i=0; i<unit_cell.atom.size(); i++){
					line_counter++;
					// declare safe temporaries for atom input
					int id=i;
					double cx=2.0, cy=2.0,cz=2.0; // coordinates - default will give an error
					int mat_id=0, lcat_id=0, hcat_id=0; // sensible defaults if omitted
					// get line
					std::string atom_line;
					getline(inputfile,atom_line);
					std::istringstream atom_iss(atom_line,std::istringstream::in);
					atom_iss >> id >> cx >> cy >> cz >> mat_id >> lcat_id >> hcat_id;
					//std::cout << id << "\t" << cx << "\t" << cy << "\t" << cz<< "\t"  << mat_id << "\t" << lcat_id << "\t" << hcat_id << std::endl;
					//inputfile >> id >> cx >> cy >> cz >> mat_id >> lcat_id >> hcat_id;
					// now check for mostly sane input
					if(cx>=0.0 && cx <=1.0) unit_cell.atom[i].x=cx;
					else{std::cerr << "Error! atom x-coordinate for atom " << id << " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl; err::vexit();}
					if(cy>=0.0 && cy <=1.0) unit_cell.atom[i].y=cy;
					else{std::cerr << "Error! atom y-coordinate for atom " << id << " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl; err::vexit();}
					if(cz>=0.0 && cz <=1.0) unit_cell.atom[i].z=cz;
					else{std::cerr << "Error! atom z-coordinate for atom " << id << " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl; err::vexit();}
					if(mat_id >=0 && mat_id<mp::num_materials) unit_cell.atom[i].mat=mat_id;
					else{ std::cerr << "Error! Requested material id " << mat_id << "for atom number " << id <<  " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is outside of valid range 1-" << mp::num_materials << ". Exiting" << std::endl; err::vexit();} 
					unit_cell.atom[i].lc=lcat_id;
					unit_cell.atom[i].hc=hcat_id;
					//std::cout << i << "\t" << id << "\t" << cx << "\t" << cy << "\t" << cz << "\t" << mat_id << "\t" << lcat_id << "\t" << hcat_id << std::endl;
				}
				break;
			case 5:
				iss >> num_interactions >> exc_type;
				//std::cout << num_interactions << "\t" << exc_type << std::endl;
				if(num_interactions>=0) unit_cell.interaction.resize(num_interactions);
				else { std::cerr << "Error! Requested number of interactions " << num_interactions << " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is less than 0. Exiting" << std::endl; err::vexit();}
				// if exchange type omitted, then assume isotropic values from material file
				//if(exc_type==-1) unit_cell.exchange_type=0;
				// loop over all interactions and read into class
				for (int i=0; i<num_interactions; i++){
					//std::cout << "setting up interaction "<< i+1<< " of " << num_interactions << " interactions" << std::endl; 
					// declare safe temporaries for interaction input
					int id=i;
					int iatom=-1,jatom=-1; // atom pairs
					int dx=0, dy=0,dz=0; // relative unit cell coordinates
					// get line
					std::string int_line;
					getline(inputfile,int_line);
					//std::cout << int_line.c_str() << std::endl;
					std::istringstream int_iss(int_line,std::istringstream::in);
					int_iss >> id >> iatom >> jatom >> dx >> dy >> dz;
					//inputfile >> id >> iatom >> jatom >> dx >> dy >> dz;
					line_counter++;
					// check for sane input
					if(iatom>=0 && iatom < unit_cell.atom.size()) unit_cell.interaction[i].i=iatom;
					else{std::cerr << "Error! iatom number "<< iatom <<" for interaction id " << id << " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is outside of valid range 0-"<< unit_cell.atom.size()-1 << ". Exiting" << std::endl; err::vexit();}
					if(iatom>=0 && jatom < unit_cell.atom.size()) unit_cell.interaction[i].j=jatom;
					else{std::cerr << "Error! jatom number "<< jatom <<" for interaction id " << id << " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is outside of valid range 0-"<< unit_cell.atom.size()-1 << ". Exiting" << std::endl; err::vexit();}
					unit_cell.interaction[i].dx=dx;
					unit_cell.interaction[i].dy=dy;
					unit_cell.interaction[i].dz=dz;
					// check for long range interactions
					if(abs(dx)>interaction_range) interaction_range=abs(dx);
					if(abs(dy)>interaction_range) interaction_range=abs(dy);
					if(abs(dz)>interaction_range) interaction_range=abs(dz);
					
					int iatom_mat = unit_cell.atom[iatom].mat;
					int jatom_mat = unit_cell.atom[jatom].mat;
					switch(exc_type){
						case -1: // assume isotropic
							unit_cell.interaction[i].Jij[0][0]=mp::material[iatom_mat].Jij_matrix[jatom_mat];
							break;
						case 0:
							int_iss >> unit_cell.interaction[i].Jij[0][0];
							//std::cout << i << "\t" << unit_cell.interaction[i].Jij[0][0] << std::endl;
							break;
						case 1:
							int_iss >> unit_cell.interaction[i].Jij[0][0] >> unit_cell.interaction[i].Jij[1][1] >> unit_cell.interaction[i].Jij[2][2];
							break;
						case 2:
							int_iss >> unit_cell.interaction[i].Jij[0][0] >> unit_cell.interaction[i].Jij[0][1] >> unit_cell.interaction[i].Jij[0][2];
							int_iss >> unit_cell.interaction[i].Jij[1][0] >> unit_cell.interaction[i].Jij[1][1] >> unit_cell.interaction[i].Jij[1][2];
							int_iss >> unit_cell.interaction[i].Jij[2][0] >> unit_cell.interaction[i].Jij[2][1] >> unit_cell.interaction[i].Jij[2][2];
							break;
						default:
							std::cerr << "Error! Requested exchange type " << exc_type << " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is outside of valid range 0-2. Exiting" << std::endl; err::vexit();
					}
					// increment number of interactions for atom i
					unit_cell.atom[iatom].ni++;
				}
				// set interaction range
				unit_cell.interaction_range=interaction_range;
				// set exchange type
				unit_cell.exchange_type=exc_type;
				break;
			default:
				std::cerr << "Error! Unknown line type on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << ". Exiting" << std::endl; err::vexit();
		}
		line_id++;
	} // end of while loop
	
	zlog << zTs() << "Done!" << std::endl;
	zlog << zTs() << "\t" << "Number of atoms read-in: " << unit_cell.atom.size() << std::endl;
	zlog << zTs() << "\t" << "Number of interactions read-in: " << unit_cell.interaction.size() << std::endl;
	zlog << zTs() << "\t" << "Exchange type: " <<  unit_cell.exchange_type << std::endl;
	zlog << zTs() << "\t" << "Calculated interaction range: " << unit_cell.interaction_range << " Unit Cells" << std::endl;

	return;
}

void unit_cell_set(unit_cell_t & unit_cell){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::unit_cell_set has been called" << std::endl;}	

	// check for read-in of unit cell
	std::string blank="";
	if(cs::unit_cell_file.c_str()!=blank){
		read_unit_cell(unit_cell, cs::unit_cell_file);
		return;
	}
	
	// global values
	unit_cell.dimensions[0]=cs::unit_cell_size[0];
	unit_cell.dimensions[1]=cs::unit_cell_size[1];
	unit_cell.dimensions[2]=cs::unit_cell_size[2];
	
	unit_cell.exchange_type=-1;
	
	unit_cell.shape[0][0]=1.0;
	unit_cell.shape[0][1]=0.0;
	unit_cell.shape[0][2]=0.0;

	unit_cell.shape[1][0]=0.0;
	unit_cell.shape[1][1]=1.0;
	unit_cell.shape[1][2]=0.0;

	unit_cell.shape[2][0]=0.0;
	unit_cell.shape[2][1]=0.0;
	unit_cell.shape[2][2]=1.0;

		// Simple Cubic
		if(cs::crystal_structure=="sc"){
			unit_cell.lcsize=1;
			unit_cell.hcsize=1;
			unit_cell.interaction_range=1;
			unit_cell.atom.resize(1);
			unit_cell.surface_threshold=6;
			//-----------------------------
			unit_cell.atom[0].x=0.0;
			unit_cell.atom[0].y=0.0;
			unit_cell.atom[0].z=0.0;
			unit_cell.atom[0].lc=0;
			unit_cell.atom[0].hc=0;
			unit_cell.atom[0].ni=6;
			//-----------------------------
			unit_cell.interaction.resize(6);
			//-----------------------------
			unit_cell.interaction[0].i=0;
			unit_cell.interaction[0].j=0;
			unit_cell.interaction[0].dx=1;
			unit_cell.interaction[0].dy=0;
			unit_cell.interaction[0].dz=0;
			//-----------------------------
			unit_cell.interaction[1].i=0;
			unit_cell.interaction[1].j=0;
			unit_cell.interaction[1].dx=-1;
			unit_cell.interaction[1].dy=0;
			unit_cell.interaction[1].dz=0;
			//-----------------------------
			unit_cell.interaction[2].i=0;
			unit_cell.interaction[2].j=0;
			unit_cell.interaction[2].dx=0;
			unit_cell.interaction[2].dy=1;
			unit_cell.interaction[2].dz=0;
			//-----------------------------
			unit_cell.interaction[3].i=0;
			unit_cell.interaction[3].j=0;
			unit_cell.interaction[3].dx=0;
			unit_cell.interaction[3].dy=-1;
			unit_cell.interaction[3].dz=0;
			//-----------------------------
			unit_cell.interaction[4].i=0;
			unit_cell.interaction[4].j=0;
			unit_cell.interaction[4].dx=0;
			unit_cell.interaction[4].dy=0;
			unit_cell.interaction[4].dz=1;
			//-----------------------------
			unit_cell.interaction[5].i=0;
			unit_cell.interaction[5].j=0;
			unit_cell.interaction[5].dx=0;
			unit_cell.interaction[5].dy=0;
			unit_cell.interaction[5].dz=-1;
			//-----------------------------
		}
		else if(cs::crystal_structure=="bcc"){
			unit_cell.lcsize=2;
			unit_cell.hcsize=2;
			unit_cell.interaction_range=1;
			unit_cell.atom.resize(2);
			unit_cell.surface_threshold=8;
			//-----------------------------
			unit_cell.atom[0].x=0;
			unit_cell.atom[0].y=0;
			unit_cell.atom[0].z=0;
			unit_cell.atom[0].lc=0;
			unit_cell.atom[0].hc=0;
			unit_cell.atom[0].ni=8;
			//-----------------------------
			unit_cell.atom[1].x=0.5;
			unit_cell.atom[1].y=0.5;
			unit_cell.atom[1].z=0.5;
			unit_cell.atom[1].lc=1;
			unit_cell.atom[1].hc=1;
			unit_cell.atom[1].ni=8;
			//-----------------------------
			unit_cell.interaction.resize(16);
			//-----------------------------
			unit_cell.interaction[0].i=0;
			unit_cell.interaction[0].j=1;
			unit_cell.interaction[0].dx=-1;
			unit_cell.interaction[0].dy=-1;
			unit_cell.interaction[0].dz=-1;
			//-----------------------------
			unit_cell.interaction[1].i=0;
			unit_cell.interaction[1].j=1;
			unit_cell.interaction[1].dx=-1;
			unit_cell.interaction[1].dy=-1;
			unit_cell.interaction[1].dz=0;
			//-----------------------------
			unit_cell.interaction[2].i=0;
			unit_cell.interaction[2].j=1;
			unit_cell.interaction[2].dx=-1;
			unit_cell.interaction[2].dy=0;
			unit_cell.interaction[2].dz=-1;
			//-----------------------------
			unit_cell.interaction[3].i=0;
			unit_cell.interaction[3].j=1;
			unit_cell.interaction[3].dx=-1;
			unit_cell.interaction[3].dy=0;
			unit_cell.interaction[3].dz=0;
			//-----------------------------
			unit_cell.interaction[4].i=0;
			unit_cell.interaction[4].j=1;
			unit_cell.interaction[4].dx=0;
			unit_cell.interaction[4].dy=-1;
			unit_cell.interaction[4].dz=-1;
			//-----------------------------
			unit_cell.interaction[5].i=0;
			unit_cell.interaction[5].j=1;
			unit_cell.interaction[5].dx=0;
			unit_cell.interaction[5].dy=-1;
			unit_cell.interaction[5].dz=0;
			//-----------------------------
			unit_cell.interaction[6].i=0;
			unit_cell.interaction[6].j=1;
			unit_cell.interaction[6].dx=0;
			unit_cell.interaction[6].dy=0;
			unit_cell.interaction[6].dz=-1;
			//-----------------------------
			unit_cell.interaction[7].i=0;
			unit_cell.interaction[7].j=1;
			unit_cell.interaction[7].dx=0;
			unit_cell.interaction[7].dy=0;
			unit_cell.interaction[7].dz=0;
			//-----------------------------
			unit_cell.interaction[8].i=1;
			unit_cell.interaction[8].j=0;
			unit_cell.interaction[8].dx=0;
			unit_cell.interaction[8].dy=0;
			unit_cell.interaction[8].dz=0;
			//-----------------------------
			unit_cell.interaction[9].i=1;
			unit_cell.interaction[9].j=0;
			unit_cell.interaction[9].dx=0;
			unit_cell.interaction[9].dy=0;
			unit_cell.interaction[9].dz=1;
			//-----------------------------
			unit_cell.interaction[10].i=1;
			unit_cell.interaction[10].j=0;
			unit_cell.interaction[10].dx=0;
			unit_cell.interaction[10].dy=1;
			unit_cell.interaction[10].dz=0;
			//-----------------------------
			unit_cell.interaction[11].i=1;
			unit_cell.interaction[11].j=0;
			unit_cell.interaction[11].dx=0;
			unit_cell.interaction[11].dy=1;
			unit_cell.interaction[11].dz=1;
			//-----------------------------
			unit_cell.interaction[12].i=1;
			unit_cell.interaction[12].j=0;
			unit_cell.interaction[12].dx=1;
			unit_cell.interaction[12].dy=0;
			unit_cell.interaction[12].dz=0;
			//-----------------------------
			unit_cell.interaction[13].i=1;
			unit_cell.interaction[13].j=0;
			unit_cell.interaction[13].dx=1;
			unit_cell.interaction[13].dy=0;
			unit_cell.interaction[13].dz=1;
			//-----------------------------
			unit_cell.interaction[14].i=1;
			unit_cell.interaction[14].j=0;
			unit_cell.interaction[14].dx=1;
			unit_cell.interaction[14].dy=1;
			unit_cell.interaction[14].dz=0;
			//-----------------------------
			unit_cell.interaction[15].i=1;
			unit_cell.interaction[15].j=0;
			unit_cell.interaction[15].dx=1;
			unit_cell.interaction[15].dy=1;
			unit_cell.interaction[15].dz=1;
		}
		else if(cs::crystal_structure=="fct"){
			unit_cell.lcsize=2;
			unit_cell.hcsize=1;
			unit_cell.interaction_range=1;
			unit_cell.atom.resize(2);
			unit_cell.surface_threshold=4;
			//-----------------------------
			unit_cell.atom[0].x=0.0;
			unit_cell.atom[0].y=0.0;
			unit_cell.atom[0].z=0.0;
			unit_cell.atom[0].lc=0;
			unit_cell.atom[0].hc=0;
			unit_cell.atom[0].ni=4;
			//-----------------------------
			unit_cell.atom[1].x=0.5;
			unit_cell.atom[1].y=0.5;
			unit_cell.atom[1].z=0.0;
			unit_cell.atom[1].lc=1;
			unit_cell.atom[1].hc=0;
			unit_cell.atom[1].ni=4;
			//-----------------------------
		}
		else if(cs::crystal_structure=="fcc"){
			unit_cell.lcsize=4;
			unit_cell.hcsize=2;
			unit_cell.interaction_range=1;
			unit_cell.atom.resize(4);
			unit_cell.surface_threshold=12;
			//-----------------------------
			unit_cell.atom[0].x=0.0;
			unit_cell.atom[0].y=0.0;
			unit_cell.atom[0].z=0.0;
			unit_cell.atom[0].lc=0;
			unit_cell.atom[0].hc=0;
			unit_cell.atom[0].ni=12;
			//-----------------------------
			unit_cell.atom[1].x=0.5;
			unit_cell.atom[1].y=0.5;
			unit_cell.atom[1].z=0.0;
			unit_cell.atom[1].lc=1;
			unit_cell.atom[1].hc=0;
			unit_cell.atom[1].ni=12;
			//-----------------------------
			unit_cell.atom[2].x=0.5;
			unit_cell.atom[2].y=0.0;
			unit_cell.atom[2].z=0.5;
			unit_cell.atom[2].lc=2;
			unit_cell.atom[2].hc=1;
			unit_cell.atom[2].ni=12;
			//-----------------------------
			unit_cell.atom[3].x=0.0;
			unit_cell.atom[3].y=0.5;
			unit_cell.atom[3].z=0.5;
			unit_cell.atom[3].lc=3;
			unit_cell.atom[3].hc=1;
			unit_cell.atom[3].ni=12;
			//-----------------------------
			unit_cell.interaction.resize(48);
			//-----------------------------
			unit_cell.interaction[0].i=0;
			unit_cell.interaction[0].j=1;
			unit_cell.interaction[0].dx=-1;
			unit_cell.interaction[0].dy=-1;
			unit_cell.interaction[0].dz=0;
			//----------------------------------------------
			unit_cell.interaction[1].i=0;
			unit_cell.interaction[1].j=2;
			unit_cell.interaction[1].dx=-1;
			unit_cell.interaction[1].dy=0;
			unit_cell.interaction[1].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[2].i=0;
			unit_cell.interaction[2].j=1;
			unit_cell.interaction[2].dx=-1;
			unit_cell.interaction[2].dy=0;
			unit_cell.interaction[2].dz=0;
			//----------------------------------------------
			unit_cell.interaction[3].i=0;
			unit_cell.interaction[3].j=2;
			unit_cell.interaction[3].dx=-1;
			unit_cell.interaction[3].dy=0;
			unit_cell.interaction[3].dz=0;
			//----------------------------------------------
			unit_cell.interaction[4].i=0;
			unit_cell.interaction[4].j=3;
			unit_cell.interaction[4].dx=0;
			unit_cell.interaction[4].dy=-1;
			unit_cell.interaction[4].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[5].i=0;
			unit_cell.interaction[5].j=1;
			unit_cell.interaction[5].dx=0;
			unit_cell.interaction[5].dy=-1;
			unit_cell.interaction[5].dz=0;
			//----------------------------------------------
			unit_cell.interaction[6].i=0;
			unit_cell.interaction[6].j=3;
			unit_cell.interaction[6].dx=0;
			unit_cell.interaction[6].dy=-1;
			unit_cell.interaction[6].dz=0;
			//----------------------------------------------
			unit_cell.interaction[7].i=0;
			unit_cell.interaction[7].j=2;
			unit_cell.interaction[7].dx=0;
			unit_cell.interaction[7].dy=0;
			unit_cell.interaction[7].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[8].i=0;
			unit_cell.interaction[8].j=3;
			unit_cell.interaction[8].dx=0;
			unit_cell.interaction[8].dy=0;
			unit_cell.interaction[8].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[9].i=0;
			unit_cell.interaction[9].j=1;
			unit_cell.interaction[9].dx=0;
			unit_cell.interaction[9].dy=0;
			unit_cell.interaction[9].dz=0;
			//----------------------------------------------
			unit_cell.interaction[10].i=0;
			unit_cell.interaction[10].j=2;
			unit_cell.interaction[10].dx=0;
			unit_cell.interaction[10].dy=0;
			unit_cell.interaction[10].dz=0;
			//----------------------------------------------
			unit_cell.interaction[11].i=0;
			unit_cell.interaction[11].j=3;
			unit_cell.interaction[11].dx=0;
			unit_cell.interaction[11].dy=0;
			unit_cell.interaction[11].dz=0;
			//----------------------------------------------
			unit_cell.interaction[12].i=1;
			unit_cell.interaction[12].j=2;
			unit_cell.interaction[12].dx=0;
			unit_cell.interaction[12].dy=0;
			unit_cell.interaction[12].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[13].i=1;
			unit_cell.interaction[13].j=3;
			unit_cell.interaction[13].dx=0;
			unit_cell.interaction[13].dy=0;
			unit_cell.interaction[13].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[14].i=1;
			unit_cell.interaction[14].j=0;
			unit_cell.interaction[14].dx=0;
			unit_cell.interaction[14].dy=0;
			unit_cell.interaction[14].dz=0;
			//----------------------------------------------
			unit_cell.interaction[15].i=1;
			unit_cell.interaction[15].j=2;
			unit_cell.interaction[15].dx=0;
			unit_cell.interaction[15].dy=0;
			unit_cell.interaction[15].dz=0;
			//----------------------------------------------
			unit_cell.interaction[16].i=1;
			unit_cell.interaction[16].j=3;
			unit_cell.interaction[16].dx=0;
			unit_cell.interaction[16].dy=0;
			unit_cell.interaction[16].dz=0;
			//----------------------------------------------
			unit_cell.interaction[17].i=1;
			unit_cell.interaction[17].j=2;
			unit_cell.interaction[17].dx=0;
			unit_cell.interaction[17].dy=1;
			unit_cell.interaction[17].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[18].i=1;
			unit_cell.interaction[18].j=0;
			unit_cell.interaction[18].dx=0;
			unit_cell.interaction[18].dy=1;
			unit_cell.interaction[18].dz=0;
			//----------------------------------------------
			unit_cell.interaction[19].i=1;
			unit_cell.interaction[19].j=2;
			unit_cell.interaction[19].dx=0;
			unit_cell.interaction[19].dy=1;
			unit_cell.interaction[19].dz=0;
			//----------------------------------------------
			unit_cell.interaction[20].i=1;
			unit_cell.interaction[20].j=3;
			unit_cell.interaction[20].dx=1;
			unit_cell.interaction[20].dy=0;
			unit_cell.interaction[20].dz=-1;
			//----------------------------------------------
			unit_cell.interaction[21].i=1;
			unit_cell.interaction[21].j=0;
			unit_cell.interaction[21].dx=1;
			unit_cell.interaction[21].dy=0;
			unit_cell.interaction[21].dz=0;
			//----------------------------------------------
			unit_cell.interaction[22].i=1;
			unit_cell.interaction[22].j=3;
			unit_cell.interaction[22].dx=1;
			unit_cell.interaction[22].dy=0;
			unit_cell.interaction[22].dz=0;
			//----------------------------------------------
			unit_cell.interaction[23].i=1;
			unit_cell.interaction[23].j=0;
			unit_cell.interaction[23].dx=1;
			unit_cell.interaction[23].dy=1;
			unit_cell.interaction[23].dz=0;
			//----------------------------------------------
			unit_cell.interaction[24].i=2;
			unit_cell.interaction[24].j=1;
			unit_cell.interaction[24].dx=0;
			unit_cell.interaction[24].dy=-1;
			unit_cell.interaction[24].dz=0;
			//----------------------------------------------
			unit_cell.interaction[25].i=2;
			unit_cell.interaction[25].j=3;
			unit_cell.interaction[25].dx=0;
			unit_cell.interaction[25].dy=-1;
			unit_cell.interaction[25].dz=0;
			//----------------------------------------------
			unit_cell.interaction[26].i=2;
			unit_cell.interaction[26].j=1;
			unit_cell.interaction[26].dx=0;
			unit_cell.interaction[26].dy=-1;
			unit_cell.interaction[26].dz=1;
			//----------------------------------------------
			unit_cell.interaction[27].i=2;
			unit_cell.interaction[27].j=0;
			unit_cell.interaction[27].dx=0;
			unit_cell.interaction[27].dy=0;
			unit_cell.interaction[27].dz=0;
			//----------------------------------------------
			unit_cell.interaction[28].i=2;
			unit_cell.interaction[28].j=1;
			unit_cell.interaction[28].dx=0;
			unit_cell.interaction[28].dy=0;
			unit_cell.interaction[28].dz=0;
			//----------------------------------------------
			unit_cell.interaction[29].i=2;
			unit_cell.interaction[29].j=3;
			unit_cell.interaction[29].dx=0;
			unit_cell.interaction[29].dy=0;
			unit_cell.interaction[29].dz=0;
			//----------------------------------------------
			unit_cell.interaction[30].i=2;
			unit_cell.interaction[30].j=0;
			unit_cell.interaction[30].dx=0;
			unit_cell.interaction[30].dy=0;
			unit_cell.interaction[30].dz=1;
			//----------------------------------------------
			unit_cell.interaction[31].i=2;
			unit_cell.interaction[31].j=1;
			unit_cell.interaction[31].dx=0;
			unit_cell.interaction[31].dy=0;
			unit_cell.interaction[31].dz=1;
			//----------------------------------------------
			unit_cell.interaction[32].i=2;
			unit_cell.interaction[32].j=3;
			unit_cell.interaction[32].dx=1;
			unit_cell.interaction[32].dy=-1;
			unit_cell.interaction[32].dz=0;
			//----------------------------------------------
			unit_cell.interaction[33].i=2;
			unit_cell.interaction[33].j=0;
			unit_cell.interaction[33].dx=1;
			unit_cell.interaction[33].dy=0;
			unit_cell.interaction[33].dz=0;
			//----------------------------------------------
			unit_cell.interaction[34].i=2;
			unit_cell.interaction[34].j=3;
			unit_cell.interaction[34].dx=1;
			unit_cell.interaction[34].dy=0;
			unit_cell.interaction[34].dz=0;
			//----------------------------------------------
			unit_cell.interaction[35].i=2;
			unit_cell.interaction[35].j=0;
			unit_cell.interaction[35].dx=1;
			unit_cell.interaction[35].dy=0;
			unit_cell.interaction[35].dz=1;
			//----------------------------------------------
			unit_cell.interaction[36].i=3;
			unit_cell.interaction[36].j=1;
			unit_cell.interaction[36].dx=-1;
			unit_cell.interaction[36].dy=0;
			unit_cell.interaction[36].dz=0;
			//----------------------------------------------
			unit_cell.interaction[37].i=3;
			unit_cell.interaction[37].j=2;
			unit_cell.interaction[37].dx=-1;
			unit_cell.interaction[37].dy=0;
			unit_cell.interaction[37].dz=0;
			//----------------------------------------------
			unit_cell.interaction[38].i=3;
			unit_cell.interaction[38].j=1;
			unit_cell.interaction[38].dx=-1;
			unit_cell.interaction[38].dy=0;
			unit_cell.interaction[38].dz=1;
			//----------------------------------------------
			unit_cell.interaction[39].i=3;
			unit_cell.interaction[39].j=2;
			unit_cell.interaction[39].dx=-1;
			unit_cell.interaction[39].dy=1;
			unit_cell.interaction[39].dz=0;
			//----------------------------------------------
			unit_cell.interaction[40].i=3;
			unit_cell.interaction[40].j=0;
			unit_cell.interaction[40].dx=0;
			unit_cell.interaction[40].dy=0;
			unit_cell.interaction[40].dz=0;
			//----------------------------------------------
			unit_cell.interaction[41].i=3;
			unit_cell.interaction[41].j=1;
			unit_cell.interaction[41].dx=0;
			unit_cell.interaction[41].dy=0;
			unit_cell.interaction[41].dz=0;
			//----------------------------------------------
			unit_cell.interaction[42].i=3;
			unit_cell.interaction[42].j=2;
			unit_cell.interaction[42].dx=0;
			unit_cell.interaction[42].dy=0;
			unit_cell.interaction[42].dz=0;
			//----------------------------------------------
			unit_cell.interaction[43].i=3;
			unit_cell.interaction[43].j=0;
			unit_cell.interaction[43].dx=0;
			unit_cell.interaction[43].dy=0;
			unit_cell.interaction[43].dz=1;
			//----------------------------------------------
			unit_cell.interaction[44].i=3;
			unit_cell.interaction[44].j=1;
			unit_cell.interaction[44].dx=0;
			unit_cell.interaction[44].dy=0;
			unit_cell.interaction[44].dz=1;
			//----------------------------------------------
			unit_cell.interaction[45].i=3;
			unit_cell.interaction[45].j=0;
			unit_cell.interaction[45].dx=0;
			unit_cell.interaction[45].dy=1;
			unit_cell.interaction[45].dz=0;
			//----------------------------------------------
			unit_cell.interaction[46].i=3;
			unit_cell.interaction[46].j=2;
			unit_cell.interaction[46].dx=0;
			unit_cell.interaction[46].dy=1;
			unit_cell.interaction[46].dz=0;
			//----------------------------------------------
			unit_cell.interaction[47].i=3;
			unit_cell.interaction[47].j=0;
			unit_cell.interaction[47].dx=0;
			unit_cell.interaction[47].dy=1;
			unit_cell.interaction[47].dz=1;
			//----------------------------------------------

		}
		else{
			std::cerr << "Error -  unknown crystal_type" << std::endl; 
			err::vexit();
		}

	}

} // end of namespace cs

