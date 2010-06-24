///
/// @file
/// @brief Contains vin and vout namespaces for file input.output in vampire. 
///
/// @details File and screen input and output are controlled via the separate namespaces.
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
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    15/01/2010
/// @internal
///	Created:		15/01/2010
///	Revision:	  ---
///=====================================================================================
///

// Headers
#include "atoms.hpp"
//#include "create_voronoi.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "units.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "stats.hpp"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

//==========================================================
// Global Output Streams
//==========================================================
std::ofstream vinfo("info");
std::ofstream vdp("vdp");
std::ofstream vmag("vmag");

/// @namespace
/// @brief Contains variables and functions for writing out program data.
/// 
/// @internal
///=====================================================================================
///
namespace vout{
/*	
	int scr(std::stringstream buff){
		if(num_processors=1){
			std::cout << buff; 
		else if(my_rank==0){
			std::cout << buff;
		}
		
		return EXIT_SUCCESS;
	}

} // end of namespace vout  */

//==========================================================
// Namespace output
//==========================================================
namespace output{
	int output_flags[10];
	int output_inc[10];
}

namespace pov{
	int counter=0;
}


int pov_file(){
//-----------------------------------------------------------------------
//
//
//
//-----------------------------------------------------------------------
		using pov::counter;

		#ifdef MPICF
		const int num_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
		#else
		const int num_atoms = atoms::num_atoms;
		#endif
		
		std::stringstream pov_file_sstr;
		pov_file_sstr << "spins.";
		pov_file_sstr << std::setfill('0') << std::setw(3) << vmpi::my_rank;
		pov_file_sstr << "." << std::setfill('0') << std::setw(5) << counter;
		pov_file_sstr << ".pin";
		//pov_file_sstr << "spins." << mpi_generic::my_rank << "." << counter << ".pov";
		std::string pov_file = pov_file_sstr.str();
		const char* pov_filec = pov_file.c_str();
		std::ofstream pov_file_ofstr;

		if(vmpi::my_rank==0){
			std::stringstream pov_hdr_sstr;
			pov_hdr_sstr << "spins." << std::setfill('0') << std::setw(3) << counter << ".pov";
			std::string pov_hdr = pov_hdr_sstr.str();
			const char* pov_hdrc = pov_hdr.c_str();
			pov_file_ofstr.open (pov_hdrc);
	
			double size, mag_vec;
			double vec[3];

			size = sqrt(material_parameters::system_dimensions[0]*material_parameters::system_dimensions[0] +
					material_parameters::system_dimensions[1]*material_parameters::system_dimensions[1] +
					material_parameters::system_dimensions[2]*material_parameters::system_dimensions[2]);

			vec[0] = (1.0/material_parameters::system_dimensions[0]);
			vec[1] = (1.0/material_parameters::system_dimensions[1]);
			vec[2] = (1.0/material_parameters::system_dimensions[2]);
			mag_vec = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
			vec[0]/=mag_vec;
			vec[1]/=mag_vec;
			vec[2]/=mag_vec;

			//---------------------------------------------------
			// Output file header (rank 0)
			//---------------------------------------------------
			pov_file_ofstr << "#include \"colors.inc\"" << std::endl;
			pov_file_ofstr << "#include \"metals.inc\""	<< std::endl;
			pov_file_ofstr << "#include \"screen.inc\""	<< std::endl;
			pov_file_ofstr << "#declare LX=" << material_parameters::system_dimensions[0]*0.5 << ";" << std::endl;
			pov_file_ofstr << "#declare LY=" << material_parameters::system_dimensions[1]*0.5 << ";" << std::endl;
			pov_file_ofstr << "#declare LZ=" << material_parameters::system_dimensions[2]*0.5 << ";" << std::endl;
			pov_file_ofstr << "#declare CX=" << size*vec[0]*6.0 << ";" << std::endl;
			pov_file_ofstr << "#declare CY=" << size*vec[1]*6.0 << ";" << std::endl;
			pov_file_ofstr << "#declare CZ=" << size*vec[2]*6.0 << ";" << std::endl;
			//pov_file_ofstr << "#declare MX=" << 2.0*mp::system_dimensions[0] << ";" << std::endl;
			//pov_file_ofstr << "#declare MY=" << -1.0*mp::system_dimensions[1] << ";" << std::endl;
			//pov_file_ofstr << "#declare MZ=" << mp::system_dimensions[2]*0.5 << ";" << std::endl;
	 		pov_file_ofstr << "#declare ref=0.4;" << std::endl;
	 		pov_file_ofstr << "#declare sscale=2.0;" << std::endl;
			pov_file_ofstr << "background { color Gray30 }" << std::endl;
			// for macrospins
			//pov_file_ofstr << "#declare LX=" <<  mp::system_dimensions[0] <<";" <<  std::endl;
			//pov_file_ofstr << "#declare LY=" <<  0.0 << ";" << std::endl;
			//pov_file_ofstr << "#declare CX=" << (size*vec[0]*12.0) << ";" << std::endl;
			//pov_file_ofstr << "#declare CY=" << (size*vec[1]*12.0)+mp::system_dimensions[1] << ";" << std::endl;
			//pov_file_ofstr << "#declare CZ=" << size*vec[2]*5.0 << ";" << std::endl;

			pov_file_ofstr << "Set_Camera(<CX,CY,CZ>, <LX,LY,LZ>, 15)" << std::endl;
			pov_file_ofstr << "Set_Camera_Aspect(4,3)" << std::endl;
			pov_file_ofstr << "Set_Camera_Sky(<0,0,1>)" << std::endl;
			pov_file_ofstr << "light_source { <2*CX, 2*CY, 2*CZ> color White}" << std::endl;

			// Axes
			//pov_file_ofstr << "cylinder{ <-100000,0,0>,<100000,0,0>," << mp::system_dimensions[2]*0.01 << " texture { pigment {color rgb <0,0,0>} finish {reflection {ref} diffuse 1 ambient 0}}}" << std::endl;
                        //pov_file_ofstr << "cylinder{ <0,-100000,0>,<0,100000,0>," << mp::system_dimensions[2]*0.01 << " texture { pigment {color rgb <0,0,0>} finish {reflection {ref} diffuse 1 ambient 0}}}" << std::endl;
                        //pov_file_ofstr << "cylinder{ <0,0,-100000>,<0,0,100000>," << mp::system_dimensions[2]*0.01 << " texture { pigment {color rgb <0,0,0>} finish {reflection {ref} diffuse 1 ambient 0}}}" << std::endl;
			// Macrospin output
			/*for(int mat=0; mat<mp::num_materials;mat++){
			  double mx =stats::sublattice_magm_array[mat]*stats::sublattice_mx_array[mat];
                          double my =stats::sublattice_magm_array[mat]*stats::sublattice_my_array[mat];
                          double mz =stats::sublattice_magm_array[mat]*stats::sublattice_mz_array[mat];
			  double red,green,blue,ireal;
			  ireal = mz;

			  if(ireal>0.8){
			    red = 0.0;
			    green = 0.0;
			    blue = 1.0;
			  }
			  else if(ireal>=0.0){
			    red = 1.0-ireal*1.2;
			    green = 1.0-ireal*1.2;
			    blue = 1.0;
			  }
			  else if(ireal>=-0.8){
			    red = 1.0;
			    green = 1.0+ireal*1.2;
			    blue = 1.0+ireal*1.2;
			  }
			  else if(ireal<-0.8){
			    red = 1.0;
			    green = 0.0;
			    blue = 0.0;
			  }
			  else{
			    red = 1.0;
			    green = 1.0;
			    blue = 1.0;
			  }

			  if(blue<0.0) blue=0.0;
			  if(red<0.0) red=0.0;
			  if(green<0.0) green=0.0;

			  pov_file_ofstr << "cone {<MX+" << mx*0.5*mp::system_dimensions[2] << ","
					 << "MY+" << my*0.5*mp::system_dimensions[2] << ","
                                         << "MZ+" << mz*0.5*mp::system_dimensions[2] << ">,0.0 <"
					 << "MX+" << mx*0.5*0.75*mp::system_dimensions[2] << ","
					 << "MY+" << my*0.5*0.75*mp::system_dimensions[2] << ","
					 << "MZ+" << mz*0.5*0.75*mp::system_dimensions[2] << ">," << mp::system_dimensions[2]*0.125
					 << "texture { pigment {color rgb <" << red << " " << green << " " << blue << ">}"
					 << "finish {reflection {ref} diffuse 1 ambient 0}}}" << std::endl;			  
                          pov_file_ofstr << "cylinder {<MX,MY,MZ>,<"
                                         << "MX+" << mx*0.5*0.75*mp::system_dimensions[2] << ","
                                         << "MY+" << my*0.5*0.75*mp::system_dimensions[2] << ","
                                         << "MZ+" << mz*0.5*0.75*mp::system_dimensions[2] << ">," << mp::system_dimensions[2]*0.066
                                         << "texture { pigment {color rgb <" << red << " " << green << " " << blue << ">}"
                                         << "finish {reflection {ref} diffuse 1 ambient 0}}}" << std::endl;
			}*/

			for(int p =0;p<vmpi::num_processors;p++){
				std::stringstream pov_sstr;
				pov_sstr << "spins." << std::setfill('0') << std::setw(3) << p << "." << std::setfill('0') << std::setw(5) << counter << ".pin";
				pov_file_ofstr << "#include \"" << pov_sstr.str() << "\"" << std::endl;
				//pov_sstr << "spins." << p << "." << counter << ".pov";
				//pov_file_ofstr << "#include \"" << pov_sstr.str() << "\"" << std::endl;
			}
			pov_file_ofstr.close();
		}


		pov_file_ofstr.open (pov_filec);

	  	for(int atom=0; atom<num_atoms; atom++){
	
			double red,green,blue,ireal;
			ireal = atoms::z_spin_array[atom];

			if(ireal>0.8){
				red = 0.0;
				green = 0.0;
				blue = 1.0;
			}
			else if(ireal>=0.0){
				red = 1.0-ireal*1.2;
				green = 1.0-ireal*1.2;
				blue = 1.0;
			}
			else if(ireal>=-0.8){
				red = 1.0;
				green = 1.0+ireal*1.2;
				blue = 1.0+ireal*1.2;
			}
			else if(ireal<-0.8){
				red = 1.0;
				green = 0.0;
				blue = 0.0;
			}
			else{
				red = 1.0;
				green = 1.0;
				blue = 1.0;
			}

			if(blue<0.0) blue=0.0;
			if(red<0.0) red=0.0;
			if(green<0.0) green=0.0;

			//#ifdef MPICF
				//double cx=mpi_create_variables::mpi_atom_global_coord_array[3*atom+0]*material_parameters::lattice_space_conversion[0];
				//double cy=mpi_create_variables::mpi_atom_global_coord_array[3*atom+1]*material_parameters::lattice_space_conversion[1];
				//double cz=mpi_create_variables::mpi_atom_global_coord_array[3*atom+2]*material_parameters::lattice_space_conversion[2];
			//#else
				double cx=atoms::x_coord_array[atom];  //*material_parameters::lattice_space_conversion[0];
				double cy=atoms::y_coord_array[atom];  //*material_parameters::lattice_space_conversion[1];
				double cz=atoms::z_coord_array[atom];  //*material_parameters::lattice_space_conversion[2];
			//#endif
			double sx=0.5*atoms::x_spin_array[atom];
			double sy=0.5*atoms::y_spin_array[atom];
			double sz=0.5*atoms::z_spin_array[atom];


	  		pov_file_ofstr << "cone {<" << cx << "+" << sx << "*sscale,"
										 << cy << "+" << sy << "*sscale,"
										 << cz << "+" << sz << "*sscale>,0.0 <"
										 << cx << "-" << sx << "*sscale,"
										 << cy << "-" << sy << "*sscale,"
										 << cz << "-" << sz << "*sscale>,sscale*0.5 "
						<< "texture { pigment {color rgb <" << red << " " << green << " " << blue << ">}"
						<< "finish {reflection {ref} diffuse 1 ambient 0}}}" << std::endl;
	  	}
	
		pov_file_ofstr.close();

		counter++;

	return 0;
	}

}
/*
if(1==0){
		ofstream spin_file;
		spin_file.open ("spins.dat");
		spin_file << atoms::num_atoms << endl;
		spin_file << "" << endl;
	  	
	  	for(int atom=0; atom<atoms::num_atoms; atom++){
	  		spin_file << material_parameters::material[atoms::type_array[atom]].element << "\t" << 
	  					atoms::x_coord_array[atom]*material_parameters::lattice_space_conversion[0] << "\t" << 
	  					atoms::y_coord_array[atom]*material_parameters::lattice_space_conversion[1] << "\t" << 
	  					atoms::z_coord_array[atom]*material_parameters::lattice_space_conversion[2] << "\t" <<
	  					atoms::x_spin_array[atom] << "\t" << 
	  					atoms::y_spin_array[atom] << "\t" << 
	  					atoms::z_spin_array[atom] << "\t" << endl;
	  	}
	
		spin_file.close();
	}
	*/
	
