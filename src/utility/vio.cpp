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
#include "voronoi.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "units.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

//==========================================================
// Global Output Streams
//==========================================================
//std::ofstream vinfo("info");
//std::ofstream vdp("vdp");
//std::ofstream vmag("vmag");

/// @namespace
/// @brief Contains variables and functions for reading in program data.
/// 
/// @internal
///=====================================================================================
///
namespace vin{
	
// Function Prototypes
//int read(string const);
int match(string const, string const, string const, string const, int const);
int read_mat_file(std::string const);
int match_create(std::string const, std::string const, int const);
int match_dimension(std::string const, std::string const, std::string const, int const);
int match_sim(std::string const, std::string const, std::string const, int const);
int match_material(string const, string const, string const, int const, int const, int const);

/// @brief Function to read in variables from a file.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.1
/// @date    18/01/2010
///
/// @param[in] filename Name of file to be opened 
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		14/01/2010
///	Revision:	  ---
///=====================================================================================
///
int read(string const filename){
	// ifstream declaration
	std::ifstream inputfile;
	
	// Open file read only
	inputfile.open(filename.c_str());
	
	// Check for opening
	if(!inputfile.is_open()){
		return EXIT_FAILURE;   // return to calling function for error checking or message
	}
	int line_counter=0;
	// Loop over all lines and pass keyword to matching function
	while (! inputfile.eof() ){
		line_counter++;
		// read in whole line
		std::string line;
		getline(inputfile,line);

		// Clear whitespace and tabs
		line.erase(remove(line.begin(), line.end(), '\t'), line.end());
		line.erase(remove(line.begin(), line.end(), ' '), line.end());

		// strip key,word,unit,value
		std::string key="";
		std::string word="";
		std::string value="";
		std::string unit="";

		// get size of string
		int linelength = line.length();
		int last=0;
		
		// set character triggers
		const char* colon=":";	// Word identifier
		const char* eq="=";		// Value identifier
		const char* exc="!";		// Unit identifier
		const char* hash="#";	// Comment identifier
		//const char* arrow=">";	// List identifier
		
		// Determine key by looping over characters in line
		for(int i=0;i<linelength;i++){
			char c=line.at(i);
			last=i;
			
			// if character is not "." or "=" or "!" or "#" interpret as key
			if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
				key.push_back(c);
			}
			else break;
		}
		const int end_key=last;
		
		// Determine the rest
		for(int i=end_key;i<linelength;i++){
			
			char c=line.at(i);
			//last=i;
				// period found - interpret as word
				if(c== *colon){
					for(int j=i+1;j<linelength;j++){
						// if character is not special add to value
						char c=line.at(j);
						if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
							word.push_back(c);
						}
						// if character is special then go back to main loop
						else{
							i=j-1;
							break;
						}
					}
				}
				// equals found - interpret as value
				else if(c== *eq){
					for(int j=i+1;j<linelength;j++){
						// if character is not special add to value
						char c=line.at(j);
						if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
							value.push_back(c);
						}
						// if character is special then go back to main loop
						else{
							i=j-1;
							break;
						}
					}
				}
				// exclaimation mark found - interpret as unit
				else if(c== *exc){
					for(int j=i+1;j<linelength;j++){
						// if character is not special add to value
						char c=line.at(j);
						if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
							unit.push_back(c);
						}
						// if character is special then go back to main loop
						else{
							i=j-1;
							break;
						}
					}
				}
				// hash found - interpret as comment
				else if(c== *hash){
					break;
				}
				//break;
		}
		string empty="";
		if(key!=empty){
		//std::cout << "\t" << "key:  " << key << std::endl;
		//std::cout << "\t" << "word: " << word << std::endl;
		//std::cout << "\t" << "value:" << value << std::endl;
		//std::cout << "\t" << "unit: " << unit << std::endl;
		int matchcheck = match(key, word, value, unit, line_counter);
		if(matchcheck==EXIT_FAILURE){
			exit(1);
		}
		}
	}
	// Close file
	inputfile.close();

	return EXIT_SUCCESS;
}

/// @brief Function to match keywords, variables and units to an initialisation variable.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.1
/// @date    18/01/2010
///
/// @param[in] keyword Unique string variable linked to an initialisation variable 
/// @param[in] value Value of keyword linked to initialisation variable 
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		15/01/2010
///	Revision:	  ---
///=====================================================================================
///
int match(string const key, string const word, string const value, string const unit, int const line){
//int match(string const key, string const word, string const value, string const unit, int const line, std::ifstream& inputfile){

	std::string test;

	//===================================================================
	// Test for create variables
	//===================================================================
	test="create";
	if(key==test){
		int frs=vin::match_create(word, value, line);
		return frs;
	}
	//===================================================================
	// Test for dimension variables
	//===================================================================
	else
	test="dimensions";
	if(key==test){
		int frs=vin::match_dimension(word, value, unit, line);
		return frs;
	}
	//===================================================================
	// Test for simulation variables
	//===================================================================
	else
	test="sim";
	if(key==test){
		int frs=vin::match_sim(word, value, unit, line);
		return frs;
	}
	//-------------------------------------------------------------------
	// Get material filename
	//-------------------------------------------------------------------
	else
	test="material";
	if(key==test){
		test="file";
		if(word==test){
			std::string matfile=value;
			// strip quotes
			matfile.erase(remove(matfile.begin(), matfile.end(), '\"'), matfile.end());
			test="";
			if(matfile!=test){
				//std::cout << matfile << std::endl;
				read_mat_file(matfile);
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - empty filename in control statement \'material:" << word << "\' on line " << line << " of input file" << std::endl;
				return EXIT_FAILURE;
			}
		}
		else{
			std::cerr << "Error - Unknown control statement \'material:" << word << "\' on line " << line << " of input file" << std::endl;
			return EXIT_FAILURE;
		}
	}
	else
		std::cerr << "Error - Unknown control statement \'" << key <<":"<< word << "\' on line " << line << " of input file" << std::endl;
		return EXIT_FAILURE;

} // end of match function

int match_create(string const word, string const value, int const line){
		//-------------------------------------------------------------------
		// system_creation_flags[1] - Set system particle shape
		//-------------------------------------------------------------------
		std::string test="full";
		if(word==test){
			mp::system_creation_flags[1]=0;
			return EXIT_SUCCESS;
		}
		else 
		test="cube";
		if(word==test){
			mp::system_creation_flags[1]=1;
			return EXIT_SUCCESS;
		}
		else 
		test="cylinder";
		if(word==test){
			mp::system_creation_flags[1]=2;
			return EXIT_SUCCESS;
		}
		else
		test="ellipsinder";
		if(word==test){
			mp::system_creation_flags[1]=3;
			return EXIT_SUCCESS;
		}
		else
		test="sphere";
		if(word==test){
			mp::system_creation_flags[1]=4;
			return EXIT_SUCCESS;
		}
		else
		test="truncated-octahedron";
		if(word==test){
			mp::system_creation_flags[1]=5;
			return EXIT_SUCCESS;
		}
		else
		//-------------------------------------------------------------------
		// system_creation_flags[2] - Set system type
		//-------------------------------------------------------------------
		test="particle";
		if(word==test){
			mp::system_creation_flags[2]=0;
			return EXIT_SUCCESS;
		}
		else
		test="particle-array";
		if(word==test){
			mp::system_creation_flags[2]=1;
			return EXIT_SUCCESS;
		}
		else
		test="hex-particle-array";
		if(word==test){
			mp::system_creation_flags[2]=2;
			return EXIT_SUCCESS;
		}
		else
		test="voronoi-film";
		if(word==test){
			mp::system_creation_flags[2]=3;
			return EXIT_SUCCESS;
		}
		//--------------------------------------------------------------------
		else
		test="voronoi-variance";
		if(word==test){
			double vsd=atof(value.c_str());
			if((vsd<0.0) || (vsd > 1.0)){
				std::cerr << "Error in input file - voronoi variance is outside of valid range (0.0-1.0)" << std::endl;
				return EXIT_FAILURE;}
			else{
				create_voronoi::voronoi_sd=vsd;
				return EXIT_SUCCESS;
			}
		}
		//--------------------------------------------------------------------
		else
		test="voronoi-parity";
		if(word==test){
			int vp=atoi(value.c_str());
			if((vp !=0) && (vp !=1)){
				std::cerr << "Error in input file - voronoi parity must be either 0 or 1, you input "<< vp << std::endl;
				return EXIT_FAILURE;}
			else{
				create_voronoi::parity=vp;
				return EXIT_SUCCESS;
			}
		}
		else
		//-------------------------------------------------------------------
		// system_creation_flags[3] - Set neighbourlist type
		//-------------------------------------------------------------------
		test="Jij-explicit";
		if(word==test){
			mp::system_creation_flags[3]=0;
			return EXIT_SUCCESS;
		}
		else
		//-------------------------------------------------------------------
		// system_creation_flags[4] - Set Multilayer Flag
		//-------------------------------------------------------------------
		/*test="multilayer";
		if(word==test){
			mp::system_creation_flags[4]=1;
			// test for multilayer lines
			test=">";
			if(value==test){
				//read_multilayer();
				std::cout << "read data" << std::endl;
				return EXIT_SUCCESS;
			}
			// test for missing parameter
			else
			test="";
			if(value==test){
				std::cerr << "Error - multilayers require additional parameters and none are specified" << std::endl;
				return EXIT_FAILURE;
			}
			else{
				//read_multilayer(value);
				std::cout << "read data from file "<< value << std::endl;
				return EXIT_SUCCESS;
			}
			
		}
		else
		test="intermixed";
		if(word==test){
			mp::system_creation_flags[4]=2;
			return EXIT_SUCCESS;
		}
		else*/
		test="particle-parity";
		if(word==test){
			int pp=atoi(value.c_str());
				mp::particle_creation_parity=pp;
				//std::cout << "ax: " << mp::lattice_constant[0] << std::endl;
				return EXIT_SUCCESS;
		}
		//--------------------------------------------------------------------
		else
		test="crystal-structure";
		if(word==test){
			// Strip quotes
			std::string cs=value;
			cs.erase(remove(cs.begin(), cs.end(), '\"'), cs.end());
			mp::crystal_structure=cs;
			return EXIT_SUCCESS;
		}
		//--------------------------------------------------------------------
		else
		test="single-spin";
		if(word==test){
			mp::single_spin=true;
			return EXIT_SUCCESS;
		}
		//--------------------------------------------------------------------
		// keyword not found
		//--------------------------------------------------------------------
		else{
			std::cerr << "Error - Unknown control statement \'create:" << word << "\' on line " << line << " of input file" << std::endl;
			return EXIT_FAILURE;
		}
		
	return EXIT_SUCCESS;
}

int match_dimension(string const word, string const value, string const unit, int const line){
		//-------------------------------------------------------------------
		// System dimension variables
		//-------------------------------------------------------------------
		std::string test="a";
		if(word==test){
			double a=atof(value.c_str());
			string unit_type;
			units::convert(unit,a,unit_type);
			string str="length";
			if(unit_type==str){
				mp::lattice_constant[0]=a;
				mp::lattice_constant[1]=a;
				mp::lattice_constant[2]=a;
				//std::cout << "ax: " << mp::lattice_constant[0] << std::endl;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << std::endl;
				exit(1);
			}
		}
		else
		test="c";
		if(word==test){
			double c=atof(value.c_str());
			string unit_type;
			units::convert(unit,c,unit_type);
			string str="length";
			if(unit_type==str){
				mp::lattice_constant[2]=c;
				//std::cout << "ax: " << mp::lattice_constant[0] << std::endl;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << std::endl;
				exit(1);
			}
		}
		else
		//--------------------------------------------------------------------
		test="ax";
		if(word==test){
			double ax=atof(value.c_str());
			string unit_type;
			units::convert(unit,ax,unit_type);
			string str="length";
			if(unit_type==str){
				mp::lattice_constant[0]=ax;
				//std::cout << "ax: " << mp::lattice_constant[0] << std::endl;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << std::endl;
				exit(1);
			}
		}
		else
		//--------------------------------------------------------------------
		test="ay";
		if(word==test){
			double ay=atof(value.c_str());
			string unit_type;
			units::convert(unit,ay,unit_type);
			string str="length";
			if(unit_type==str){
				mp::lattice_constant[1]=ay;
				//std::cout << "ax: " << mp::lattice_constant[0] << std::endl;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << std::endl;
				exit(1);
			}
		}
		else
		//--------------------------------------------------------------------
		test="az";
		if(word==test){
			double az=atof(value.c_str());
			string unit_type;
			units::convert(unit,az,unit_type);
			string str="length";
			if(unit_type==str){
				mp::lattice_constant[2]=az;
				//std::cout << "az: " << mp::lattice_constant[0] << std::endl;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << std::endl;
				exit(1);
			}
		}
		else
		//--------------------------------------------------------------------
		test="d";
		if(word==test){
			double d=atof(value.c_str());
			string unit_type;
			units::convert(unit,d,unit_type);
			string str="length";
			if(unit_type==str){
				mp::system_dimensions[0]=d;
				mp::system_dimensions[1]=d;
				mp::system_dimensions[2]=d;
				//std::cout << "ax: " << mp::lattice_constant[0] << std::endl;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << std::endl;
				exit(1);
			}
		}
		else
		//--------------------------------------------------------------------
		test="dx";
		if(word==test){
			double dx=atof(value.c_str());
			string unit_type;
			units::convert(unit,dx,unit_type);
			string str="length";
			if(unit_type==str){
				mp::system_dimensions[0]=dx;
				//std::cout << "ax: " << mp::lattice_constant[0] << std::endl;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << std::endl;
				exit(1);
			}
		}
		else
		//--------------------------------------------------------------------
		test="dy";
		if(word==test){
			double dy=atof(value.c_str());
			string unit_type;
			units::convert(unit,dy,unit_type);
			string str="length";
			if(unit_type==str){
				mp::system_dimensions[1]=dy;
				//std::cout << "ax: " << mp::lattice_constant[0] << std::endl;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << std::endl;
				exit(1);
			}
		}
		else
		//--------------------------------------------------------------------
		test="dz";
		if(word==test){
			double dz=atof(value.c_str());
			string unit_type;
			units::convert(unit,dz,unit_type);
			string str="length";
			if(unit_type==str){
				mp::system_dimensions[2]=dz;
				//std::cout << "ax: " << mp::lattice_constant[0] << std::endl;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << std::endl;
				exit(1);
			}
		}
		else
		//--------------------------------------------------------------------
		test="particle-size";
		if(word==test){
			double psize=atof(value.c_str());
			string unit_type;
			units::convert(unit,psize,unit_type);
			string str="length";
			if(unit_type==str){
				mp::particle_scale=psize;
				//std::cout << "particle_size: " << mp::particle_scale << std::endl;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << std::endl;
				exit(1);
			}
		}
		else
		//--------------------------------------------------------------------
		test="particle-spacing";
		if(word==test){
			double pspacing=atof(value.c_str());
			string unit_type;
			units::convert(unit,pspacing,unit_type);
			string str="length";
			if(unit_type==str){
				mp::particle_spacing=pspacing;
				//std::cout << "particle_spacing: " << mp::particle_spacing << std::endl;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << std::endl;
				exit(1);
			}
		}
		//else
		//--------------------------------------------------------------------
		else{
			std::cerr << "Error - Unknown control statement \'dimension:"<< word << "\' on line " << line << " of input file" << std::endl;
			return EXIT_FAILURE;
		}
	
	
	return EXIT_SUCCESS;
}

int match_sim(string const word, string const value, string const unit, int const line){
		//-------------------------------------------------------------------
		// System simulation variables
		//-------------------------------------------------------------------
		std::string test="integrator";
		if(word==test){
			test="LLG-Heun";
			if(value==test){
				sim::integrator=0;
				return EXIT_SUCCESS;
			}
			test="Monte-Carlo";
			if(value==test){
				sim::integrator=1;
				return EXIT_SUCCESS;
			}
			test="LLG-Midpoint";
			if(value==test){
				sim::integrator=2;
				return EXIT_SUCCESS;
			}
			test="Constrained-Monte-Carlo";
			if(value==test){
				sim::integrator=3;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
				std::cerr << "\t\"LLG-Heun\"" << std::endl;
				std::cerr << "\t\"Monte-Carlo\"" << std::endl;
				err::vexit();
			}
		}
		//-------------------------------------------------------------------
		test="program";
		if(word==test){
			test="Benchmark";
			if(value==test){
				sim::program=0;
				return EXIT_SUCCESS;
			}
			test="Time-Series";
			if(value==test){
				sim::program=1;
				return EXIT_SUCCESS;
			}
			test="Hysteresis-Loop";
			if(value==test){
				sim::program=2;
				return EXIT_SUCCESS;
			}
			test="Static-Hysteresis-Loop";
			if(value==test){
				sim::program=3;
				return EXIT_SUCCESS;
			}
			test="Curie-Temperature";
			if(value==test){
				sim::program=4;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
				std::cerr << "\t\"Benchmark\"" << std::endl;
				std::cerr << "\t\"Time-Series\"" << std::endl;
				std::cerr << "\t\"Hysteresis-Loop\"" << std::endl;
				std::cerr << "\t\"Static-Hysteresis-Loop\"" << std::endl;
				std::cerr << "\t\"Curie-Temperature\"" << std::endl;
				err::vexit();
			}
		}
		//-------------------------------------------------------------------
		test="exchange";
		if(word==test){
			test="true";
			if(value==test){
				sim::hamiltonian_simulation_flags[0]=1;
				return EXIT_SUCCESS;
			}
			test="false";
			if(value==test){
				sim::hamiltonian_simulation_flags[0]=0;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - value for \'sim:" << word << "\' must be either \"true\" or \"false\"" << std::endl;
				exit(1);
			}
		}
		//-------------------------------------------------------------------
		test="anisotropy";
		if(word==test){
			test="true";
			if(value==test){
				sim::hamiltonian_simulation_flags[1]=1;
				return EXIT_SUCCESS;
			}
			test="false";
			if(value==test){
				sim::hamiltonian_simulation_flags[1]=0;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - value for \'sim:" << word << "\' must be either \"true\" or \"false\"" << std::endl;
				exit(1);
			}
		}
		//-------------------------------------------------------------------
		test="applied";
		if(word==test){
			test="true";
			if(value==test){
				sim::hamiltonian_simulation_flags[2]=1;
				return EXIT_SUCCESS;
			}
			test="false";
			if(value==test){
				sim::hamiltonian_simulation_flags[2]=0;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - value for \'sim:" << word << "\' must be either \"true\" or \"false\"" << std::endl;
				exit(1);
			}
		}
		//-------------------------------------------------------------------
		test="thermal";
		if(word==test){
			test="true";
			if(value==test){
				sim::hamiltonian_simulation_flags[3]=1;
				return EXIT_SUCCESS;
			}
			test="false";
			if(value==test){
				sim::hamiltonian_simulation_flags[3]=0;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - value for \'sim:" << word << "\' must be either \"true\" or \"false\"" << std::endl;
				exit(1);
			}
		}
		//-------------------------------------------------------------------
		test="dipolar";
		if(word==test){
			test="true";
			if(value==test){
				sim::hamiltonian_simulation_flags[4]=1;
				return EXIT_SUCCESS;
			}
			test="false";
			if(value==test){
				sim::hamiltonian_simulation_flags[4]=0;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - value for \'sim:" << word << "\' must be either \"true\" or \"false\"" << std::endl;
				exit(1);
			}
		}
		//--------------------------------------------------------------------
		test="dt";
		if(word==test){
			double dt=atof(value.c_str());
			// Test for valid range
			if((dt>=1.0E-20) && (dt<1.0E-6)){
				mp::dt_SI=dt;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - sim:" << word << " on line " << line << " of input file must be in the range 1E-20 - 1.0E-6" << std::endl;
				err::vexit();
			}
		}
		//--------------------------------------------------------------------
		test="total-time";
		if(word==test){
			int tt=atoi(value.c_str());
			sim::total_time=tt;
			return EXIT_SUCCESS;
		}
		//--------------------------------------------------------------------
		test="loop-time";
		if(word==test){
			int tt=atoi(value.c_str());
			sim::loop_time=tt;
			return EXIT_SUCCESS;
		}
		//--------------------------------------------------------------------
		test="partial-time";
		if(word==test){
			int tt=atoi(value.c_str());
			sim::partial_time=tt;
			return EXIT_SUCCESS;
		}
		//--------------------------------------------------------------------
		test="equilibration-time";
		if(word==test){
			int tt=atoi(value.c_str());
			sim::equilibration_time=tt;
			return EXIT_SUCCESS;
		}
		//--------------------------------------------------------------------
		test="maximum-temperature";
		if(word==test){
			double T=atof(value.c_str());
			// Test for valid range
			if((T>=0.0) && (T<1.0E10)){
				sim::Tmax=T;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - sim:" << word << " on line " << line << " of input file must be in the range 0 - 1.0E10" << std::endl;
				exit(1);
			}
		}
		//--------------------------------------------------------------------
		test="minimum-temperature";
		if(word==test){
			double T=atof(value.c_str());
			// Test for valid range
			if((T>=0.0) && (T<1.0E10)){
				sim::Tmin=T;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - sim:" << word << " on line " << line << " of input file must be in the range 0 - 1.0E10" << std::endl;
				exit(1);
			}
		}
		//--------------------------------------------------------------------
		test="temperature";
		if(word==test){
			double T=atof(value.c_str());
			// Test for valid range
			if((T>=0.0) && (T<1.0E10)){
				sim::temperature=T;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - sim:" << word << " on line " << line << " of input file must be in the range 0 - 1.0E10" << std::endl;
				exit(1);
			}
		}
		//--------------------------------------------------------------------
		test="delta-temperature";
		if(word==test){
			double T=atof(value.c_str());
			// Test for valid range
			if((T>=0.0) && (T<1.0E10)){
				sim::delta_temperature=T;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - sim:" << word << " on line " << line << " of input file must be in the range 0 - 1.0E10" << std::endl;
				exit(1);
			}
		}
		//--------------------------------------------------------------------
		test="H-applied";
		if(word==test){
			double H=atof(value.c_str());
			string unit_type="field";
			// if no unit given, assume internal
			if(unit.size() != 0){
				units::convert(unit,H,unit_type);
			}
			string str="field";
			if(unit_type==str){
				sim::H_applied=H;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'sim:" << word << "\'"<< std::endl;
				exit(1);
			}
		}
		//--------------------------------------------------------------------
		test="Hmin";
		if(word==test){
			double H=atof(value.c_str());
			string unit_type="field";
			// if no unit given, assume internal
			if(unit.size() != 0){
				units::convert(unit,H,unit_type);
			}
			string str="field";
			if(unit_type==str){
				sim::Hmin=H;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'sim:" << word << "\'"<< std::endl;
				exit(1);
			}
		}
		//--------------------------------------------------------------------
		test="Hmax";
		if(word==test){
			double H=atof(value.c_str());
			string unit_type="field";
			// if no unit given, assume internal
			if(unit.size() != 0){
				units::convert(unit,H,unit_type);
			}
			string str="field";
			if(unit_type==str){
				sim::Hmax=H;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'sim:" << word << "\'"<< std::endl;
				exit(1);
			}
		}
		//--------------------------------------------------------------------
		test="Hinc";
		if(word==test){
			double H=atof(value.c_str());
			string unit_type="field";
			// if no unit given, assume internal
			if(unit.size() != 0){
				units::convert(unit,H,unit_type);
			}
			string str="field";
			if(unit_type==str){
				sim::Hinc=H;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'sim:" << word << "\'"<< std::endl;
				exit(1);
			}
		}
	//double H_vec[3]={0.0,0.0,1.0};
		//--------------------------------------------------------------------
		else{
			std::cerr << "Error - Unknown control statement \'sim:"<< word << "\' on line " << line << " of input file" << std::endl;
			return EXIT_FAILURE;
		}
	
	
	return EXIT_SUCCESS;
}

// temporary array of materials for reading in material data
//std::cout << "here" << std::endl;
  std::vector<mp::materials_t> read_material(0);

int read_mat_file(std::string const matfile){
	
	// Declare input stream
	std::ifstream inputfile;
	
	// resize temporary materials array for storage of variables
	read_material.resize(mp::max_materials);

	// Open file read only
	inputfile.open(matfile.c_str());
	
	// Check for opening
	if(!inputfile.is_open()){
		std::cerr << "Error opening file " << matfile << "- file does not exist!" << std::endl; 
		exit(1);   // return to calling function for error checking or message
	}
	//-------------------------------------------------------
	// Material 0
	//-------------------------------------------------------


	
	int line_counter=0;
	// Loop over all lines and pass keyword to matching function
	while (! inputfile.eof() ){
		line_counter++;
		// read in whole line
		std::string line;
		getline(inputfile,line);

		// Clear whitespace, quotes and tabs
		line.erase(remove(line.begin(), line.end(), '\t'), line.end());
		line.erase(remove(line.begin(), line.end(), ' '), line.end());
		line.erase(remove(line.begin(), line.end(), '\"'), line.end());

		// strip key,word,unit,value
		std::string key="";
		std::string word="";
		std::string value="";
		std::string unit="";
		std::string index="";
		int super_index=0;
		int sub_index=0;

		// get size of string
		int linelength = line.length();
		int last=0;
		
		// set character triggers
		const char* colon=":";	// Word identifier
		const char* eq="=";		// Value identifier
		const char* exc="!";		// Unit identifier
		const char* hash="#";	// Comment identifier
		const char* si="[";		// Index identifier
		const char* ei="]";		// End index identifier
		
		// Determine key and super index by looping over characters in line
		for(int i=0;i<linelength;i++){
			char c=line.at(i);
			last=i;
			
			// if character is not ":" or "=" or "!" or "#" interpret as key
			if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash) && (c != *si) && (c != *ei)){
				key.push_back(c);
			}
			// Check for number of materials statement
			else if(c == *eq){
				// break to read in value
				break;
			}
			// Check for superindex
			else if(c ==*si){
				const int old=last;
				// Get super index
				for(int j=old+1;j<linelength;j++){
					c=line.at(j);
					if(c != *ei){
						index.push_back(c);
					}
					else{
						break;
					}
					last=j;
				}
				
				// check for valid index
				super_index = atoi(index.c_str());
				if((super_index>=0) && (super_index<mp::max_materials)){
					break;
				}
				else{
					std::cerr << "Invalid index number " << index << " on line " << line_counter << " in material input file" << std::endl;
					std::cerr << "Causes could be invalid character or outside of range, ie less than zero or greater than max_materials=" << mp::max_materials << ", exiting" << std::endl;
					exit(1);
				}
				
			}
			// For anything else
			else break;
		}
		const int end_key=last;
		
		//
		//exit(1);
		// Determine the rest
		for(int i=end_key;i<linelength;i++){
			
			char c=line.at(i);
			// colon found - interpret as word
			if(c== *colon){
				for(int j=i+1;j<linelength;j++){
					// if character is not special add to value
					char c=line.at(j);
					if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
						// check for sub-index
						if(c == *si){
							index="";
							while(line.at(j+1) != *ei){
								j++;
								index.push_back(line.at(j));
							}
							sub_index=atoi(index.c_str());
							// Check for valid index
							if((super_index<0) || (super_index>=mp::max_materials)){
								std::cerr << "Invalid sub-index number " << index << " on line " << line_counter << " in material input file" << std::endl;
								std::cerr << "Causes could be invalid character or outside of range, ie less than zero or greater than max_materials=" << mp::max_materials << ", exiting" << std::endl;
								exit(1);
							}
							// end of word
							break;
						}
						else word.push_back(c);
					}
					// if character is special then go back to main loop
					else{
						i=j-1;
						break;
					}
				}
			}
			// equals found - interpret as value
			else if(c== *eq){
				for(int j=i+1;j<linelength;j++){
					// if character is not special add to value
					char c=line.at(j);
					if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
						value.push_back(c);
					}
					// if character is special then go back to main loop
					else{
						i=j-1;
						break;
					}
				}
			}
			// exclaimation mark found - interpret as unit
			else if(c== *exc){
				for(int j=i+1;j<linelength;j++){
					// if character is not special add to value
					char c=line.at(j);
					if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
						unit.push_back(c);
					}
					// if character is special then go back to main loop
					else{
						i=j-1;
						break;
					}
				}
			}
			// hash found - interpret as comment
			else if(c== *hash){
				break;
			}
			//break;
		}
		string empty="";
		if(key!=empty){
			//std::cout << key << "[" << super_index << "]:" << word << "[" << sub_index << "]=" << value << " !" << unit << std::endl;
			//std::cout << "\t" << "key:  " << key << std::endl;
			//std::cout << "\t" << "word: " << word << std::endl;
			//std::cout << "\t" << "value:" << value << std::endl;
			//std::cout << "\t" << "unit: " << unit << std::endl;
			int matchcheck = vin::match_material(word, value, unit, line_counter, super_index, sub_index);
			if(matchcheck==EXIT_FAILURE){
				exit(1);
			}
		}
	}
	
	// resize global material array
	mp::material.resize(mp::num_materials);
	
	// Copy data to global material array
	for(int mat=0;mat<mp::num_materials;mat++){
		mp::material[mat]=read_material[mat];
	}
	
	// Resize read array to zero
	read_material.resize(0);
	
	
	// Close file
	inputfile.close();

	return EXIT_SUCCESS;
	
}

int match_material(string const word, string const value, string const unit, int const line, int const super_index, int const sub_index){
		//-------------------------------------------------------------------
		// system_creation_flags[1] - Set system particle shape
		//-------------------------------------------------------------------
		std::string test="num-materials";
		if(word==test){
			mp::num_materials=atoi(value.c_str());
			return EXIT_SUCCESS;
		}
		else 
		test="name";
		if(word==test){
			read_material[super_index].name=value;
			return EXIT_SUCCESS;
		}
		else 
		test="alpha";
		if(word==test){
			read_material[super_index].alpha=atof(value.c_str());
			return EXIT_SUCCESS;
		}
		else
		test="Jij_matrix";
		if(word==test){
			double Jij=atof(value.c_str());
			string unit_type="energy";
			// if no unit given, assume internal
			if(unit.size() != 0){
				units::convert(unit,Jij,unit_type);
			}
			string str="energy";
			if(unit_type==str){
				read_material[super_index].Jij_matrix_SI[sub_index]=Jij;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'"<< std::endl;
				exit(1);
			}
			//return EXIT_SUCCESS;
		}
		else
		test="mu_s";
		if(word==test){
			double mu_s=atof(value.c_str());
			string unit_type="moment";
			// if no unit given, assume internal
			if(unit.size() != 0){
				units::convert(unit,mu_s,unit_type);
			}
			string str="moment";
			if(unit_type==str){
				// Set moment flag
				read_material[super_index].moment_flag=true;
				read_material[super_index].mu_s_SI=mu_s;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'"<< std::endl;
				exit(1);
			}
		}
		else
		test="magnetisation";
		if(word==test){
			double mag=atof(value.c_str());
			string unit_type="magnetisation";
			// if no unit given, assume internal
			if(unit.size() != 0){
				units::convert(unit,mag,unit_type);
			}
			string str="magnetisation";
			if(unit_type==str){
				// Set moment flag
				read_material[super_index].moment_flag=false;
				read_material[super_index].magnetisation=mag;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'"<< std::endl;
				exit(1);
			}
		}
		else
		test="Ku1";
		if(word==test){
			double K=atof(value.c_str());
			string unit_type="anisotropy";
			// if no unit given, assume internal
			if(unit.size() != 0){
				units::convert(unit,K,unit_type);
				//read_material[super_index].anis_flag=false;
				//std::cout << "setting flag to false" << std::endl;
			}
			string str="anisotropy";
			if(unit_type==str){
				// Set moment flag
				read_material[super_index].Ku1_SI=K;
				return EXIT_SUCCESS;
			}
			else{
				std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'"<< std::endl;
				exit(1);
			}
		}
		else
		test="gamma-rel";
		if(word==test){
			read_material[super_index].gamma_rel=atof(value.c_str());
			return EXIT_SUCCESS;
		}
		//------------------------------------------------------------
		else
		test="Sx";
		if(word==test){
			read_material[super_index].initial_spin[0]=atof(value.c_str());
			read_material[super_index].random_spins=false;
			return EXIT_SUCCESS;
		}
		//------------------------------------------------------------
		else
		test="Sy";
		if(word==test){
			read_material[super_index].initial_spin[1]=atof(value.c_str());
			read_material[super_index].random_spins=false;
			return EXIT_SUCCESS;
		}
		//------------------------------------------------------------
		else
		test="Sz";
		if(word==test){
			read_material[super_index].initial_spin[2]=atof(value.c_str());
			read_material[super_index].random_spins=false;
			return EXIT_SUCCESS;
		}
		//------------------------------------------------------------
		else
		test="random-spins";
		if(word==test){
			string t="true";
			string f="false";
			if(value==t){
				read_material[super_index].random_spins=true;
				return EXIT_SUCCESS;
			}
			else if(value==f){
				read_material[super_index].random_spins=false;
				return EXIT_SUCCESS;
			}
			else {
				std::cerr << "Error in input file - material[" << super_index << "]:random-spins must be either true or false" << std::endl;
				return EXIT_FAILURE;
			}
		}
		//------------------------------------------------------------
		else
		test="hamiltonian";
		if(word==test){
			read_material[super_index].hamiltonian_type=value;
			return EXIT_SUCCESS;
		}
		else
		test="element";
		if(word==test){
			// Test for 3 characters
			if(value.length()>3){
				std::cerr << "Error - element identifier on line "<<  line << " of material file must be a maximum of three characters long" << std::endl;
			}
			else{
				// pad value to be equal to 3 characters
				string tmp="   ";
				for(unsigned int i=0;i<3;i++){
					if(i<value.length()){
						tmp.at(i)=value.at(i);
					}
				}
				read_material[super_index].element=tmp;
				return EXIT_SUCCESS;
			}
		}
		//--------------------------------------------------------------------
		else
		test="crystal-structure";
		if(word==test){
			read_material[super_index].crystal_structure=value;
			return EXIT_SUCCESS;
		}
		else
		test="geometry";
		if(word==test){
			// Open geometry file
			std::ifstream gfile(value.c_str());
			if(!gfile.is_open()){
				std::cerr << "Error - geometry file " << value.c_str() << " not found, exiting!" << std::endl;
				return EXIT_FAILURE;
			}
			gfile >> read_material[super_index].geometry;
			if((read_material[super_index].geometry<3) || (read_material[super_index].geometry>100)){
				std::cerr << "Error in geometry input file " << value.c_str() << " - first number must be non zero integer in the range 3-100"<< std::endl;
				return EXIT_FAILURE;
			}
			//std::cout << "ngp " << read_material[super_index].geometry << std::endl;
			for(int c=0;c<read_material[super_index].geometry;c++){
				for(int xy=0;xy<2;xy++){
					double var;
					gfile >> var;
					if(gfile.eof()){
						std::cerr << "Error in geometry input file " << value.c_str() << " end of file reached before reading all coordinates" << std::endl;
						return EXIT_FAILURE;
					}
					read_material[super_index].geometry_coords[c][xy];
					if((var<0.0) || (var > 1.0)){
						std::cerr << "Error in geometry input file " << value.c_str() << " value is outside of valid range (0.0-1.0)" << std::endl;
						return EXIT_FAILURE;
					}
					else read_material[super_index].geometry_coords[c][xy]=var;
				}
				//std::cout << read_material[super_index].geometry_coords[c][0] << "\t" << read_material[super_index].geometry_coords[c][1] << std::endl;
			}
			//double min=atof(value.c_str());
			//if((min<-0.11) || (min > 1.11)){
			//	std::cerr << "Error in input file - material[" << super_index << "]:min is outside of valid range (0.0-1.0)" << std::endl;
			//	return EXIT_FAILURE;}
			//else{
			//	read_material[super_index].min=min;
				return EXIT_SUCCESS;
			//}
		}
		//--------------------------------------------------------------------
		else
		test="alloy-master";
		if(word==test){
			read_material[super_index].alloy_master=true;
			return EXIT_SUCCESS;
		}
		//--------------------------------------------------------------------
		else
		test="alloy-class";
		if(word==test){
			int ac=atoi(value.c_str());
			if((ac<0) || (ac > 0)){
				std::cerr << "Error in input file - material[" << super_index << "]:alloy-class is outside of valid range (0-0)" << std::endl;
				return EXIT_FAILURE;
			}
			else{
				read_material[super_index].alloy_class=ac;
				return EXIT_SUCCESS;
			}
		}
		//--------------------------------------------------------------------
		else
		test="alloy";
		if(word==test){
			double a=atof(value.c_str());
			if((a < 0.0) || (a > 1.0)){
				std::cerr << "Error in input file - material[" << super_index << "]:alloy["<< sub_index << "] is outside of valid range (0.0-1.0)" << std::endl;
				return EXIT_FAILURE;
			}
			else{
				read_material[super_index].alloy[sub_index]=a;
				return EXIT_SUCCESS;
			}
			//return EXIT_SUCCESS;
		}
		//--------------------------------------------------------------------
		else
		test="min";
		if(word==test){
			double min=atof(value.c_str());
			if((min<0.0) || (min > 1.0)){
				std::cerr << "Error in input file - material[" << super_index << "]:min is outside of valid range (0.0-1.0)" << std::endl;
				return EXIT_FAILURE;}
			else{
				read_material[super_index].min=min;
				return EXIT_SUCCESS;
			}
		}
		//--------------------------------------------------------------------
		else
		test="max";
		if(word==test){
			double max=atof(value.c_str());
			if((max<0.0) || (max > 1.0)){
				std::cerr << "Error in input file - material[" << super_index << "]:max is outside of valid range (0.0-1.0)" << std::endl;
				return EXIT_FAILURE;}
			else{
				read_material[super_index].max=max;
				return EXIT_SUCCESS;
			}
		}
		else
		//--------------------------------------------------------------------
		test="core-shell-size";
		if(word==test){
			double css=atof(value.c_str());
			if((css<0.0) || (css > 1.0)){
				std::cerr << "Error in input file - material[" << super_index << "]:core-shell-size is outside of valid range (0.0-1.0)" << std::endl;
				return EXIT_FAILURE;}
			else{
				read_material[super_index].core_shell_size=css;
				return EXIT_SUCCESS;
			}
		}
		//-------------------------------------------------------------------
		else
		test="interface-roughness";
		if(word==test){
			double ir=atof(value.c_str());
			if((ir<0.0) || (ir > 1.0)){
				std::cerr << "Error in input file - material[" << super_index << "]:interface-roughness is outside of valid range (0.0-1.0)" << std::endl;
				return EXIT_FAILURE;}
			else{
				read_material[super_index].interface_roughness=ir;
				return EXIT_SUCCESS;
			}
		}
		else
		//-------------------------------------------------------------------
		test="density";
		if(word==test){
			double d=atof(value.c_str());
			if((d<0.0) || (d > 1.0)){
				std::cerr << "Error in input file - material[" << super_index << "]:density is outside of valid range (0.0-1.0)" << std::endl;
				return EXIT_FAILURE;}
			else{
				read_material[super_index].density=d;
				return EXIT_SUCCESS;
			}
		}
		else
		test="continuous";
		if(word==test){
			string t="true";
			string f="false";
			if(value==t){
				read_material[super_index].continuous=true;
				return EXIT_SUCCESS;
			}
			else if(value==f){
				read_material[super_index].continuous=false;
				return EXIT_SUCCESS;
			}
			else {
				std::cerr << "Error in input file - material[" << super_index << "]:continuous must be either true or false" << std::endl;
				return EXIT_FAILURE;
			}
		}
		else
		//-------------------------------------------------------------------
		test="intermixing";
		if(word==test){
			double i=atof(value.c_str());
			if((i<0.0) || (i > 1.0)){
				std::cerr << "Error in input file - material[" << super_index << "]:intermixing[" << sub_index <<"] is outside of valid range (0.0-1.0)" << std::endl;
				return EXIT_FAILURE;}
			else{
				read_material[super_index].intermixing[sub_index]=i;
				return EXIT_SUCCESS;
			}
		}
		//--------------------------------------------------------------------

		//--------------------------------------------------------------------
		// keyword not found
		//--------------------------------------------------------------------
		else{
			std::cerr << "Error - Unknown control statement \'material[" << super_index << "]:" << word << "\' on line " << line << " of material file" << std::endl;
			return EXIT_FAILURE;
		}
		
	return EXIT_SUCCESS;
}



} // end of namespace vin 


namespace vout{
  std::ofstream errfile;
  null_streambuf nullbuf;

  void redirect(std::ostream& strm, std::string filename) {
    errfile.open(filename.c_str());
    // redirect ouput into the file
    strm.rdbuf (errfile.rdbuf());
  }

  void nullify(std::ostream& strm){
    strm.rdbuf(&nullbuf);
  }
  
} // end of namespace vout  */
//==========================================================
// Namespace output
//==========================================================

namespace pov{
	int counter=0;
}


int output_pov_file(){
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
		pov_file_sstr << ".pov";
		//pov_file_sstr << "spins." << mpi_generic::my_rank << "." << counter << ".pov";
		std::string pov_file = pov_file_sstr.str();
		const char* pov_filec = pov_file.c_str();
		std::ofstream pov_file_ofstr;

		if(vmpi::my_rank==0){
			std::stringstream pov_hdr_sstr;
			pov_hdr_sstr << "spins." << counter << ".pov";
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
			//pov_file_ofstr << "#declare CX=" << material_parameters::system_dimensions[0]*6.0 << ";" << std::endl;
			//pov_file_ofstr << "#declare CY=" << material_parameters::system_dimensions[1]*6.0 << ";" << std::endl;
			//pov_file_ofstr << "#declare CZ=" << material_parameters::system_dimensions[2]*6.0 << ";" << std::endl;
	 		pov_file_ofstr << "#declare ref=0.4;" << std::endl;
	 		pov_file_ofstr << "#declare sscale=2.0;" << std::endl;
			pov_file_ofstr << "background { color Gray30 }" << std::endl;
			pov_file_ofstr << "Set_Camera(<CX,CY,CZ>, <LX,LY,LZ>, 15)" << std::endl;
			pov_file_ofstr << "Set_Camera_Aspect(4,3)" << std::endl;
			pov_file_ofstr << "Set_Camera_Sky(<0,0,1>)" << std::endl;
	   	pov_file_ofstr << "light_source { <2*CX, 2*CY, 2*CZ> color White}" << std::endl;
			for(int p =0;p<vmpi::num_processors;p++){
				std::stringstream pov_sstr;
				pov_sstr << "spins." << std::setfill('0') << std::setw(3) << p << "." << std::setfill('0') << std::setw(5) << counter << ".pov";
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
			//	double cx=mpi_create_variables::mpi_atom_global_coord_array[3*atom+0]*material_parameters::lattice_space_conversion[0];
			//	double cy=mpi_create_variables::mpi_atom_global_coord_array[3*atom+1]*material_parameters::lattice_space_conversion[1];
			//	double cz=mpi_create_variables::mpi_atom_global_coord_array[3*atom+2]*material_parameters::lattice_space_conversion[2];
			//#else
				double cx=atoms::x_coord_array[atom]; //*material_parameters::lattice_space_conversion[0];
				double cy=atoms::y_coord_array[atom]; //*material_parameters::lattice_space_conversion[1];
				double cz=atoms::z_coord_array[atom]; //*material_parameters::lattice_space_conversion[2];
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
	
