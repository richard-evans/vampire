/// \file main.cpp
/// Vampire Program Entry

#include <iostream>
#include <vector>
#include <stdio.h>

#include "create.hpp"
#include "material.hpp"
#include "public.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
//#include <types.h>
//#include <wait.h>

// comment
// another comment
using namespace std;

//int initialise_variables(std::string const);
//int create_system();
int simulate_system();
//int handle_opengl(int, char**);
#ifdef MPICF
int initialise_mpi();
int finalise_mpi();
int mpi_hosts();
#endif 

//int initialise_system();

//int main (int argc, char** argv)
/// Main function for vampire
/// Prints out program header and calls main program routines
int main(int argc, char* argv[]){

	//=============================================================
	// Check for valid command-line arguments
	//=============================================================
	std::string infile="vinput";
	
	for(int arg = 1; arg < argc; arg++){
		std::string sw=argv[arg];
		// input file
		if(sw=="-f"){
			// check number of args not exceeded
			if(arg+1 < argc){
				arg++;
				infile=string(argv[arg]);
			}
			else{
				std::cerr << "Error - no file specified for \'-f\' command line option" << std::endl;
				return EXIT_FAILURE;
			}
		}
		else{
			std::cerr << "Error - unknown command line parameter \'" << sw << "\'" << std::endl;
			return EXIT_FAILURE;
		}
	}
	
	//=============================================================
	// For parallel execution intialise MPI
	//=============================================================	
	
	#ifdef MPICF
	initialise_mpi();
	#endif 
	
   //=============================================================
	//      Output Program Header
	//=============================================================
	if(vmpi::my_rank==0){
		cout << "                   __     __                    _ " << endl;          
		cout << "                   \\ \\   / /_ _ _ __ ___  _ __ (_)_ __ ___ " << endl;
		cout << "                    \\ \\ / / _` | '_ ` _ \\| '_ \\| | '__/ _ \\ " << endl;
		cout << "                     \\ V / (_| | | | | | | |_) | | | |  __/ " << endl;
		cout << "                      \\_/ \\__,_|_| |_| |_| .__/|_|_|  \\___| " << endl;
		cout << "                                         |_|                " << endl;
		cout << " " << endl;
		cout << "          Visual Atomistic and Micromagnetic Parallel IntegratoR Engine" << endl;
		cout << " " << endl;
		cout << "           Contributors: Richard F L Evans, Joe Barker, Thomas Ostler" << endl;
		cout << "                         Weijia Fan, Roy W Chantrell" << endl;
		cout << " " << endl;
		cout << "                      Version 1.0 " << __DATE__ << " " << __TIME__ << endl;
		#ifdef COMP	
		cout << "                        Compiled with " << COMP << endl;
		#endif 
		cout << " " << endl;
		cout << "================================================================================" << endl;
		int sysstat = system ("date");
		if(sysstat!=0) std::cerr << "Error retrieving date from system" << std::endl;
	}

	#ifdef MPICF
	mpi_hosts();
	#endif 

	//if(mpi_generic::my_rank==0){
	//	cout << "================================================================================" << endl;
	//	cout << " " << endl;
	//	cout << "Starting Simulation..." << endl;
	//}
    //=============================================================
	//      Initialise material parameters and atomistic variables
	//=============================================================

	if(vmpi::my_rank==0){
		cout << "================================================================================" << endl;
		cout << " " << endl;
		cout << "Initialising system variables" << endl;
	}
	//cout << "Initialising system variables"<< endl;
	mp::initialise(infile);

	//=============================================================
	//      Create atomistic system, neighbourlist etc
	//=============================================================
     
	if(vmpi::my_rank==0){
		cout << "Creating system" << endl;
	}
	//  cout << "Creating atomic system"<< endl;
	cs::create();
	//create_system();

	//initialise_mpi();

	//cout << mpi_generic::my_rank << endl;

		//initialise_system();
        
		//for (int time=0;time<10001;time++){
	if(vmpi::my_rank==0){
		cout << "Starting Simulation..." << endl;
	}
		//for(;;){
	simulate_system();
		//}
     
	//=============================================================
	//      Deallocate allocated arrays close output files etc
	//=============================================================

	//cout << "Cleaning Up Atomistic System"<< endl;
	//  cleanup_system();

	//=============================================================
	//      Finalise MPI
	//=============================================================

	#ifdef MPICF
	finalise_mpi();
	#endif

  return 0;
}

/// \mainpage Vampire
/// 
/// \section intro_sec Introduction
/// 
/// Vampire is an anacronym for Visual Atomistic and Micromagnetic Parallel IntegratoR Engine and 
/// simulates the magnetic properties of materials using a classical spin model. Models can be 
/// defined from the atomistic scale up to micrometer scale using a variety of methods.The code is
/// open source and has been developed by the Computational Magnetism Group at The University of
/// York, United Kingdom.
///
/// \section features Program Features
/// \subsection features_ss1 System Generation
/// \arg Common crystal structures including Simple Cubic, FCC, BCC, Hexagonal close packed and user-defined
/// \arg Single Nanoparticles with spherical, truncated octahedron, cylindrical, and cubic geometries
/// \arg Regular arrays of nanoparticles, eg bit-patterned recording media
/// \arg Voronoi thin films generated by qvoronoi
/// \arg Multilayered materials including interface roughness and intermixing
/// \subsection features_ss2 Material Properties
/// \arg Generic Heisenberg Hamiltonians (Exchange, Applied Field)
/// \arg Uniaxial, Cubic, 2-ion and Neel Surface Anisotropies
/// \arg User-definable Hamiltonians (greater than nearest-neighbour)
/// \arg Support for Ab-initio input data (Exchange, Anisotropies)
/// \subsection features_ss3 Atomistic Spin Dynamics
/// \arg Landau-Lifshitz-Gilbert Equation of motion with Langevin Dynamics
/// \arg Landau-Lifshitz Equation of motion with Coloured Noise
/// \arg Heun and Semi-analytical numerical solvers
/// \subsection features_ss4 Monte Carlo Methods
/// \arg Basic Monte Carlo integration
/// \arg Constrained Monte Carlo
/// \subsection features_ss5 Energy Minimisation
/// \arg LaGrange Multiplier Energy Minimisation
/// \subsection features_ss6 Parallel Features
/// Vampire has been designed to run efficiently on parallel supercomputers 
/// and supports the following parallel modes:
/// \arg Parallel Statistics
/// \arg 2D-decomposition (Thin films)
/// \arg CUDA Acceleration
/// \subsection features_ss5 Visualisation
/// \arg Output to PovRAY Format
/// \arg Realtime OpenGL visualisation (Coming Soon!)
/// \section install_sec Installation
/// Vampire is distributed as both executable (serial version) and source code (serial, parallel, and CUDA versions). 
/// Compilation has been tested on a wide range of C++ compilers, including GNU,
/// Intel and Pathscale. 
/// \subsection exec Installation for binary distribution
/// The executables, libraries and scripts come packaged in a vampire.x.x.xxx.tar.gz file. To install the software, 
/// first unpack the archive with: \n
/// \verbatim tar -xzf vampire.x.x.xxx.tar.gz \endverbatim
/// Then change into the unpacked directory \n
/// \verbatim cd vampire \endverbatim
/// Finally run the install script as super user \n
/// \verbatim sudo ./install.sh \endverbatim
///
/// \subsection bin Installation for source code distribution
/// Essentially three bits of information are required to compile 
///  
/// etc...
/// 
