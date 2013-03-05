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
#include <iostream>
#include <vector>
#include <sstream>

#include "create.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

int simulate_system();

/// Main function for vampire
/// Prints out program header and calls main program routines
int main(int argc, char* argv[]){


	//=============================================================
	// Check for valid command-line arguments
	//=============================================================
	std::string infile="input";
	
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
      if(sw=="--help"){
         // Output mini-header
         std::cout << "vampire 2.0.1 http://vampire.york.ac.uk/" << std::endl;
         std::cout << "Copyright (c) 2009-2013 R F L Evans" << std::endl << std::endl;
         // check number of args not exceeded
         if(arg+1 < argc){
            arg++;
            std::string keyword=string(argv[arg]);
            vio::keyword_help(keyword);
            return EXIT_SUCCESS;
            exit(EXIT_SUCCESS);
         }
         else{
            std::cout << "\x1b[1m" << "NAME" << "\x1b[0m" << std::endl;
            std::cout << "        vampire - atomistic simulator for magnetic materials" << std::endl << std::endl;
            std::cout << "\x1b[1m" << "SYNOPSIS" << "\x1b[0m" << std::endl;
            std::cout << "        vampire [--version] [--help<keyword>] [-f<input_file>]" << std::endl << std::endl;
            std::cout << "\x1b[1m" << "DESCRIPTION" << "\x1b[0m" << std::endl;
            std::cout << "        Vampire is an open source software simulation package for atomistic \n"
                         "        simulation of magnetic materials developed at the University of York.\n"
                         "        The objective of vampire is to provide a community standard tool for \n"
                         "        atomistic simulation of magnetic materials with the highest possible\n"
                         "        performance." << std::endl << std::endl;
            std::cout << "        Using a variety of common simulation methods it can calculate the \n"
                         "        equilibrium and dynamic magnetic properties of a wide variety of \n"
                         "        magnetic materials and phenomena, including ferro, ferri and \n"
                         "        antiferromagnets, core-shell nanoparticles, ultrafast spin dynamics,\n"
                         "        magnetic recording media, heat assisted magnetic recording, exchange\n"
                         "        bias, magnetic multilayer films and complete devices." << std::endl << std::endl;
            std::cout << "        Formatted and hyperlinked documentation, as well as tutorials, \n"
                         "        references and more is available at http://vampire.york.ac.uk/ ." << std::endl << std::endl;
            std::cout << "\x1b[1m" << "HELP TOPICS" << "\x1b[0m" << std::endl;
            std::cout << "        The following list provides an overview of the help documentation.\n"
                         "        Keywords are underlined and detailed information about each one can\n"
                         "        be accessed using vampire --help<keyword>." << std::endl << std::endl;
            std::cout << "        \x1b[4moverview\x1b[0m - overview of the program features and capabilities" << std::endl;
            std::cout << "        \x1b[4mdimensions\x1b[0m - control of system dimensions" << std::endl;
            std::cout << "        \x1b[4mcreate\x1b[0m - control of structural attributes, eg shape, crystal structure" << std::endl;
            std::cout << "        \x1b[4msimulate\x1b[0m - control of simulation attributes, eg integrator, timestep" << std::endl;
            std::cout << "        \x1b[4moutput\x1b[0m - control of output data and statistics" << std::endl;
            std::cout << "        \x1b[4mmaterials\x1b[0m - definition of magnetic material properties" << std::endl;
            std::cout << std::endl;
            std::cout << "\x1b[1m" << "AUTHORS" << "\x1b[0m" << std::endl;
            std::cout << "        Vampire is principally developed by Dr. Richard Evans, however it would not\n"
                         "        have been possible without the contributions from many members of the \n"
                         "        Computational Magnetism Group at York. Special thanks go to Weijia Fan,\n"
                         "        Phanwadee Chureemart, Thomas Ostler, Joe Barker and Roy Chantrell for \n"
                         "        their support, development and beta testing of the code." << std::endl << std::endl;
            std::cout << "\x1b[1m" << "REPORTING BUGS" << "\x1b[0m" << std::endl;
            std::cout << "        Please report bugs and issues directly to the author via email to\n"
                         "        richard.evans@york.ac.uk or on the github bug tracker at\n"
                         "        https://github.com/richard-evans/vampire/issues/" << std::endl << std::endl;
            return EXIT_SUCCESS;
         }
      }
		else{
			std::cerr << "Error - unknown command line parameter \'" << sw << "\'" << std::endl;
			return EXIT_FAILURE;
		}
	}
	
	// For parallel execution intialise MPI
	#ifdef MPICF
		vmpi::initialise();
	#endif 
	
	// Initialise log file
   vout::zLogTsInit(std::string(argv[0]));
   vio::initialise_vlog_timestamp(std::string(argv[0]));

	// Output Program Header
	if(vmpi::my_rank==0){
		std::cout << "                                                _          " << std::endl;
		std::cout << "                                               (_)         " << std::endl;
		std::cout << "                    __   ____ _ _ __ ___  _ __  _ _ __ ___ " << std::endl;
		std::cout << "                    \\ \\ / / _` | '_ ` _ \\| '_ \\| | '__/ _ \\" << std::endl;
		std::cout << "                     \\ V / (_| | | | | | | |_) | | | |  __/" << std::endl;
		std::cout << "                      \\_/ \\__,_|_| |_| |_| .__/|_|_|  \\___|" << std::endl;
		std::cout << "                                         | |               " << std::endl;
		std::cout << "                                         |_|               " << std::endl;
		std::cout << std::endl;
		std::cout << "                       Version 2.0 " << __DATE__ << " " << __TIME__ << std::endl;
		std::cout << std::endl;

		std::cout << "  Licensed under the GNU Public License(v2). See licence file for details." << std::endl;
		std::cout << std::endl;
		std::cout << "  Contributors: Richard F L Evans, Weijia Fan, Joe Barker, " << std::endl;
		std::cout << "                Thomas Ostler, Phanwadee Chureemart, Roy W Chantrell" << std::endl;
		std::cout << " " << std::endl;
		#ifdef COMP	
		std::cout << "                Compiled with:  " << COMP << std::endl;
		#endif 
		std::cout << "                Compiler Flags: ";
		#ifdef CUDA
		std::cout << "CUDA ";
		#endif
		#ifdef MPICF
		std::cout << "MPI ";
		#endif
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "================================================================================" << std::endl;
		time_t rawtime = time(NULL);
		struct tm * timeinfo = localtime(&rawtime);
		std::cout<<asctime(timeinfo);
	}


	#ifdef MPICF
		vmpi::hosts();
	#endif 
  
#ifdef MPICF
  // nullify non root cout stream
  if(vmpi::my_rank!=0){
    vout::nullify(std::cout);
  }
  #endif
  // redirect std::err to file
  //std::stringstream ss;
  //ss << "vampire."<<vmpi::my_rank<<".err";
  //vout::redirect(std::cerr,ss.str());
	

	// Initialise system
	mp::initialise(infile);

	// Create system
	cs::create();

	// Simulate system
	sim::run();

	// Finalise MPI
	#ifdef MPICF
		vmpi::finalise();
		// concatenate log, sort, and appen departure message.
		if(vmpi::num_processors!=1) system("ls zlog.* | xargs cat | sort -n > zlog");
	#endif

	zlog << zTs() << "Program ended gracefully. Exiting." << std::endl;

	return EXIT_SUCCESS;
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
