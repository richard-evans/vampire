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
      	    terminaltextcolor(RED);
            std::cerr << "Error - no file specified for \'-f\' command line option" << std::endl;
            terminaltextcolor(WHITE);
            return EXIT_FAILURE;
         }
      }
      else{
         terminaltextcolor(RED);
         std::cerr << "Error - unknown command line parameter \'" << sw << "\'" << std::endl;
         terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }
   }

   // For parallel execution intialise MPI
   #ifdef MPICF
      vmpi::initialise();
   #endif

   // Initialise log file
   vout::zLogTsInit(std::string(argv[0]));
   

   // Output Program Header
   if(vmpi::my_rank==0){
      terminaltextcolor(YELLOW);
      std::cout << "                                                _          " << std::endl;
      std::cout << "                                               (_)         " << std::endl;
      std::cout << "                    __   ____ _ _ __ ___  _ __  _ _ __ ___ " << std::endl;
      std::cout << "                    \\ \\ / / _` | '_ ` _ \\| '_ \\| | '__/ _ \\" << std::endl;
      std::cout << "                     \\ V / (_| | | | | | | |_) | | | |  __/" << std::endl;
      std::cout << "                      \\_/ \\__,_|_| |_| |_| .__/|_|_|  \\___|" << std::endl;
      std::cout << "                                         | |               " << std::endl;
      std::cout << "                                         |_|               " << std::endl;
      std::cout << std::endl;
      std::cout << "                      Version 3.0.2 " << __DATE__ << " " << __TIME__ << std::endl;
      std::cout << std::endl;
	   terminaltextcolor(WHITE);
      std::cout << "  Licensed under the GNU Public License(v2). See licence file for details." << std::endl;
      std::cout << std::endl;
      std::cout << "  Lead Developer: Richard F L Evans <richard.evans@york.ac.uk>" << std::endl;
      std::cout << std::endl;
      std::cout << "  Contributors: Weijia Fan, Phanwadee Chureemart, Joe Barker, " << std::endl;
      std::cout << "                Thomas Ostler, Andreas Biternas, Roy W Chantrell" << std::endl;
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

   // Initialise system
   mp::initialise(infile);

   // Create system
   cs::create();

   // Simulate system
   sim::run();

   // Finalise MPI
   #ifdef MPICF
      vmpi::finalise();
      // concatenate log, sort, and append departure message.
      #ifdef WIN_COMPILE
         if(vmpi::num_processors!=1 && vmpi::my_rank==0) system("type log.* 2>NUL | sort > log");
      #else
         if(vmpi::num_processors!=1 && vmpi::my_rank==0) system("ls log.* | xargs cat | sort -n > log");
      #endif
   #endif

   zlog << zTs() << "Simulation ended gracefully." << std::endl;
   terminaltextcolor(GREEN);
   std::cout << "Simulation ended gracefully." << std::endl;
   terminaltextcolor(WHITE);   

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
