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
#ifndef VIO_H_
#define VIO_H_

#include <fstream>
#include <string>
#include <iostream>
#include <time.h>
#include <sys/types.h>
#ifdef WIN_COMPILE
  #include <process.h>
  #include <windows.h>
#else
  #include <unistd.h>
#endif

#include <vmpi.hpp>

#ifdef MPICF
struct null_streambuf 
: public std::streambuf 
{ 
  void overflow(char c) 
  { 
  } 
}; 
#endif

// Global Output Streams
extern std::ofstream zinfo;
extern std::ofstream zmag;
extern std::ofstream zgrain;
extern std::ofstream zlog;

enum textcolor { 
		RED=31,
		GREEN=32,
		YELLOW=33, 
		BLUE=34,
		PURPLE=35,
		WHITE=0

};
void terminaltextcolor(enum textcolor );

// global timestap function
extern std::string zTs();

namespace vin{
	extern int read(std::string const);
}

namespace vout{
	
	extern std::vector<unsigned int> file_output_list;
	extern std::vector<unsigned int> screen_output_list;
	extern std::vector<unsigned int> grain_output_list;
	
	extern int output_grain_rate;
   extern bool gnuplot_array_format;
	
	extern bool output_atoms_config;
	extern int output_atoms_config_rate;
	
	extern double atoms_output_min[3];
	extern double atoms_output_max[3];
	
	extern bool output_cells_config;
	extern int output_cells_config_rate;
	
	extern bool output_grains_config;
	extern int output_config_grain_rate;
	
	//extern bool output_povray;
	//extern int output_povray_rate;

	//extern bool output_povray_cells;
	//extern int output_povray_cells_rate;
	
	extern void data();
	extern void config();
	extern void zLogTsInit(std::string);
	
	//extern int pov_file();

	void redirect(std::ostream& strm, std::string filename);
	void nullify(std::ostream& strm);  

}
#endif /*VIO_H_*/
