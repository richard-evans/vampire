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
#include "cells.hpp"
#include "demag.hpp"
#include "errors.hpp"
#include "gpu.hpp"
#include "grains.hpp"
#include "ltmp.hpp"
#include "voronoi.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
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

// Global output filestreams
std::ofstream zinfo;
std::ofstream zlog;
std::ofstream zmag;
std::ofstream zgrain;


void terminaltextcolor(enum textcolor color){
#ifdef WIN_COMPILE
 int fincolor=15;
 if(color==RED) fincolor=12; if(color==GREEN) fincolor=10; if(color==YELLOW) fincolor=14;
 if(color==BLUE) fincolor=9; if(color==PURPLE) fincolor=13;
 SetConsoleTextAttribute(GetStdHandle( STD_OUTPUT_HANDLE ), fincolor);
#else
  std::ostringstream fincolor;
  fincolor<< color;
   std::cout << "\033["<<fincolor.str()<<"m";
   std::cerr << "\033["<<fincolor.str()<<"m";
#endif
};

namespace vout{

   std::string zLogProgramName; /// Program Name
   std::string zLogHostName; /// Host Name
   bool        zLogInitialised=false; /// Initialised flag
   #ifdef WIN_COMPILE
      int      zLogPid; /// Process ID
   #else
      pid_t    zLogPid; /// Process ID
   #endif

   void zLogTsInit(std::string tmp){

      // Get program name and process ID
      std::string tmprev;
      int linelength = tmp.length();

      // set character triggers
      const char* key="/";	/// Word identifier

      // copy characters after last /
      for(int i=linelength-1;i>=0;i--){

         char c=tmp.at(i);

         if(c != *key){
            tmprev.push_back(c);
         }
         else break;
      }

      //reverse read into program name
      linelength=tmprev.size();
      for(int i=linelength-1;i>=0;i--){
         char c=tmprev.at(i);
         zLogProgramName.push_back(c);
      }

      // Get hostname
      char loghostname [80];
      #ifdef WIN_COMPILE
         DWORD sizelhn = sizeof ( loghostname );
         int GHS=!GetComputerName(loghostname, &sizelhn); //GetComputerName returns true when retrieves hostname
      #else
         int GHS=gethostname(loghostname, 80);
      #endif
	  terminaltextcolor(YELLOW);
      if(GHS!=0) std::cerr << "Warning: Unable to retrieve hostname for zlog file." << std::endl;
	  terminaltextcolor(WHITE);
      zLogHostName = loghostname;

      // Now get process ID
      #ifdef WIN_COMPILE
         zLogPid = _getpid();
      #else
         zLogPid = getpid();
      #endif

      // Set unique filename for log if num_procs > 1
      std::stringstream logfn;
      if(vmpi::num_processors==1) logfn << "log";
      else logfn << "log."<<vmpi::my_rank;

      // Open log filename
      std::string log_file = logfn.str();
      const char* log_filec = log_file.c_str();
      zlog.open(log_filec);

      // Mark as initialised;
      zLogInitialised=true;

      zlog << zTs() << "Logfile opened" << std::endl;

      return;
   }

}

/// @brief Function to output timestamp to stream
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2012. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    19/04/2012
///
/// @return TS
///
/// @internal
///	Created:		19/04/2012
///	Revision:	  ---
///=====================================================================================
///
std::string zTs(){

  std::string NullString;
  NullString="";

	if(vout::zLogInitialised==true){
		std::ostringstream Ts;

		// varibale for time
		time_t seconds;

		// get current time
		seconds = time (NULL);
		struct tm * timeinfo;
		char logtime [80];

		timeinfo = localtime ( &seconds );
		// Format time string
		strftime (logtime,80,"%Y-%m-%d %X ",timeinfo);

		Ts << logtime << vout::zLogProgramName << " [" << vout::zLogHostName << ":" << vout::zLogPid << ":"<< vmpi::my_rank << "] ";

		return Ts.str();

	}
	else{
		terminaltextcolor(RED);
		std::cerr << "Error! - zlog not initialised, exiting" << std::endl;
		// This can be recursive - vexit calls zTs()
      //err::vexit();
      // Exit manually
      std::cerr << "Fatal error: Aborting program. See log file for details." << std::endl;
	  terminaltextcolor(WHITE);
      exit(EXIT_FAILURE);
	}

	return NullString;
}

///-------------------------------------------------------
/// Function to write header information about simulation
///-------------------------------------------------------
void write_output_file_header(std::ofstream& ofile, std::vector<unsigned int>& file_output_list){

   //------------------------------------
   // Determine current time
   //------------------------------------
   time_t seconds;

   // get time now
   seconds = time (NULL);
   struct tm * timeinfo;
   char oftime [80];

   timeinfo = localtime ( &seconds );
   // format time string
   strftime (oftime,80,"%Y-%m-%d %X ",timeinfo);

   //------------------------------------
   // Determine current directory
   //------------------------------------
   char directory [256];

   #ifdef WIN_COMPILE
      _getcwd(directory, sizeof(directory));
   #else
      getcwd(directory, sizeof(directory));
   #endif

   //------------------------------------
   // Output output file header
   //------------------------------------
   ofile << "#----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
   ofile << "# " << "Output file for vampire simulation" << std::endl;
   ofile << "# " << "  time       : " << oftime << "    process id : " << vout::zLogPid << std::endl;
   ofile << "# " << "  hostname   : " << vout::zLogHostName << std::endl;
   ofile << "# " << "  path       : " << directory << std::endl;
   ofile << "#----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
   //ofile << "# time" << "\t" << "temperature" << "\t" <<  "|m|" << "\t" << "..." << std::endl; // to be concluded...

   return;

}

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
  int read_mat_file(std::string const, int const);
int match_create(std::string const, std::string const, std::string const, int const);
int match_dimension(std::string const, std::string const, std::string const, int const);
int match_sim(std::string const, std::string const, std::string const, int const);
int match_vout_list(std::string const, std::string const, int const, std::vector<unsigned int> &);
int match_vout_grain_list(std::string const, std::string const, int const, std::vector<unsigned int> &);
int match_material(string const, string const, string const, int const, int const, int const);
int match_config(string const, string const, int const);

// Function to extract all variables from a string and return a vector
std::vector<double> DoublesFromString(std::string value){

	// array for storing variables
	std::vector<double> array(0);

	// set source for ss
	std::istringstream source(value);

	// double variable to store values
	double temp = 0.0;

	// string to store text
	std::string field;

	// loop over all comma separated values
	while(getline(source,field,',')){

		// convert string to ss
		std::stringstream fs(field);

		// read in variable
		fs >> temp;

		// push data value back to array
		array.push_back(temp);

	}

	// return values to calling function
	return array;

}

///
/// Function to check for correct unit type and valid variable range
///-----------------------------------------------------------------------
///
void check_for_valid_value(double& value, /// value of variable as in input file
									std::string word, /// input file keyword
									int line, /// input file line
									std::string prefix, /// input file prefix
									std::string unit, /// unit specified in input file
									std::string unit_type, /// expected unit type
									double range_min, /// acceptable minimum value for variable
									double range_max, /// acceptable maximum value for variable
									std::string input_file_type, ///input file name
									std::string range_text) /// customised text
{

	// Define test unit
	std::string test_unit_type=unit_type;

	// Define integer for unit conversion status
	int convert_status=0;

	// If no unit given, assume internal, otherwise convert to internal units
	if(unit.size() != 0) convert_status = units::convert(unit,value,test_unit_type);

	// Test for valid conversion
	if(convert_status==EXIT_FAILURE){
		terminaltextcolor(RED);
		std::cerr << "Error: Unit \'" << unit << "\' specified on line " << line << " of " << input_file_type << " file is not a valid unit." << std::endl;
		terminaltextcolor(WHITE);
		zlog << zTs() << "Error: Unit \'" << unit << "\' specified on line " << line << " of " << input_file_type << " file is not a valid unit." << std::endl;
		err::vexit();
	}

	// Test for change in unit type in case of wrong unit type
	if(unit_type!=test_unit_type){
		terminaltextcolor(RED);
		std::cerr << "Error: Unit \'" << unit << "\' of type \'" << test_unit_type << "\' specified on line " << line << " of " << input_file_type << " is invalid for parameter " << prefix << word << "."<< std::endl;
		terminaltextcolor(WHITE);
		zlog << zTs() << "Error: Unit \'" << unit << "\' of type \'" << test_unit_type << "\' specified on line " << line << " of " << input_file_type << " is invalid for parameter " << prefix << word << "."<< std::endl;
		err::vexit();
	}

	// Check for valid range
	if((fabs(value)<range_min) || (fabs(value)>range_max)){
		terminaltextcolor(RED);
		std::cerr << "Error: " << prefix << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
		terminaltextcolor(WHITE);
		zlog << zTs() << "Error: " << prefix << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
		err::vexit();
	}

	// Success - input is sane!
	return;

}

///
/// Function to check for valid int variable range
///-----------------------------------------------------------------------
///
void check_for_valid_int(  int& value, /// value of variable as in input file
                           std::string word, /// input file keyword
                           int line, /// input file line
                           std::string prefix, /// input file prefix
                           int range_min, /// acceptable minimum value for variable
                           int range_max, /// acceptable maximum value for variable
                           std::string input_file_type, ///input file name
                           std::string range_text) /// customised text
{

   // Check for valid range
   if((value<range_min) || (value>range_max)){
	  terminaltextcolor(RED);
      std::cerr << "Error: " << prefix << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
      terminaltextcolor(WHITE);
	  zlog << zTs() << "Error: " << prefix << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
      err::vexit();
   }

   // Success - input is sane!
   return;

}

///
/// Overloaded function to check for valid uint variable range
///-----------------------------------------------------------------------
///
void check_for_valid_int(  unsigned int& value, /// value of variable as in input file
                           std::string word, /// input file keyword
                           int line, /// input file line
                           std::string prefix, /// input file prefix
                           unsigned int range_min, /// acceptable minimum value for variable
                           unsigned int range_max, /// acceptable maximum value for variable
                           std::string input_file_type, ///input file name
                           std::string range_text) /// customised text
{

   // Check for valid range
   if((value<range_min) || (value>range_max)){
	  terminaltextcolor(RED);
      std::cerr << "Error: " << prefix << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
      terminaltextcolor(WHITE);
	  zlog << zTs() << "Error: " << prefix << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
      err::vexit();
   }

   // Success - input is sane!
   return;

}

///-----------------------------------------------------------------------
/// Function to check for valid boolean
///
/// (c) R F L Evans 2013
///
/// If input is invalid, then function will output error message and
/// program will exit from here. Otherwise returns a sanitised bool.
///
///-----------------------------------------------------------------------
bool check_for_valid_bool( std::string value, /// variable as in input file
                           std::string word, /// input file keyword
                           int line, /// input file line
                           std::string prefix, /// input file prefix
                           std::string input_file_type) ///input file name
{
   // Define string constants
   const std::string t="true";
   const std::string f="false";
   const std::string b="";

   // Check for three possible correct answers
   if(value==t) return true;
   if(value==f) return false;
   if(value==b) return true;

   // Invalid input - print error and exit
   terminaltextcolor(RED);
   std::cerr << "Error: " << prefix << word << " on line " << line << " of " << input_file_type << " file must be true or false." << std::endl;
   terminaltextcolor(WHITE);
   zlog << zTs() << "Error: " << prefix << word << " on line " << line << " of " << input_file_type << " file must be true or false." << std::endl;
   err::vexit();

   return false;
}

///
/// Function to check for correct 3-component vector and ensure length of 1
///-------------------------------------------------------------------------
///
void check_for_valid_unit_vector(std::vector<double>& u, /// unit vector
                           std::string word, /// input file keyword
                           int line, /// input file line
                           std::string prefix, /// input file prefix
                           std::string input_file_type) ///input file name
{

   // check size
   if(u.size()!=3){
	  terminaltextcolor(RED);
      std::cerr << "Error: unit-vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must have three values." << std::endl;
      terminaltextcolor(WHITE);
	  zlog << zTs() << "Error: unit-vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must have three values." << std::endl;
      err::vexit();
   }

   // Normalise
   double ULength=sqrt(u.at(0)*u.at(0)+u.at(1)*u.at(1)+u.at(2)*u.at(2));

   // Check for correct length unit vector
   if(ULength < 1.0e-9){
	  terminaltextcolor(RED);
      std::cerr << "Error: unit-vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must be normalisable (possibly all zero)." << std::endl;
      terminaltextcolor(WHITE);
	  zlog << zTs() << "Error: unit-vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must be normalisable (possibly all zero)." << std::endl;
      err::vexit();
   }
   u.at(0)/=ULength;
   u.at(1)/=ULength;
   u.at(2)/=ULength;

   // Success - input is sane!
   return;

}

///
/// Function to check for correct 3-component vector and ensure length of 1
///-------------------------------------------------------------------------
///
void check_for_valid_vector(std::vector<double>& u, /// unit vector
                           std::string word, /// input file keyword
                           int line, /// input file line
                           std::string prefix, /// input file prefix
                           std::string input_file_type) ///input file name
{

   // check size
   if(u.size()!=3){
	  terminaltextcolor(RED);
      std::cerr << "Error: vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must have three values." << std::endl;
      terminaltextcolor(WHITE);
	  zlog << zTs() << "Error: vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must have three values." << std::endl;
      err::vexit();
   }
   // Check for valid range
   if(fabs(u.at(0)) >1.e10){
	  terminaltextcolor(RED);
      std::cerr << "Error: first element of vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must be between +/- 1e10." << std::endl;
      terminaltextcolor(WHITE);
	  zlog << zTs() << "Error: first element of vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must be between +/- 1e10." << std::endl;
      err::vexit();
   }
   if(fabs(u.at(1)) >1.e10){
	  terminaltextcolor(RED);
      std::cerr << "Error: second element of vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must be between +/- 1e10." << std::endl;
      terminaltextcolor(WHITE);
	  zlog << zTs() << "Error: second element of vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must be between +/- 1e10." << std::endl;
      err::vexit();
   }
   if(fabs(u.at(2)) >1.e10){
	  terminaltextcolor(RED);
      std::cerr << "Error: third element of vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must be between +/- 1e10." << std::endl;
      terminaltextcolor(WHITE);
	  zlog << zTs() << "Error: third element of vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must be between +/- 1e10." << std::endl;
      err::vexit();
   }

   // Success - input is sane!
   return;

}
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

	// Print informative message to zlog file
	zlog << zTs() << "Opening main input file \"" << filename << "\"." << std::endl;

	// Open file read only
	inputfile.open(filename.c_str());

	// Check for opening
	if(!inputfile.is_open()){
	  terminaltextcolor(RED);
	  std::cerr << "Error opening main input file \"" << filename << "\". File does not exist!" << std::endl;
	  terminaltextcolor(WHITE);
	  zlog << zTs() << "Error: Main input file \"" << filename << "\" cannot be opened or does not exist." << std::endl;
	  zlog << zTs() << "If file exists then check file permissions to ensure it is readable by the user." << std::endl;
	  err::vexit();   // return to calling function for error checking or message
	}

        // Print informative message to zlog file
	zlog << zTs() << "Parsing system parameters from main input file." << std::endl;

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

		// clear carriage return for dos formatted files
		line.erase(remove(line.begin(), line.end(), '\r'), line.end());

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

			// if character is not ":" or "=" or "!" or "#" interpret as key
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
			err::vexit();
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

   //-------------------------------------------------------------------
   // Call module input parameters
   //-------------------------------------------------------------------
   if(ltmp::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
   else if(gpu::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;

	//===================================================================
	// Test for create variables
	//===================================================================
	test="create";
	if(key==test){
		int frs=vin::match_create(word, value, unit, line);
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
	//===================================================================
	// Test for data file output
	//===================================================================
	else
	test="output";
	if(key==test){
		int frs=vin::match_vout_list(word, value, line, vout::file_output_list);
		return frs;
	}
	//===================================================================
	// Test for screen output
	//===================================================================
	else
	test="screen";
	if(key==test){
		int frs=vin::match_vout_list(word, value, line, vout::screen_output_list);
		return frs;
	}
	//===================================================================
	// Test for grain output
	//===================================================================
	else
	test="grain";
	if(key==test){
		int frs=vin::match_vout_grain_list(word, value, line, vout::grain_output_list);
		return frs;
	}
	//===================================================================
	// Test for config output
	//===================================================================
	else
	test="config";
	if(key==test){
		int frs=vin::match_config(word, value, line);
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
			  read_mat_file(matfile,line);
				return EXIT_SUCCESS;
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error - empty filename in control statement \'material:" << word << "\' on line " << line << " of input file" << std::endl;
				terminaltextcolor(WHITE);
				return EXIT_FAILURE;
			}
		}
		//-------------------------------------------------------------------
		// Get unit cell filename
		//-------------------------------------------------------------------
		test="unit-cell-file";
		if(word==test){
			std::string matfile=value;
			// strip quotes
			matfile.erase(remove(matfile.begin(), matfile.end(), '\"'), matfile.end());
			test="";
			if(matfile!=test){
				//std::cout << matfile << std::endl;
				cs::unit_cell_file=matfile;
				return EXIT_SUCCESS;
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error - empty filename in control statement \'material:" << word << "\' on line " << line << " of input file" << std::endl;
				terminaltextcolor(WHITE);
				return EXIT_FAILURE;
			}
		}
		else{
			terminaltextcolor(RED);
			std::cerr << "Error - Unknown control statement \'material:" << word << "\' on line " << line << " of input file" << std::endl;
			terminaltextcolor(WHITE);
			return EXIT_FAILURE;
		}
	}
	else
		terminaltextcolor(RED);
		std::cerr << "Error - Unknown control statement \'" << key <<":"<< word << "\' on line " << line << " of input file" << std::endl;
		terminaltextcolor(WHITE);
		return EXIT_FAILURE;

} // end of match function

int match_create(string const word, string const value, string const unit, int const line){
   ///-------------------------------------------------------------------
   /// system_creation_flags[1] - Set system particle shape
   ///-------------------------------------------------------------------

   std::string prefix="create:";

   // cs::system_creation_flags needs refactoring for readability and bug resistance
   std::string test="full";
   if(word==test){
      cs::system_creation_flags[1]=0;
      return EXIT_SUCCESS;
   }
   else
   //-------------------------------------------------------------------
   test="cube";
   if(word==test){
      cs::system_creation_flags[1]=1;
      return EXIT_SUCCESS;
   }
   else
   //-------------------------------------------------------------------
   test="cylinder";
   if(word==test){
      cs::system_creation_flags[1]=2;
      return EXIT_SUCCESS;
   }
   else
   //-------------------------------------------------------------------
   test="ellipsoid";
   if(word==test){
      cs::system_creation_flags[1]=3;
      return EXIT_SUCCESS;
   }
   else
   //-------------------------------------------------------------------
   test="sphere";
   if(word==test){
      cs::system_creation_flags[1]=4;
      return EXIT_SUCCESS;
   }
   else
   //-------------------------------------------------------------------
   test="truncated-octahedron";
   if(word==test){
      cs::system_creation_flags[1]=5;
      return EXIT_SUCCESS;
   }
   else
   //-------------------------------------------------------------------
   test="tear-drop";
   if(word==test){
      cs::system_creation_flags[1]=6;
      return EXIT_SUCCESS;
   }
   else
   //-------------------------------------------------------------------
   // system_creation_flags[2] - Set system type
   //-------------------------------------------------------------------
   test="particle";
   if(word==test){
      cs::system_creation_flags[2]=0;
      return EXIT_SUCCESS;
   }
   else
   test="particle-array";
   if(word==test){
      cs::system_creation_flags[2]=1;
      return EXIT_SUCCESS;
   }
   else
   test="hexagonal-particle-array";
   if(word==test){
      cs::system_creation_flags[2]=2;
      return EXIT_SUCCESS;
   }
   else
   test="voronoi-film";
   if(word==test){
      cs::system_creation_flags[2]=3;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   else
   test="voronoi-size-variance";
   if(word==test){
      double vsd=atof(value.c_str());
      check_for_valid_value(vsd, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
      create_voronoi::voronoi_sd=vsd;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   else
   test="voronoi-row-offset";
   if(word==test){
      create_voronoi::parity=1;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   else
   test="voronoi-random-seed";
   if(word==test){
      int vs=atoi(value.c_str());
      check_for_valid_int(vs, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
      mtrandom::voronoi_seed=vs;
         return EXIT_SUCCESS;
   }
   else
   test="voronoi-rounded-grains";
   if(word==test){
      create_voronoi::rounded=true;
      return EXIT_SUCCESS;
   }
   else
   //-------------------------------------------------------------------
   test="voronoi-rounded-grains-area";
   if(word==test){
      double vsd=atof(value.c_str());
      check_for_valid_value(vsd, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
      create_voronoi::area_cutoff=vsd;
      return EXIT_SUCCESS;
   }
   else
   //-------------------------------------------------------------------
   test="particle-centre-offset"; //parity
   if(word==test){
      cs::particle_creation_parity=1;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   else
   test="crystal-structure";
   if(word==test){
      // Strip quotes
      std::string cs=value;
      cs.erase(remove(cs.begin(), cs.end(), '\"'), cs.end());
      cs::crystal_structure=cs;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   else
   test="single-spin";
   if(word==test){
      cs::single_spin=true;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   else
   test="periodic-boundaries-x";
   if(word==test){
      cs::pbc[0]=true;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   else
   test="periodic-boundaries-y";
   if(word==test){
      cs::pbc[1]=true;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   else
   test="periodic-boundaries-z";
   if(word==test){
      cs::pbc[2]=true;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   else
   test="select-material-by-height";
   if(word==test){
      cs::SelectMaterialByZHeight=true; // default
      // also check for value
      std::string VFalse="false";
      if(value==VFalse){
         cs::SelectMaterialByZHeight=false;
      }
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   else
   test="select-material-by-geometry";
   if(word==test){
      cs::SelectMaterialByGeometry=true; // default
      // also check for value
      std::string VFalse="false";
      if(value==VFalse){
         cs::SelectMaterialByGeometry=false;
      }
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="fill-core-shell-particles";
   if(word==test){
      cs::fill_core_shell=true; // default
      // also check for value
      std::string VFalse="false";
      if(value==VFalse){
         cs::fill_core_shell=false;
      }
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="interfacial-roughness";
   if(word==test){
      cs::interfacial_roughness=true; // default
      // also check for value
      std::string VFalse="false";
      if(value==VFalse){
         cs::interfacial_roughness=false;
      }
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="material-interfacial-roughness";
   if(word==test){
      cs::interfacial_roughness_local_height_field=true; // default
      // also check for value
      std::string VFalse="false";
      if(value==VFalse){
         cs::interfacial_roughness_local_height_field=false;
      }
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="interfacial-roughness-random-seed";
   if(word==test){
      unsigned int vs=atoi(value.c_str());
      cs::interfacial_roughness_random_seed=vs;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="interfacial-roughness-number-of-seed-points";
   if(word==test){
      int sc=atoi(value.c_str());
      check_for_valid_int(sc, word, line, prefix, 0, 100000,"input","0 - 100,000");
      cs::interfacial_roughness_seed_count=sc;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="interfacial-roughness-type";
   if(word==test){
      std::string loctest="peaks";
      if(value==loctest){
         cs::interfacial_roughness_type=1;
         return EXIT_SUCCESS;
      }
      else
      loctest="troughs";
      if(value==loctest){
         cs::interfacial_roughness_type=-1;
         return EXIT_SUCCESS;
      }
      else{
         cs::interfacial_roughness_type=0;
         return EXIT_SUCCESS;
      }
   }
   //--------------------------------------------------------------------
   test="interfacial-roughness-seed-radius";
   if(word==test){
      double irsr=atof(value.c_str());
      // Test for valid range
      check_for_valid_value(irsr, word, line, prefix, unit, "length", 0.0, 10000.0,"input","0.0 - 1 micrometre");
      cs::interfacial_roughness_mean_seed_radius=irsr;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="interfacial-roughness-seed-radius-variance";
   if(word==test){
      double irsrv=atof(value.c_str());
      // Test for valid range
      check_for_valid_value(irsrv, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
      cs::interfacial_roughness_seed_radius_variance=irsrv;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="interfacial-roughness-mean-height";
   if(word==test){
      double irmh=atof(value.c_str());
      // Test for valid range
      check_for_valid_value(irmh, word, line, prefix, unit, "length", 0.1, 100.0,"input","0.1 Angstroms - 10 nanometres");
      cs::interfacial_roughness_mean_seed_height=irmh;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="interfacial-roughness-maximum-height";
   if(word==test){
      double shm=atof(value.c_str());
      // Test for valid range
      check_for_valid_value(shm, word, line, prefix, unit, "length", 0.1, 100.0,"input","0.1 Angstroms - 10 nanometres");
      cs::interfacial_roughness_seed_height_max=shm;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="interfacial-roughness-height-field-resolution";
   if(word==test){
      double irhfr=atof(value.c_str());
      // Test for valid range
      check_for_valid_value(irhfr, word, line, prefix, unit, "length", 0.1, 100.0,"input","0.1 Angstroms - 10 nanometres");
      cs::interfacial_roughness_height_field_resolution=irhfr;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="multilayers";
   if(word==test){
      int nmul=atoi(value.c_str());
      // Test for valid range
      check_for_valid_int(nmul, word, line, prefix, 1, 100,"input","1 - 100, specifying the number of multilayers to be generated");
      cs::multilayers = true;
      cs::num_multilayers = nmul;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="height-categorization";
   if(word==test){
      // Test for different options
      test="default";
      if(value==test){
         // do nothing
         return EXIT_SUCCESS;
      }
      test="multilayers";
      if(value==test){
         cs::multilayer_height_category = true;
         return EXIT_SUCCESS;
      }
      else{
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'create:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"default\"" << std::endl;
         std::cerr << "\t\"multilayers\"" << std::endl;
         zlog << zTs() << "Error - value for \'create:" << word << "\' must be one of:" << std::endl;
         zlog << zTs() << "\t\"default\"" << std::endl;
         zlog << zTs() << "\t\"multilayers\"" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
   }
   //--------------------------------------------------------------------
   // keyword not found
   //--------------------------------------------------------------------
   else{
	  terminaltextcolor(RED);
      std::cerr << "Error - Unknown control statement \'create:" << word << "\' on line " << line << " of input file" << std::endl;
      terminaltextcolor(WHITE);
	  return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}

int match_dimension(string const word, string const value, string const unit, int const line){
   //-------------------------------------------------------------------
   // System dimension variables
   //-------------------------------------------------------------------
   std::string prefix="dimensions:";

   std::string test="unit-cell-size";
   if(word==test){
      double a=atof(value.c_str());
      check_for_valid_value(a, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
      cs::unit_cell_size[0]=a;
      cs::unit_cell_size[1]=a;
      cs::unit_cell_size[2]=a;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="unit-cell-size-x";
   if(word==test){
      double ax=atof(value.c_str());
      check_for_valid_value(ax, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
      cs::unit_cell_size[0]=ax;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="unit-cell-size-y";
   if(word==test){
      double ay=atof(value.c_str());
      check_for_valid_value(ay, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
      cs::unit_cell_size[1]=ay;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="unit-cell-size-z";
   if(word==test){
      double az=atof(value.c_str());
      check_for_valid_value(az, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
      cs::unit_cell_size[2]=az;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="system-size";
   if(word==test){
      double d=atof(value.c_str());
      check_for_valid_value(d, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
      cs::system_dimensions[0]=d;
      cs::system_dimensions[1]=d;
      cs::system_dimensions[2]=d;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="system-size-x";
   if(word==test){
      double dx=atof(value.c_str());
      check_for_valid_value(dx, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
      cs::system_dimensions[0]=dx;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="system-size-y";
   if(word==test){
      double dy=atof(value.c_str());
      check_for_valid_value(dy, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
      cs::system_dimensions[1]=dy;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="system-size-z";
   if(word==test){
      double dz=atof(value.c_str());
      check_for_valid_value(dz, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
      cs::system_dimensions[2]=dz;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="particle-size";
   if(word==test){
      double psize=atof(value.c_str());
      check_for_valid_value(psize, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
      cs::particle_scale=psize;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="particle-spacing";
   if(word==test){
      double pspacing=atof(value.c_str());
      check_for_valid_value(pspacing, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
      cs::particle_spacing=pspacing;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="particle-shape-factor-x";
   if(word==test){
      double sfx=atof(value.c_str());
      check_for_valid_value(sfx, word, line, prefix, unit, "none", 0.001, 1.0,"input","0.001 - 1.0");
      cs::particle_shape_factor_x=sfx;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="particle-shape-factor-y";
   if(word==test){
      double sfy=atof(value.c_str());
      check_for_valid_value(sfy, word, line, prefix, unit, "none", 0.001, 1.0,"input","0.001 - 1.0");
      cs::particle_shape_factor_y=sfy;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="particle-shape-factor-z";
   if(word==test){
      double sfz=atof(value.c_str());
      check_for_valid_value(sfz, word, line, prefix, unit, "none", 0.001, 1.0,"input","0.001 - 1.0");
      cs::particle_shape_factor_z=sfz;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="particle-array-offset-x";
   if(word==test){
      double paox=atof(value.c_str());
      check_for_valid_value(paox, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.0 - 1.0 millimetre");
      cs::particle_array_offset_x=paox;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="particle-array-offset-y";
   if(word==test){
      double paoy=atof(value.c_str());
      check_for_valid_value(paoy, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.0 - 1.0 millimetre");
      cs::particle_array_offset_y=paoy;
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="macro-cell-size";
   if(word==test){
      double cs=atof(value.c_str());
      check_for_valid_value(cs, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.0 - 1.0 millimetre");
      cells::size=cs;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   else{
	  terminaltextcolor(RED);
      std::cerr << "Error - Unknown control statement \'dimensions:"<< word << "\' on line " << line << " of input file" << std::endl;
      terminaltextcolor(WHITE);
	  return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}

int match_sim(string const word, string const value, string const unit, int const line){

   std::string prefix="sim:";

   //-------------------------------------------------------------------
   // System simulation variables
   //-------------------------------------------------------------------
   std::string test="integrator";
   if(word==test){
      test="llg-heun";
      if(value==test){
         sim::integrator=0;
         return EXIT_SUCCESS;
      }
      test="monte-carlo";
      if(value==test){
         sim::integrator=1;
         return EXIT_SUCCESS;
      }
      test="llg-midpoint";
      if(value==test){
         sim::integrator=2;
         return EXIT_SUCCESS;
      }
      test="constrained-monte-carlo";
      if(value==test){
         sim::integrator=3;
         return EXIT_SUCCESS;
      }
      test="hybrid-constrained-monte-carlo";
      if(value==test){
         sim::integrator=4;
         return EXIT_SUCCESS;
      }
      else{
		 terminaltextcolor(RED);
         std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"llg-heun\"" << std::endl;
         std::cerr << "\t\"llg-midpoint\"" << std::endl;
         std::cerr << "\t\"monte-carlo\"" << std::endl;
         std::cerr << "\t\"constrained-monte-carlo\"" << std::endl;
		 terminaltextcolor(WHITE);
         err::vexit();
      }
   }
   //-------------------------------------------------------------------
   test="program";
   if(word==test){
      test="benchmark";
      if(value==test){
         sim::program=0;
         return EXIT_SUCCESS;
      }
      test="time-series";
      if(value==test){
         sim::program=1;
         return EXIT_SUCCESS;
      }
      test="hysteresis-loop";
      if(value==test){
         sim::program=2;
         return EXIT_SUCCESS;
      }
      test="static-hysteresis-loop";
      if(value==test){
         sim::program=3;
         return EXIT_SUCCESS;
      }
      test="curie-temperature";
      if(value==test){
         sim::program=4;
         return EXIT_SUCCESS;
      }
      test="field-cool";
      if(value==test){
         sim::program=5;
         return EXIT_SUCCESS;
      }
      test="laser-pulse";
      if(value==test){
         sim::program=6;
         return EXIT_SUCCESS;
      }
      test="hamr-simulation";
      if(value==test){
         sim::program=7;
         return EXIT_SUCCESS;
      }
      test="cmc-anisotropy";
      if(value==test){
         sim::program=8;
         return EXIT_SUCCESS;
      }
      test="hybrid-cmc";
      if(value==test){
         sim::program=9;
         return EXIT_SUCCESS;
      }
      test="reverse-hybrid-cmc";
      if(value==test){
         sim::program=10;
         return EXIT_SUCCESS;
      }
      test="LaGrange-Multiplier";
      if(value==test){
         sim::program=11;
         return EXIT_SUCCESS;
      }
      test="partial-hysteresis-loop";
      if(value==test){
         sim::program=12;
         return EXIT_SUCCESS;
      }
      test="localised-temperature-pulse";
      if(value==test){
         sim::program=13;
         return EXIT_SUCCESS;
      }
      test="effective-damping";
      if(value==test){
         sim::program=14;
         return EXIT_SUCCESS;
      }
      test="diagnostic-boltzmann";
      if(value==test){
         sim::program=50;
         return EXIT_SUCCESS;
      }
      else{
		 terminaltextcolor(RED);
         std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"benchmark\"" << std::endl;
         std::cerr << "\t\"time-series\"" << std::endl;
         std::cerr << "\t\"hysteresis-loop\"" << std::endl;
         std::cerr << "\t\"static-hysteresis-loop\"" << std::endl;
         std::cerr << "\t\"curie-temperature\"" << std::endl;
         std::cerr << "\t\"field-cool\"" << std::endl;
         std::cerr << "\t\"laser-pulse\"" << std::endl;
         std::cerr << "\t\"cmc-anisotropy\"" << std::endl;
         std::cerr << "\t\"hybrid-cmc\"" << std::endl;
         std::cerr << "\t\"reverse-hybrid-cmc\"" << std::endl;
         std::cerr << "\t\"localised-temperature-pulse\"" << std::endl;
         terminaltextcolor(WHITE);
		 err::vexit();
      }
   }
   //-------------------------------------------------------------------
   test="enable-dipole-fields";
   if(word==test){
      sim::hamiltonian_simulation_flags[4]=1;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="enable-fmr-field";
   if(word==test){
      sim::hamiltonian_simulation_flags[5]=1;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="enable-fast-dipole-fields";
   if(word==test){
      demag::fast=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="dipole-field-update-rate";
   if(word==test){
      int dpur=atoi(value.c_str());
      check_for_valid_int(dpur, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
      demag::update_rate=dpur;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="enable-surface-anisotropy";
   if(word==test){
      sim::surface_anisotropy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="surface-anisotropy-threshold";
   if(word==test){
      // test for native keyword
      test="native";
      if(value==test){
         sim::NativeSurfaceAnisotropyThreshold=true;
         return EXIT_SUCCESS;
      }
      int sat=atoi(value.c_str());
      // Test for valid range
      check_for_valid_int(sat, word, line, prefix, 0, 1000000000,"input","0 - 1,000,000,000");
      sim::surface_anisotropy_threshold=sat;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="surface-anisotropy-nearest-neighbour-range";
   if(word==test){
      // Test for valid range
      double r=atof(value.c_str());
      check_for_valid_value(r, word, line, prefix, unit, "length", 0.0, 1.0e9,"input","0.0 - 1,000,000,000");
      sim::nearest_neighbour_distance=r;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="time-step";
   if(word==test){
      double dt=atof(value.c_str());
      check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
      mp::dt_SI=dt;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="total-time-steps";
   if(word==test){
      int tt=atoi(value.c_str());
      check_for_valid_int(tt, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
      sim::total_time=tt;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="loop-time-steps";
   if(word==test){
      int tt=atoi(value.c_str());
      check_for_valid_int(tt, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
      sim::loop_time=tt;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="partial-time-steps";
   if(word==test){
      int tt=atoi(value.c_str());
      check_for_valid_int(tt, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
      terminaltextcolor(YELLOW);
      std::cout << "Warning: Keyword \'partial-time-steps\' is deprecated and may be removed in a future release. Please use \'time-steps-increment\' instead." << std::endl;
      terminaltextcolor(WHITE);
	  sim::partial_time=tt;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="time-steps-increment";
   if(word==test){
      int tt=atoi(value.c_str());
      check_for_valid_int(tt, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
      sim::partial_time=tt;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="equilibration-time-steps";
   if(word==test){
      int tt=atoi(value.c_str());
      check_for_valid_int(tt, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
      sim::equilibration_time=tt;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="simulation-cycles";
   if(word==test){
      int r=atoi(value.c_str());
      check_for_valid_int(r, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
      sim::runs=r;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="maximum-temperature";
   if(word==test){
      double T=atof(value.c_str());
      check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
      sim::Tmax=T;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="minimum-temperature";
   if(word==test){
      double T=atof(value.c_str());
      check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
      sim::Tmin=T;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="equilibration-temperature";
   if(word==test){
      double T=atof(value.c_str());
      check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
      sim::Teq=T;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="temperature";
   if(word==test){
      double T=atof(value.c_str());
      check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
      sim::temperature=T;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="temperature-increment";
   if(word==test){
      double T=atof(value.c_str());
      check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
      sim::delta_temperature=T;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="cooling-time";
   if(word==test){
      double T=atof(value.c_str());
      check_for_valid_value(T, word, line, prefix, unit, "time", 1.0e-18, 1.0,"input","1 attosecond - 1 s");
      sim::cooling_time=T;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="laser-pulse-temporal-profile";
   if(word==test){
      test="square";
      if(value==test){
         sim::pump_function=square;
         return EXIT_SUCCESS;
      }
      test="two-temperature";
      if(value==test){
         sim::pump_function=two_temperature;
         return EXIT_SUCCESS;
      }
      test="double-pulse-two-temperature";
      if(value==test){
         sim::pump_function=double_pump_two_temperature;
         return EXIT_SUCCESS;
      }
      test="double-pulse-square";
      if(value==test){
         sim::pump_function=double_pump_square;
         return EXIT_SUCCESS;
      }
      else{
		 terminaltextcolor(RED);
         std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"square\"" << std::endl;
         std::cerr << "\t\"double-pulse-square\"" << std::endl;
         std::cerr << "\t\"two-temperature\"" << std::endl;
         std::cerr << "\t\"double-pulse-two-temperature\"" << std::endl;
		 terminaltextcolor(WHITE);
         err::vexit();
      }
   }
   //--------------------------------------------------------------------
   test="laser-pulse-time";
   if(word==test){
      double pt=atof(value.c_str());
      check_for_valid_value(pt, word, line, prefix, unit, "time", 1.0e-18, 1.0,"input","1 attosecond - 1 s");
      sim::pump_time=pt;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="laser-pulse-power";
   if(word==test){
      double pp=atof(value.c_str());
      check_for_valid_value(pp, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
      sim::pump_power=pp;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="second-laser-pulse-time";
   if(word==test){
      double pt=atof(value.c_str());
      check_for_valid_value(pt, word, line, prefix, unit, "time", 1.0e-18, 1.0,"input","1 attosecond - 1 s");
      sim::double_pump_time=pt;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="second-laser-pulse-power";
   if(word==test){
      double pp=atof(value.c_str());
      check_for_valid_value(pp, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
      sim::double_pump_power=pp;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="second-laser-pulse-maximum-temperature";
   if(word==test){
      double T=atof(value.c_str());
      check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
      sim::double_pump_Tmax=T;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="second-laser-pulse-delay-time";
   if(word==test){
      double pd=atof(value.c_str());
      check_for_valid_value(pd, word, line, prefix, unit, "time", 0.0, 1.0,"input","0 - 1 s");
      sim::double_pump_delay=pd;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="two-temperature-heat-sink-coupling";
   if(word==test){
      double hscc=atof(value.c_str());
      check_for_valid_value(hscc, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
      sim::HeatSinkCouplingConstant=hscc;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="two-temperature-electron-heat-capacity";
   if(word==test){
      double hscc=atof(value.c_str());
      check_for_valid_value(hscc, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
      sim::TTCe=hscc;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="two-temperature-phonon-heat-capacity";
   if(word==test){
      double hscc=atof(value.c_str());
      check_for_valid_value(hscc, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
      sim::TTCl=hscc;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="two-temperature-electron-phonon-coupling";
   if(word==test){
      double hscc=atof(value.c_str());
      check_for_valid_value(hscc, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
      sim::TTG=hscc;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="cooling-function";
   if(word==test){
      test="exponential";
      if(value==test){
         sim::cooling_function_flag=0;
         return EXIT_SUCCESS;
      }
      test="gaussian";
      if(value==test){
         sim::cooling_function_flag=1;
         return EXIT_SUCCESS;
      }
      test="double-gaussian";
      if(value==test){
         sim::cooling_function_flag=2;
         return EXIT_SUCCESS;
      }
      test="linear";
      if(value==test){
         sim::cooling_function_flag=3;
         return EXIT_SUCCESS;
      }
      else{
		 terminaltextcolor(RED);
         std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"exponential\"" << std::endl;
         std::cerr << "\t\"gaussian\"" << std::endl;
         std::cerr << "\t\"double-gaussian\"" << std::endl;
         std::cerr << "\t\"linear\"" << std::endl;
		 terminaltextcolor(WHITE);
         err::vexit();
      }
   }
   //--------------------------------------------------------------------
   test="applied-field-strength";
   if(word==test){
      double H=atof(value.c_str());
      check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
      sim::H_applied=H;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="minimum-applied-field-strength";
   if(word==test){
      double H=atof(value.c_str());
      check_for_valid_value(H, word, line, prefix, unit, "field", 0.0, 1.0e3,"input","0 - 1,000 T");
      sim::Hmin=H;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="maximum-applied-field-strength";
   if(word==test){
      double H=atof(value.c_str());
      check_for_valid_value(H, word, line, prefix, unit, "field", 0.0, 1.0e3,"input","0 - 1,000 T");
      sim::Hmax=H;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="equilibration-applied-field-strength";
   if(word==test){
      double H=atof(value.c_str());
      check_for_valid_value(H, word, line, prefix, unit, "field", 0.0, 1.0e3,"input","0 - 1,000 T");
      sim::Heq=H;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="applied-field-strength-increment";
   if(word==test){
      double H=atof(value.c_str());
      check_for_valid_value(H, word, line, prefix, unit, "field", 1.0e-6, 1.0e3,"input","1 uT - 1,000 T");
      sim::Hinc=H;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="applied-field-angle-theta";
   if(word==test){
      double angle=atof(value.c_str());
      check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
      sim::applied_field_angle_theta=angle;
      sim::applied_field_set_by_angle=true;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="applied-field-angle-phi";
   if(word==test){
      double angle=atof(value.c_str());
      check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
      sim::applied_field_angle_phi=angle;
      sim::applied_field_set_by_angle=true;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="applied-field-unit-vector";
   if(word==test){
      std::vector<double> u(3);
      u=DoublesFromString(value);
      check_for_valid_unit_vector(u, word, line, prefix, "input");
      sim::H_vec[0]=u.at(0);
      sim::H_vec[1]=u.at(1);
      sim::H_vec[2]=u.at(2);
      sim::applied_field_set_by_angle=false;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="demagnetisation-factor";
   if(word==test){
      std::vector<double> u(3);
      u=DoublesFromString(value);
      check_for_valid_vector(u, word, line, prefix, "input");
      // Extra check for demagnetisation-factor Nx+Ny+Nz=1
      double sum=u.at(0)+u.at(1)+u.at(2);
      if(fabs(1.0-sum)>1.e-4){
		 terminaltextcolor(RED);
         std::cerr << "Error: sum of all elements of variable " << prefix << word << " on line " << line << " of input file must equal 1." << std::endl;
         terminaltextcolor(WHITE);
		 zlog << zTs() << "Error: sum of all elements of variable " << prefix << word << " on line " << line << " of input file must equal 1." << std::endl;
         err::vexit();
      }
      sim::demag_factor[0]=u.at(0);
      sim::demag_factor[1]=u.at(1);
      sim::demag_factor[2]=u.at(2);
      sim::ext_demag=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="mpi-mode";
   if(word==test){
      test="geometric-decomposition";
      if(value==test){
         vmpi::mpi_mode=0;
         return EXIT_SUCCESS;
      }
      test="replicated-data";
      if(value==test){
         vmpi::mpi_mode=1;
         vmpi::replicated_data_staged=false;
         return EXIT_SUCCESS;
      }
      test="replicated-data-staged";
      if(value==test){
         vmpi::mpi_mode=1;
         vmpi::replicated_data_staged=true;
         return EXIT_SUCCESS;
      }
      else{
		 terminaltextcolor(RED);
         std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"geometric-decomposition\"" << std::endl;
         std::cerr << "\t\"replicated-data\"" << std::endl;
         std::cerr << "\t\"replicated-data-staged\"" << std::endl;
		 terminaltextcolor(WHITE);
         err::vexit();
      }
   }
   //--------------------------------------------------------------------
   test="mpi-ppn";
   if(word==test){
      int ppn=atoi(value.c_str());
      check_for_valid_int(ppn, word, line, prefix, 1, 1024,"input","1 - 1024");
      vmpi::ppn=ppn;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="integrator-random-seed";
   if(word==test){
      int is=atoi(value.c_str());
      check_for_valid_int(is, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
      mtrandom::integration_seed=is;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="constraint-rotation-update";
   if(word==test){
      sim::constraint_rotation=true;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="constraint-angle-theta";
   if(word==test){
      double angle=atof(value.c_str());
      check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
      sim::constraint_theta=angle;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="constraint-angle-theta-minimum";
   if(word==test){
      double angle=atof(value.c_str());
      check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
      sim::constraint_theta_min=angle;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="constraint-angle-theta-maximum";
   if(word==test){
      double angle=atof(value.c_str());
      check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
      sim::constraint_theta_max=angle;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="constraint-angle-theta-increment";
   if(word==test){
      double angle=atof(value.c_str());
      check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
      sim::constraint_theta_delta=angle;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="constraint-angle-phi";
   if(word==test){
      double angle=atof(value.c_str());
      check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
      sim::constraint_phi=angle;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="constraint-angle-phi-minimum";
   if(word==test){
      double angle=atof(value.c_str());
      check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
      sim::constraint_phi_min=angle;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="constraint-angle-phi-maximum";
   if(word==test){
      double angle=atof(value.c_str());
      check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
      sim::constraint_phi_max=angle;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="constraint-angle-phi-increment";
   if(word==test){
      double angle=atof(value.c_str());
      check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
      sim::constraint_phi_delta=angle;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="monte-carlo-algorithm";
   if(word==test){
      // include namesapce here to access enum values
      using namespace sim;
      test="spin-flip";
      if(value==test){
         sim::mc_algorithm=spin_flip;
         return EXIT_SUCCESS;
      }
      test="uniform";
      if(value==test){
         sim::mc_algorithm=uniform;
         return EXIT_SUCCESS;
      }
      test="angle";
      if(value==test){
         sim::mc_algorithm=angle;
         return EXIT_SUCCESS;
      }
      test="hinzke-nowak";
      if(value==test){
         sim::mc_algorithm=hinzke_nowak;
         return EXIT_SUCCESS;
      }
      else{
		 terminaltextcolor(RED);
         std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"spin-flip\"" << std::endl;
         std::cerr << "\t\"uniform\"" << std::endl;
         std::cerr << "\t\"angle\"" << std::endl;
         std::cerr << "\t\"hinzke-nowak\"" << std::endl;
		 terminaltextcolor(WHITE);
         err::vexit();
      }
   }
   //-------------------------------------------------------------------
   test="save-checkpoint";
   if(word==test){
      test="end";
      if(value==test){
         sim::save_checkpoint_flag=true; // Save checkpoint
         sim::save_checkpoint_continuous_flag=false; // do not save checkpoints during simulation
         return EXIT_SUCCESS;
      }
      test="continuous";
      if(value==test){
         sim::save_checkpoint_flag=true; // Save checkpoint
         sim::save_checkpoint_continuous_flag=true; // save checkpoints during simulation
         return EXIT_SUCCESS;
      }
      else{
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"end\"" << std::endl;
         std::cerr << "\t\"continuous\"" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
   }
   //--------------------------------------------------------------------
   test="save-checkpoint-rate";
   if(word==test){
      int scr=atoi(value.c_str());
      check_for_valid_int(scr, word, line, prefix, 1, 2000000000,"input","1 - 2,000,000,000");
      sim::save_checkpoint_rate=scr;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="load-checkpoint";
   if(word==test){
      test="restart";
      if(value==test){
         sim::load_checkpoint_flag=true; // Load spin configurations
         sim::load_checkpoint_continue_flag=false; // Restart simulation with checkpoint configuration
         return EXIT_SUCCESS;
      }
      test="continue";
      if(value==test){
         sim::load_checkpoint_flag=true; // Load spin configurations
         sim::load_checkpoint_continue_flag=true; // Continue simulation from saved time with checkpoint configuration
         return EXIT_SUCCESS;
      }
      else{
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"restart\"" << std::endl;
         std::cerr << "\t\"continue\"" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
   }
   //--------------------------------------------------------------------
   else{
	  terminaltextcolor(RED);
      std::cerr << "Error - Unknown control statement \'sim:"<< word << "\' on line " << line << " of input file" << std::endl;
      terminaltextcolor(WHITE);
	  return EXIT_FAILURE;
   }


   return EXIT_SUCCESS;
}

int match_config(string const word, string const value, int const line){

   std::string prefix="config:";

   // System output config variables
   std::string test="atoms";
   if(word==test){
      vout::output_atoms_config=true;
      return EXIT_SUCCESS;
   }
   //-----------------------------------------
   test="atoms-output-rate";
   if(word==test){
      int i=atoi(value.c_str());
      check_for_valid_int(i, word, line, prefix, 1, 1000000,"input","1 - 1,000,000");
      vout::output_atoms_config_rate=i;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="atoms-minimum-x";
   if(word==test){
      double x=atof(value.c_str());
      check_for_valid_value(x, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
      vout::atoms_output_min[0]=x;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="atoms-minimum-y";
   if(word==test){
      double y=atof(value.c_str());
      check_for_valid_value(y, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
      vout::atoms_output_min[1]=y;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="atoms-minimum-z";
   if(word==test){
      double z=atof(value.c_str());
      check_for_valid_value(z, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
      vout::atoms_output_min[2]=z;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="atoms-maximum-x";
   if(word==test){
      double x=atof(value.c_str());
      check_for_valid_value(x, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
      vout::atoms_output_max[0]=x;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="atoms-maximum-y";
   if(word==test){
      double y=atof(value.c_str());
      check_for_valid_value(y, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
      vout::atoms_output_max[1]=y;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="atoms-maximum-z";
   if(word==test){
      double z=atof(value.c_str());
      check_for_valid_value(z, word, line, prefix, "", "none", 0.0, 1.0,"input","0.0 - 1.0");
      vout::atoms_output_max[2]=z;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="macro-cells";
   if(word==test){
      vout::output_cells_config=true;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="macro-cells-output-rate";
   if(word==test){
      int i=atoi(value.c_str());
      check_for_valid_int(i, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
      vout::output_cells_config_rate=i;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="identify-surface-atoms";
   if(word==test){
      sim::identify_surface_atoms=true;
      return EXIT_SUCCESS;
   }
   //-----------------------------------------
   else{
	  terminaltextcolor(RED);
      std::cerr << "Error - Unknown control statement \'config:"<< word << "\' on line " << line << " of input file" << std::endl;
	  terminaltextcolor(WHITE);
      return EXIT_FAILURE;
   }
}

int match_vout_list(string const word, string const value, int const line, std::vector<unsigned int> & output_list){

   std::string prefix="output:";

   std::string test="time-steps";
   if(word==test){
      output_list.push_back(0);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="real-time";
   if(word==test){
      output_list.push_back(1);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="temperature";
   if(word==test){
      output_list.push_back(2);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="applied-field-strength";
   if(word==test){
      output_list.push_back(3);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="applied-field-unit-vector";
   if(word==test){
      output_list.push_back(4);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="applied-field-alignment";
   if(word==test){
      stats::calculate_system_magnetization=true;
      output_list.push_back(12);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="magnetisation";
   if(word==test){
      stats::calculate_system_magnetization=true;
      output_list.push_back(5);
      return EXIT_SUCCESS;
   }
   else
   //-------------------------------------------------------------------
   test="magnetisation-length";
   if(word==test){
      stats::calculate_system_magnetization=true;
      output_list.push_back(6);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="mean-magnetisation-length";
   if(word==test){
      stats::calculate_system_magnetization=true;
      output_list.push_back(7);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="material-magnetisation";
   if(word==test){
      stats::calculate_material_magnetization=true;
      output_list.push_back(8);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="material-mean-magnetisation-length";
   if(word==test){
      stats::calculate_material_magnetization=true;
      output_list.push_back(9);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="total-torque";
   if(word==test){
      stats::calculate_torque=true;
      output_list.push_back(14);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="mean-total-torque";
   if(word==test){
      stats::calculate_torque=true;
      output_list.push_back(15);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="constraint-phi";
   if(word==test){
      stats::calculate_torque=true;
      output_list.push_back(16);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="constraint-theta";
   if(word==test){
      stats::calculate_torque=true;
      output_list.push_back(17);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="material-constraint-phi";
   if(word==test){
      stats::calculate_torque=true;
      output_list.push_back(18);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="material-constraint-theta";
   if(word==test){
      stats::calculate_torque=true;
      output_list.push_back(19);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="material-mean-torque";
   if(word==test){
      stats::calculate_torque=true;
      output_list.push_back(20);
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="mean-susceptibility";
   if(word==test){
      // Set flags for calculations of susceptibility and magnetization
      stats::calculate_system_susceptibility=true;
      stats::calculate_system_magnetization=true;
      output_list.push_back(21);
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="electron-temperature"; // identical to temperature
   if(word==test){
      output_list.push_back(2);
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="phonon-temperature";
   if(word==test){
      output_list.push_back(22);
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="material-temperature";
   if(word==test){
      output_list.push_back(23);
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="material-applied-field-strength";
   if(word==test){
      output_list.push_back(24);
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="material-fmr-field-strength";
   if(word==test){
      output_list.push_back(25);
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="material-applied-field-alignment";
   if(word==test){
      stats::calculate_material_magnetization=true;
      output_list.push_back(26);
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="total-energy";
   if(word==test){
      output_list.push_back(27);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="mean-total-energy";
   if(word==test){
      output_list.push_back(28);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="anisotropy-energy";
   if(word==test){
      output_list.push_back(29);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="mean-anisotropy-energy";
   if(word==test){
      output_list.push_back(30);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="cubic-anisotropy-energy";
   if(word==test){
      output_list.push_back(31);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="mean-cubic-anisotropy-energy";
   if(word==test){
      output_list.push_back(32);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="surface-anisotropy-energy";
   if(word==test){
      output_list.push_back(33);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="mean-surface-anisotropy-energy";
   if(word==test){
      output_list.push_back(34);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="exchange-energy";
   if(word==test){
      output_list.push_back(35);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="mean-exchange-energy";
   if(word==test){
      output_list.push_back(36);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="applied-field-energy";
   if(word==test){
      output_list.push_back(37);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="mean-applied-field-energy";
   if(word==test){
      output_list.push_back(38);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="magnetostatic-energy";
   if(word==test){
      output_list.push_back(39);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="mean-magnetostatic-energy";
   if(word==test){
      output_list.push_back(40);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="second-order-uniaxial-anisotropy-energy";
   if(word==test){
      output_list.push_back(41);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="mean-second-order-uniaxial-anisotropy-energy";
   if(word==test){
      output_list.push_back(42);
      stats::calculate_energy=true;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="height-magnetisation-normalised";
   if(word==test){
      stats::calculate_height_magnetization=true;
      output_list.push_back(43);
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="material-height-magnetisation-normalised";
   if(word==test){
      stats::calculate_material_height_magnetization=true;
      output_list.push_back(44);
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="height-magnetisation";
   if(word==test){
      stats::calculate_height_magnetization=true;
      output_list.push_back(45);
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="material-height-magnetisation";
   if(word==test){
      stats::calculate_material_height_magnetization=true;
      output_list.push_back(46);
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="mpi-timings";
   if(word==test){
      vmpi::DetailedMPITiming=true;
      output_list.push_back(60);
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="gnuplot-array-format";
   if(word==test){
      vout::gnuplot_array_format=true;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   test="output-rate";
   if(word==test){
      int r=atoi(value.c_str());
      check_for_valid_int(r, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
      vout::output_rate=r;
      return EXIT_SUCCESS;
   }

   //--------------------------------------------------------------------
   // keyword not found
   //--------------------------------------------------------------------
   else{
	   terminaltextcolor(RED);
      std::cerr << "Error - Unknown control statement "<< prefix << word << "\' on line " << line << " of input file" << std::endl;
      terminaltextcolor(WHITE);
	  return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
}
int match_vout_grain_list(string const word, string const value, int const line, std::vector<unsigned int> & output_list){

   std::string prefix="grain:";

   std::string test="time-steps";
   if(word==test){
      output_list.push_back(0);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="real-time";
   if(word==test){
      output_list.push_back(1);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="temperature";
   if(word==test){
      output_list.push_back(2);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="applied-field-strength";
   if(word==test){
      output_list.push_back(3);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="applied-field-unit-vector";
   if(word==test){
      output_list.push_back(4);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="magnetisation";
   if(word==test){
      output_list.push_back(10);
      return EXIT_SUCCESS;
   }
   else
   //-------------------------------------------------------------------
   test="mag-m";
   if(word==test){
      output_list.push_back(11);
      return EXIT_SUCCESS;
   }
   else
   //--------------------------------------------------------------------
   test="material-magnetisation";
   if(word==test){
      output_list.push_back(13);
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="electron-temperature"; // identical to temperature
   if(word==test){
      output_list.push_back(2);
      return EXIT_SUCCESS;
   }
   //-------------------------------------------------------------------
   test="phonon-temperature";
   if(word==test){
      output_list.push_back(22);
      return EXIT_SUCCESS;
   }
   else
   //-------------------------------------------------------------------
   test="output-rate";
   if(word==test){
      int r=atoi(value.c_str());
      check_for_valid_int(r, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
      vout::output_grain_rate=r;
      return EXIT_SUCCESS;
   }
   //--------------------------------------------------------------------
   // keyword not found
   //--------------------------------------------------------------------
   else{
	  terminaltextcolor(RED);
      std::cerr << "Error - Unknown control statement \'grain:" << word << "\' on line " << line << " of input file" << std::endl;
      terminaltextcolor(WHITE);
	  return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
}
// temporary array of materials for reading in material data
//std::cout << "here" << std::endl;
  std::vector<mp::materials_t> read_material(0);

int read_mat_file(std::string const matfile, int const LineNumber){

	// Declare input stream
	std::ifstream inputfile;

	// resize temporary materials array for storage of variables
	read_material.resize(mp::max_materials);
	cmc::cmc_mat.resize(mp::max_materials);

        // Print informative message to zlog file
	zlog << zTs() << "Opening material file \"" << matfile << "\"." << std::endl;

	// Open file read only
	inputfile.open(matfile.c_str());

	// Check for opening
	if(!inputfile.is_open()){
		terminaltextcolor(RED);
		std::cerr << "Error opening material file " << matfile << ". File does not exist!" << std::endl;
		terminaltextcolor(WHITE);
		zlog << zTs() << "Error: Material file \"" << matfile << "\" on line number " << LineNumber << " of input file cannot be opened or does not exist." << std::endl;
		zlog << zTs() << "If file exists then check file permissions to ensure it is readable by the user." << std::endl;
		err::vexit();   // return to calling function for error checking or message
	}
	//-------------------------------------------------------
	// Material 0
	//-------------------------------------------------------


        // Print informative message to zlog file
	zlog << zTs() << "Parsing material file for parameters." << std::endl;

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

		// remove carriage returns for dos formatted files
                line.erase(remove(line.begin(), line.end(), '\r'), line.end());

		// strip key,word,unit,value
		std::string key="";
		std::string word="";
		std::string value="";
		std::string unit="";
		std::string index="";
		int super_index=1; // Inital values *as would be read from input file*
		int sub_index=1;

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
				if((super_index>=1) && (super_index<mp::max_materials+1)){
					break;
				}
				else{
					std::cerr << "Invalid index number " << index << " on line " << line_counter << " in material input file" << std::endl;
					std::cerr << "Causes could be invalid character or outside of range, ie less than 1 or greater than max_materials=" << mp::max_materials << ", exiting" << std::endl;
					err::vexit();
				}

			}
			// For anything else
			else break;
		}
		const int end_key=last;

		//
		//err::vexit();
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
							if((sub_index<1) || (sub_index>=mp::max_materials)){
								std::cerr << "Invalid sub-index number " << index << " on line " << line_counter << " in material input file" << std::endl;
								std::cerr << "Causes could be invalid character or outside of range, ie less than 1 or greater than max_materials=" << mp::max_materials << ", exiting" << std::endl;
								err::vexit();
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
			int matchcheck = vin::match_material(word, value, unit, line_counter, super_index-1, sub_index-1);
			if(matchcheck==EXIT_FAILURE){
				err::vexit();
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

///-------------------------------------------------------------------
/// Function to match material key words
///-------------------------------------------------------------------
int match_material(string const word, string const value, string const unit, int const line, int const super_index, int const sub_index){
      std::string prefix="material:";
      //------------------------------------------------------------
      std::string test="num-materials";
      if(word==test){
         unsigned int nm = atoi(value.c_str());
         check_for_valid_int(nm, word, line, prefix, 1, 100,"material","1 - 100");
         mp::num_materials=nm;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="material-name";
      if(word==test){
         read_material[super_index].name=value;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="damping-constant";
      if(word==test){
         double damping=atof(value.c_str());
         check_for_valid_value(damping, word, line, prefix, unit, "none", 0.0, 10.0,"material","0.0 - 10.0");
         read_material[super_index].alpha=damping;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="exchange-matrix";
      if(word==test){
         double Jij=atof(value.c_str());
         check_for_valid_value(Jij, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e18");
         read_material[super_index].Jij_matrix_SI[sub_index]=-Jij; // Import exchange as field, *-1
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="atomic-spin-moment";
      if(word==test){
         double mu_s=atof(value.c_str());
         check_for_valid_value(mu_s, word, line, prefix, unit, "moment", 0.1*9.24e-24, 1e8*9.24e-24,"material","0.1 - 1e8 mu_B");
         read_material[super_index].moment_flag=true;
         read_material[super_index].mu_s_SI=mu_s;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="uniaxial-anisotropy-constant";
      if(word==test){
         double K=atof(value.c_str());
         check_for_valid_value(K, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e-18 J/atom");
         read_material[super_index].Ku1_SI=-K; // Import anisotropy as field, *-1
         // enable global anisotropy flag
         sim::UniaxialScalarAnisotropy=true;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="second-uniaxial-anisotropy-constant";
      if(word==test){
         double K=atof(value.c_str());
         check_for_valid_value(K, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e-18 J/atom");
         read_material[super_index].Ku2_SI=-K; // Import anisotropy as field, *-1
         sim::second_order_uniaxial_anisotropy=true;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="third-uniaxial-anisotropy-constant";
      if(word==test){
         double K=atof(value.c_str());
         check_for_valid_value(K, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e-18 J/atom");
         read_material[super_index].Ku3_SI=-K; // Import anisotropy as field, *-1
         sim::sixth_order_uniaxial_anisotropy=true;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="second-order-harmonic-anisotropy-constant";
      if(word==test){
         double K=atof(value.c_str());
         check_for_valid_value(K, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e-18 J/atom");
         read_material[super_index].sh2=K;
         sim::spherical_harmonics=true;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="fourth-order-harmonic-anisotropy-constant";
      if(word==test){
         double K=atof(value.c_str());
         check_for_valid_value(K, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e-18 J/atom");
         read_material[super_index].sh4=K;
         sim::spherical_harmonics=true;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="sixth-order-harmonic-anisotropy-constant";
      if(word==test){
         double K=atof(value.c_str());
         check_for_valid_value(K, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e-18 J/atom");
         read_material[super_index].sh6=K;
         sim::spherical_harmonics=true;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="lattice-anisotropy-constant";
      if(word==test){
         double Klatt=atof(value.c_str());
         // Test for valid range
         check_for_valid_value(Klatt, word, line, prefix, unit, "energy", -1.0e-18, 1.0e18,"material","-1e18 - 1e18");
         read_material[super_index].Klatt_SI=-Klatt; // Import anisotropy as field, *-1
         sim::lattice_anisotropy_flag=true;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="cubic-anisotropy-constant";
      if(word==test){
         double K=atof(value.c_str());
         // Test for valid range
         check_for_valid_value(K, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e-18 J/atom");
         read_material[super_index].Kc1_SI=-K; // Import anisotropy as field, *-1
         sim::CubicScalarAnisotropy=true;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="uniaxial-anisotropy-direction";
      if(word==test){
         // temporary storage container
         std::vector<double> u(3);

         // read values from string
         u=DoublesFromString(value);

         // check for sane input and normalise if necessary
         check_for_valid_unit_vector(u, word, line, prefix, "material");

         // Copy sanitised unit vector to material
         read_material[super_index].UniaxialAnisotropyUnitVector=u;

         // Enable global tensor anisotropy flag
         sim::TensorAnisotropy=true;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="uniaxial-anisotropy-tensor";
      if(word==test){
         std::vector<double> K;
         // read values from string
         K=DoublesFromString(value);
         // check size
         if(K.size()!=9){
			terminaltextcolor(RED);
            std::cerr << "Error in input file - material[" << super_index << "]:uniaxial-anisotropy-tensor must have nine values." << std::endl;
            terminaltextcolor(WHITE);
			zlog << zTs() << "Error in input file - material[" << super_index << "]:uniaxial-anisotropy-tensor must have nine values." << std::endl;
            return EXIT_FAILURE;
         }

         string unit_type="energy";
         // if no unit given, assume internal
         if(unit.size() != 0){
            units::convert(unit,K,unit_type);
            //read_material[super_index].anis_flag=false;
            //std::cout << "setting flag to false" << std::endl;
         }
         string str="energy";
         if(unit_type==str){
            // Copy anisotropy vector to material
            read_material[super_index].KuVec_SI=K;
            sim::TensorAnisotropy=true;
            return EXIT_SUCCESS;
         }
         else{
			 terminaltextcolor(RED);
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimensions:" << word << "\'"<< std::endl;
            terminaltextcolor(WHITE);
			err::vexit();
         }
      }
      //------------------------------------------------------------
      else
      test="surface-anisotropy-constant";
      if(word==test){
         double K=atof(value.c_str());
         // Test for valid range
         check_for_valid_value(K, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e-18 J/atom");
         read_material[super_index].Ks_SI=-K;// Import anisotropy as field, *-1
         sim::surface_anisotropy=true; // enable surface anisotropy
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="relative-gamma";
      if(word==test){
         double gr = atof(value.c_str());
         // Test for valid range
         check_for_valid_value(gr, word, line, prefix, unit, "none", 0.01, 100.0,"material"," 0.01 - 100.0");
         read_material[super_index].gamma_rel=gr;
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="initial-spin-direction";
      if(word==test){
         // first test for random spins
         test="random";
         if(value==test){
            read_material[super_index].random_spins=true;
         }
         else{
            // temporary storage container
            std::vector<double> u(3);

            // read values from string
            u=DoublesFromString(value);

            // check for sane input and normalise if necessary
            check_for_valid_unit_vector(u, word, line, prefix, "material");

            // Copy sanitised unit vector to material
            read_material[super_index].initial_spin[0]=u.at(0);
            read_material[super_index].initial_spin[1]=u.at(1);
            read_material[super_index].initial_spin[2]=u.at(2);

            // ensure random spins is unset
            read_material[super_index].random_spins=false;
         }
         // return
         return EXIT_SUCCESS;
      }
      //------------------------------------------------------------
      else
      test="material-element";
      if(word==test){
         // Test for 3 characters
         if(value.length()>3){
			terminaltextcolor(RED);
            std::cerr << "Error - element identifier on line "<<  line << " of material file must be a maximum of three characters long" << std::endl;
            terminaltextcolor(WHITE);
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
      test="geometry-file";
      if(word==test){
         // Open geometry file
         std::ifstream gfile(value.c_str());
         if(!gfile.is_open()){
			terminaltextcolor(RED);
            std::cerr << "Error - geometry file " << value.c_str() << " not found, exiting!" << std::endl;
            terminaltextcolor(WHITE);
			return EXIT_FAILURE;
         }
         gfile >> read_material[super_index].geometry;
         if((read_material[super_index].geometry<3) || (read_material[super_index].geometry>100)){
            terminaltextcolor(RED);
			std::cerr << "Error in geometry input file " << value.c_str() << " - first number must be non zero integer in the range 3-100"<< std::endl;
            terminaltextcolor(WHITE);
			return EXIT_FAILURE;
         }
         //std::cout << "ngp " << read_material[super_index].geometry << std::endl;
         for(int c=0;c<read_material[super_index].geometry;c++){
            for(int xy=0;xy<2;xy++){
               double var;
               gfile >> var;
               if(gfile.eof()){
				  terminaltextcolor(RED);
                  std::cerr << "Error in geometry input file " << value.c_str() << " end of file reached before reading all coordinates" << std::endl;
                  terminaltextcolor(WHITE);
				  return EXIT_FAILURE;
               }
               if((var<0.0) || (var > 1.0)){
				  terminaltextcolor(RED);
                  std::cerr << "Error in geometry input file " << value.c_str() << " value is outside of valid range (0.0-1.0)" << std::endl;
                  terminaltextcolor(WHITE);
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
      test="lattice-anisotropy-file";
      if(word==test){

         // Open file and check for success
         std::ifstream latt_file(value.c_str());
         if(!latt_file.is_open()){
            std::cerr << "Error: lattice file " << value.c_str() << " specified on line " << line << " of material file not found. Exiting." << std::endl;
            zlog << zTs() << "Error: lattice file " << value.c_str() << " specified on line " << line << " of material file not found. Exiting." << std::endl;
            return EXIT_FAILURE;
         }

         // specify number of points to be read
         int num_pts=0;

         // Read in number of temperature points
         latt_file >> num_pts;

         // Check for valid number of points
         if(num_pts<=1){
            std::cerr << "Error in lattice-anisotropy-file " << value.c_str() << " on line " << line << " of material file. The first number must be an integer greater than 1. Exiting." << std::endl;
            zlog << zTs() << "Error in lattice-anisotropy-file " << value.c_str() << " on line " << line << " of material file. The first number must be an integer greater than 1. Exiting." << std::endl;
            return EXIT_FAILURE;
         }

         // Loop over all lines in file
         for(int c=0;c<num_pts;c++){

            // temporary variables
            double T;
            double k;

            // Read in points and add them to material
            latt_file >> T >> k;

            // Check for premature end of file
            if(latt_file.eof()){
               std::cerr << "Error in lattice anisotropy-file " << value.c_str() << " on line " << line << " of material file. End of file reached before reading all values. Exiting" << std::endl;
               zlog << zTs() << "Error in lattice anisotropy-file " << value.c_str() << " on line " << line << " of material file. End of file reached before reading all values. Exiting" << std::endl;
               return EXIT_FAILURE;
            }
            read_material[super_index].lattice_anisotropy.add_point(T,k);
         }

         return EXIT_SUCCESS;

      }
      //--------------------------------------------------------------------
      else
      test="alloy-host"; // determines host material
      if(word==test){
         read_material[super_index].alloy_master=true; // if this keyword is set, then atoms of this type will be scanned for alloy materials
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      else
      test="alloy-class"; // determines unit cell category id for ordered alloys
      if(word==test){
         int ac=atoi(value.c_str());
         // test for 'disordered'
         std::string dis="disordered";
         if(value==dis){
            read_material[super_index].alloy_class=-1; // value for random alloy
            return EXIT_SUCCESS;
         }
         // test for valid ordered alloy, value of -1 will be deprecated
         if((ac<-1) || (ac > 3)){
			terminaltextcolor(RED);
            std::cerr << "Error in input file - material[" << super_index+1 << "]:alloy-class is outside of valid range (0-3)" << std::endl;
            terminaltextcolor(WHITE);
			return EXIT_FAILURE;
         }
         else{
            read_material[super_index].alloy_class=ac;
            return EXIT_SUCCESS;
         }
      }
      //--------------------------------------------------------------------
      else
      test="alloy-fraction"; // determines %mixing for disordered alloys
      if(word==test){
         double a=atof(value.c_str());
         if((a < 0.0) || (a > 1.0)){
			terminaltextcolor(RED);
            std::cerr << "Error in input file - material[" << super_index+1 << "]:alloy["<< sub_index+1 << "] is outside of valid range (0.0-1.0)" << std::endl;
            terminaltextcolor(WHITE);
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
      test="minimum-height";
      if(word==test){
         double min=atof(value.c_str());
         check_for_valid_value(min, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
         cs::SelectMaterialByZHeight=true;
         read_material[super_index].min=min;
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      else
      test="maximum-height";
      if(word==test){
         double max=atof(value.c_str());
         check_for_valid_value(max, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
         cs::SelectMaterialByZHeight=true;
         read_material[super_index].max=max;
         return EXIT_SUCCESS;
      }
      else
      //--------------------------------------------------------------------
      test="core-shell-size";
      if(word==test){
         double css=atof(value.c_str());
         check_for_valid_value(css, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
         read_material[super_index].core_shell_size=css;
         return EXIT_SUCCESS;
      }
      //-------------------------------------------------------------------
      else
      test="interface-roughness";
      if(word==test){
         double ir=atof(value.c_str());
         check_for_valid_value(ir, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
         read_material[super_index].interface_roughness=ir;
         return EXIT_SUCCESS;
      }
      else
      //-------------------------------------------------------------------
      test="density";
      if(word==test){
         double d=atof(value.c_str());
         check_for_valid_value(d, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
         read_material[super_index].density=d;
         return EXIT_SUCCESS;
      }
      else
      test="continuous";
      if(word==test){
         read_material[super_index].continuous=true;
         return EXIT_SUCCESS;
      }
      else
      //-------------------------------------------------------------------
      test="intermixing";
      if(word==test){
         double i=atof(value.c_str());
         //check_for_valid_value(i, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
         if((i<0.0) || (i > 1.0)){
			terminaltextcolor(RED);
            std::cerr << "Error in input file - material[" << super_index+1 << "]:intermixing[" << sub_index+1 <<"] is outside of valid range (0.0-1.0)" << std::endl;
			terminaltextcolor(WHITE);
            return EXIT_FAILURE;}
         else{
            read_material[super_index].intermixing[sub_index]=i;
            return EXIT_SUCCESS;
         }
      }
		//--------------------------------------------------------------------
		else
		test="constrained"; // determines use of alternate integrator
		if(word==test){
			string t="true";
			string f="false";
			if(value==t){
				read_material[super_index].constrained=true;
				return EXIT_SUCCESS;
			}
			else if(value==f){
				read_material[super_index].constrained=false;
				return EXIT_SUCCESS;
			}
			else {
				terminaltextcolor(RED);
				std::cerr << "Error in input file - material[" << super_index+1 << "]:constrained must be either true or false" << std::endl;
				terminaltextcolor(WHITE);
				return EXIT_FAILURE;
			}
		}
		//--------------------------------------------------------------------
		test="constraint-angle-theta";
		if(word==test){
			double angle=atof(value.c_str());
			// Test for valid range
			if((angle>=0.0) && (angle<=360.0)){
				cmc::cmc_mat[super_index].constraint_theta=angle;
				return EXIT_SUCCESS;
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 360.0" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		//--------------------------------------------------------------------
		test="constraint-angle-theta-minimum";
		if(word==test){
			double angle=atof(value.c_str());
			// Test for valid range
			if((angle>=0.0) && (angle<=360.0)){
				cmc::cmc_mat[super_index].constraint_theta_min=angle;
				return EXIT_SUCCESS;
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 360.0" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		//--------------------------------------------------------------------
		test="constraint-angle-theta-maximum";
		if(word==test){
			double angle=atof(value.c_str());
			// Test for valid range
			if((angle>=0.0) && (angle<=360.0)){
				cmc::cmc_mat[super_index].constraint_theta_max=angle;
				return EXIT_SUCCESS;
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 360.0" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		//--------------------------------------------------------------------
		test="constraint-angle-theta-increment";
		if(word==test){
			double angle=atof(value.c_str());
			// Test for valid range
			if((angle>=0.0) && (angle<=360.0)){
				cmc::cmc_mat[super_index].constraint_theta_delta=angle;
				return EXIT_SUCCESS;
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 360.0" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		//--------------------------------------------------------------------
		test="constraint-angle-phi-minimum";
		if(word==test){
			double angle=atof(value.c_str());
			// Test for valid range
			if((angle>=0.0) && (angle<=180.0)){
				cmc::cmc_mat[super_index].constraint_phi_min=angle;
				return EXIT_SUCCESS;
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 180.0" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		//--------------------------------------------------------------------
		test="constraint-angle-phi";
		if(word==test){
			double angle=atof(value.c_str());
			// Test for valid range
			if((angle>=0.0) && (angle<=180.0)){
				cmc::cmc_mat[super_index].constraint_phi=angle;
				return EXIT_SUCCESS;
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 180.0" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		//--------------------------------------------------------------------
		test="constraint-angle-phi-maximum";
		if(word==test){
			double angle=atof(value.c_str());
			// Test for valid range
			if((angle>=0.0) && (angle<=180.0)){
				cmc::cmc_mat[super_index].constraint_phi_max=angle;
				return EXIT_SUCCESS;
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 180.0" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		//--------------------------------------------------------------------
		test="constraint-angle-phi-increment";
		if(word==test){
			double angle=atof(value.c_str());
			// Test for valid range
			if((angle>=0.0) && (angle<=180.0)){
				cmc::cmc_mat[super_index].constraint_phi_delta=angle;
				return EXIT_SUCCESS;
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 180.0" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		//--------------------------------------------------------------------
		test="temperature";
		if(word==test){
			double T=atof(value.c_str());
			// Test for valid range
			if((T>=0.0) && (T<1.0E5)){
				read_material[super_index].temperature=T;
				// set local temperature flag
				sim::local_temperature=true;
				return EXIT_SUCCESS;
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 1.0E5" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
      //--------------------------------------------------------------------
      test="use-phonon-temperature";
      /*
        logical use-phonon-temperature
           This flag enables specific materials to couple to the phonon temperature
           of the system for simulations using the two temperature model. The default
           is for all materials to use the electron temperature. Valid values are true,
           false or (blank) [same as true].
       */
      if(word==test){
         // Test for sane input
         //bool sanitised_bool=check_for_valid_bool(value, word, line, prefix,"material");

         // set flag
         read_material[super_index].couple_to_phonon_temperature=true; //sanitised_bool;

         // enable local temperature flag
         sim::local_temperature=true;
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      test="fill-space";
      /*
        logical fill-space [false]
           This flag causes the material to be a fill material, instead of
           have atoms generated as part of the usual structure, for example
           particle shpaes, particle arrays and voronoi films. The default
           value is false for all materials. Valid values are true,
           false or (blank) [same as true].
       */
      if(word==test){
         // Test for sane input
         //bool sanitised_bool=check_for_valid_bool(value, word, line, prefix,"material");

         // set flag
         read_material[super_index].fill=true; //sanitised_bool;

         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
		test="applied-field-strength";
		if(word==test){
			double H=atof(value.c_str());
			// test for unit
			string unit_type="field";
			// if no unit given, assume internal
			if(unit.size() != 0){
				units::convert(unit,H,unit_type);
			}
			string str="field";
			if(unit_type==str){
				// Test for valid range
				if((H>=0.0) && (H<1.0E5)){
					read_material[super_index].applied_field_strength=H;
					// set local applied field flag
					sim::local_applied_field=true;
					return EXIT_SUCCESS;
				}
				else{
					terminaltextcolor(RED);
					std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 1.0E5" << std::endl;
					terminaltextcolor(WHITE);
					err::vexit();
				}
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - unit type \'" << unit_type << "\' is invalid for parameter material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 1.0E5" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		//------------------------------------------------------------
		test="applied-field-unit-vector";
		if(word==test){
			// temporary storage container
			std::vector<double> u(3);

			// read values from string
			u=DoublesFromString(value);

			// check size
			if(u.size()!=3){
				terminaltextcolor(RED);
				std::cerr << "Error in input file - material[" << super_index+1 << "]:"<< word << " must have three values." << std::endl;
				terminaltextcolor(WHITE);
				zlog << zTs() << "Error in input file - material[" << super_index+1 << "]:"<< word << " must have three values." << std::endl;
				return EXIT_FAILURE;
			}

			// Normalise
			double ULength=sqrt(u.at(0)*u.at(0)+u.at(1)*u.at(1)+u.at(2)*u.at(2));

			// Check for correct length unit vector
			if(ULength < 1.0e-9){
				terminaltextcolor(RED);
				std::cerr << "Error in input file - material[" << super_index+1 << "]:"<< word << " must be normalisable (possibly all zero)." << std::endl;
				terminaltextcolor(WHITE);
				zlog << zTs() << "Error in input file - material[" << super_index+1 << "]:"<< word << " must be normalisable (possibly all zero)." << std::endl;
				return EXIT_FAILURE;
			}
			u.at(0)/=ULength;
			u.at(1)/=ULength;
			u.at(2)/=ULength;

			// Copy anisotropy direction to material
			read_material[super_index].applied_field_unit_vector=u;

			// set local applied field flag
			sim::local_applied_field=true;

			return EXIT_SUCCESS;

		}
		//--------------------------------------------------------------------
		test="fmr-field-strength";
		if(word==test){
			double H=atof(value.c_str());
			// test for unit
			string unit_type="field";
			// if no unit given, assume internal
			if(unit.size() != 0){
				units::convert(unit,H,unit_type);
			}
			string str="field";
			if(unit_type==str){
				// Test for valid range
				if((H>=0.0) && (H<1.0E5)){
					read_material[super_index].applied_field_strength=H;
					// set local fmr flag
					sim::local_fmr_field=true;
					return EXIT_SUCCESS;
				}
				else{
					terminaltextcolor(RED);
					std::cerr << "Error - sim:" << word << " on line " << line << " of input file must be in the range 0 - 1.0E5" << std::endl;
					terminaltextcolor(WHITE);
					err::vexit();
				}
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - unit type \'" << unit_type << "\' is invalid for parameter material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 1.0E5" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		//--------------------------------------------------------------------
		test="fmr-field-frequency";
		if(word==test){
			double f=atof(value.c_str());
			// Test for valid range
			if((f>=0.0) && (f<1.0E20)){
				read_material[super_index].fmr_field_frequency=f;
				// set local fmr flag
				sim::local_fmr_field=true;
				return EXIT_SUCCESS;
			}
			else{
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 1.0E20" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		//------------------------------------------------------------
		test="fmr-field-unit-vector";
		if(word==test){
			// temporary storage container
			std::vector<double> u(3);

			// read values from string
			u=DoublesFromString(value);

			// check size
			if(u.size()!=3){
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " must have three values." << std::endl;
				terminaltextcolor(WHITE);
				zlog << zTs() << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " must have three values." << std::endl;
				return EXIT_FAILURE;
			}

			// Normalise
			double ULength=sqrt(u.at(0)*u.at(0)+u.at(1)*u.at(1)+u.at(2)*u.at(2));

			// Check for correct length unit vector
			if(ULength < 1.0e-9){
				terminaltextcolor(RED);
				std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " must be normalisable (possibly all zero)." << std::endl;
				terminaltextcolor(WHITE);
				zlog << zTs() << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " must be normalisable (possibly all zero)." << std::endl;
				return EXIT_FAILURE;
			}
			u.at(0)/=ULength;
			u.at(1)/=ULength;
			u.at(2)/=ULength;

			// Copy anisotropy direction to material
			read_material[super_index].fmr_field_unit_vector=u;

			// set local applied field flag
			sim::local_fmr_field=true;

			return EXIT_SUCCESS;

		}
      //--------------------------------------------------------------------
      else
      test="temperature-rescaling-exponent";
      if(word==test){
         double alpha=atof(value.c_str());
         check_for_valid_value(alpha, word, line, prefix, unit, "none", 0.0, 10.0,"material"," 0.0 - 10.0");
         read_material[super_index].temperature_rescaling_alpha=alpha;
         return EXIT_SUCCESS;
      }
      //--------------------------------------------------------------------
      else
      test="temperature-rescaling-curie-temperature";
      if(word==test){
         double Tc=atof(value.c_str());
         check_for_valid_value(Tc, word, line, prefix, unit, "none", 0.0, 10000.0,"material"," 0 - 10000 K");
         read_material[super_index].temperature_rescaling_Tc=Tc;
         return EXIT_SUCCESS;
      }
		//--------------------------------------------------------------------
		// keyword not found
		//--------------------------------------------------------------------
		else{
			terminaltextcolor(RED);
			std::cerr << "Error - Unknown control statement \'material[" << super_index+1 << "]:" << word << "\' on line " << line << " of material file" << std::endl;
			terminaltextcolor(WHITE);
			return EXIT_FAILURE;
		}

	return EXIT_SUCCESS;
}



} // end of namespace vin

namespace vout{

	// Namespace variable declarations
	std::vector<unsigned int> file_output_list(0);
	std::vector<unsigned int> screen_output_list(0);
	std::vector<unsigned int> grain_output_list(0);

   // Variables to control rate of data output to screen, output file and grain file
   int output_rate=1;
   int output_grain_rate=1;
   //int output_screen_rate=1; needs to be implemented

   bool gnuplot_array_format=false;

	std::ofstream errfile;

  #ifdef MPICF
	null_streambuf nullbuf;

	void redirect(std::ostream& strm, std::string filename) {
		errfile.open(filename.c_str());
		// redirect ouput into the file
		strm.rdbuf (errfile.rdbuf());
	}

	void nullify(std::ostream& strm){
		strm.rdbuf(&nullbuf);
	}
	#endif

	// Output Function 0
	void time(std::ostream& stream){
		stream << sim::time << "\t";
	}

	// Output Function 1
	void real_time(std::ostream& stream){
		stream << sim::time*mp::dt_SI << "\t";
	}

	// Output Function 2
	void temperature(std::ostream& stream){
		stream << sim::temperature << "\t";
	}

	// Output Function 3
	void Happ(std::ostream& stream){
		stream << sim::H_applied << "\t";
	}

	// Output Function 4
	void Hvec(std::ostream& stream){
		stream << sim::H_vec[0] << "\t"<< sim::H_vec[1] << "\t"<< sim::H_vec[2] << "\t";
	}

	// Output Function 5
   void mvec(std::ostream& stream){
      stream << stats::system_magnetization.output_normalized_magnetization();
   }

	// Output Function 6
   void magm(std::ostream& stream){
      stream << stats::system_magnetization.output_normalized_magnetization_length() << "\t";
   }

	// Output Function 7
   void mean_magm(std::ostream& stream){
      stream << stats::system_magnetization.output_normalized_mean_magnetization_length();
   }

	// Output Function 8
   void mat_mvec(std::ostream& stream){
      stream << stats::material_magnetization.output_normalized_magnetization();
   }

	// Output Function 9
   void mat_mean_magm(std::ostream& stream){
      stream << stats::material_magnetization.output_normalized_mean_magnetization_length();
   }

	// Output Function 10
	void grain_mvec(std::ostream& stream){

		unsigned int id=0; // grain id (excluding grains with zero atoms)

		// loop over all grains
		for(int grain=0;grain<grains::num_grains;grain++){
			// check for grains with zero atoms
			if(grains::grain_size_array[grain]!=0){
				stream << grains::x_mag_array[grain] << "\t";
				stream << grains::y_mag_array[grain] << "\t";
				stream << grains::z_mag_array[grain] << "\t";
				stream << grains::mag_m_array[grain] << "\t";
				id++;
			}
		}
	}

	// Output Function 11
	void grain_magm(std::ostream& stream){

		unsigned int id=0; // grain id (excluding grains with zero atoms)

		// loop over all grains
		for(int grain=0;grain<grains::num_grains;grain++){
			// check for grains with zero atoms
			if(grains::grain_size_array[grain]!=0){
				stream << grains::mag_m_array[grain] << "\t";
				id++;
			}
		}
	}

	// Output Function 12
	void mdoth(std::ostream& stream){
      // initialise vector of H
      std::vector<double> H(&sim::H_vec[0], &sim::H_vec[0]+3);
      stream << stats::system_magnetization.output_normalized_magnetization_dot_product(H);
	}

	// Output Function 13
	void grain_mat_mvec(std::ostream& stream){

		grains::output_mat_mag(stream);

	}

	// Output Function 14
	void systorque(std::ostream& stream){
		stream << stats::total_system_torque[0] << "\t";
		stream << stats::total_system_torque[1] << "\t";
		stream << stats::total_system_torque[2] << "\t";
	}

	// Output Function 15
	void mean_systorque(std::ostream& stream){
		stream << stats::total_mean_system_torque[0]/stats::torque_data_counter << "\t";
		stream << stats::total_mean_system_torque[1]/stats::torque_data_counter << "\t";
		stream << stats::total_mean_system_torque[2]/stats::torque_data_counter << "\t";
	}

	// Output Function 16
	void constraint_phi(std::ostream& stream){
		stream << sim::constraint_phi << "\t";
	}

	// Output Function 17
	void constraint_theta(std::ostream& stream){
		stream << sim::constraint_theta << "\t";
	}

	// Output Function 18
	void material_constraint_phi(std::ostream& stream){
		for(int mat=0;mat<mp::num_materials;mat++){
			stream << cmc::cmc_mat[mat].constraint_phi << "\t";
		}
	}

	// Output Function 19
	void material_constraint_theta(std::ostream& stream){
		for(int mat=0;mat<mp::num_materials;mat++){
			stream << cmc::cmc_mat[mat].constraint_theta << "\t";
		}
	}

	// Output Function 20
	void material_mean_systorque(std::ostream& stream){
		for(int mat=0;mat<mp::num_materials;mat++){
			stream << stats::sublattice_mean_torque_x_array[mat]/stats::torque_data_counter << "\t";
			stream << stats::sublattice_mean_torque_y_array[mat]/stats::torque_data_counter << "\t";
			stream << stats::sublattice_mean_torque_z_array[mat]/stats::torque_data_counter << "\t";
		}
	}

   // Output Function 21
   void mean_system_susceptibility(std::ostream& stream){
      stream << stats::system_susceptibility.output_mean_susceptibility(sim::temperature);
   }

	// Output Function 22
	void phonon_temperature(std::ostream& stream){
		stream << sim::TTTp << "\t";
	}

	// Output Function 23
	void material_temperature(std::ostream& stream){
		for(int mat=0;mat<mp::material.size();mat++){
			stream << mp::material[mat].temperature << "\t";
		}
	}

	// Output Function 24
	void material_applied_field_strength(std::ostream& stream){
		for(int mat=0;mat<mp::material.size();mat++){
			stream << mp::material[mat].applied_field_strength << "\t";
		}
	}

	// Output Function 25
	void material_fmr_field_strength(std::ostream& stream){
		const double real_time=sim::time*mp::dt_SI;

		for(int mat=0;mat<mp::material.size();mat++){
			const double Hsinwt_local=mp::material[mat].fmr_field_strength*sin(2.0*M_PI*real_time*mp::material[mat].fmr_field_frequency);
			stream << Hsinwt_local << "\t";
		}
	}

	// Output Function 26
	void mat_mdoth(std::ostream& stream){
      // initialise vector of H
      std::vector<double> H(&sim::H_vec[0], &sim::H_vec[0]+3);
      stream << stats::material_magnetization.output_normalized_magnetization_dot_product(H);
	}

   // Output Function 27
   void total_energy(std::ostream& stream){
      stats::output_energy(stream, stats::all, stats::total);
   }

   // Output Function 28
   void mean_total_energy(std::ostream& stream){
      stats::output_energy(stream, stats::all, stats::mean);
   }

   // Output Function 29
   void total_anisotropy_energy(std::ostream& stream){
      stats::output_energy(stream, stats::anisotropy, stats::total);
   }

   // Output Function 30
   void mean_total_anisotropy_energy(std::ostream& stream){
      stats::output_energy(stream, stats::anisotropy, stats::mean);
   }

   // Output Function 31
   void total_cubic_anisotropy_energy(std::ostream& stream){
      stats::output_energy(stream, stats::cubic_anisotropy, stats::total);
   }

   // Output Function 32
   void mean_total_cubic_anisotropy_energy(std::ostream& stream){
      stats::output_energy(stream, stats::cubic_anisotropy, stats::mean);
   }

   // Output Function 33
   void total_surface_anisotropy_energy(std::ostream& stream){
      stats::output_energy(stream, stats::surface_anisotropy, stats::total);
   }

   // Output Function 34
   void mean_total_surface_anisotropy_energy(std::ostream& stream){
      stats::output_energy(stream, stats::surface_anisotropy, stats::mean);
   }

   // Output Function 35
   void total_exchange_energy(std::ostream& stream){
      stats::output_energy(stream, stats::exchange, stats::total);
   }

   // Output Function 36
   void mean_total_exchange_energy(std::ostream& stream){
      stats::output_energy(stream, stats::exchange, stats::mean);
   }

   // Output Function 37
   void total_applied_field_energy(std::ostream& stream){
      stats::output_energy(stream, stats::applied_field, stats::total);
   }

   // Output Function 38
   void mean_total_applied_field_energy(std::ostream& stream){
      stats::output_energy(stream, stats::applied_field, stats::mean);
   }

   // Output Function 39
   void total_magnetostatic_energy(std::ostream& stream){
      stats::output_energy(stream, stats::magnetostatic, stats::total);
   }

   // Output Function 40
   void mean_total_magnetostatic_energy(std::ostream& stream){
      stats::output_energy(stream, stats::magnetostatic, stats::mean);
   }

   // Output Function 41
   void total_so_anisotropy_energy(std::ostream& stream){
      stats::output_energy(stream, stats::second_order_anisotropy, stats::total);
   }

   // Output Function 42
   void mean_total_so_anisotropy_energy(std::ostream& stream){
      stats::output_energy(stream, stats::second_order_anisotropy, stats::mean);
   }

   // Output Function 43
   void height_mvec(std::ostream& stream){
      stream << stats::height_magnetization.output_normalized_magnetization();
   }

   // Output Function 44
   void material_height_mvec(std::ostream& stream){
      stream << stats::material_height_magnetization.output_normalized_magnetization();
   }

   // Output Function 45
   void height_mvec_actual(std::ostream& stream){
      stream << stats::height_magnetization.output_magnetization();
   }

   // Output Function 46
   void material_height_mvec_actual(std::ostream& stream){
      stream << stats::material_height_magnetization.output_magnetization();
   }

   // Output Function 60
	void MPITimings(std::ostream& stream){

		stream << vmpi::AverageComputeTime+vmpi::AverageWaitTime << "\t" << vmpi::AverageComputeTime << "\t" << vmpi::AverageWaitTime;
		stream << "\t" << vmpi::MaximumComputeTime << "\t" << vmpi::MaximumWaitTime << "\t";
	}

	// Data output wrapper function
	void data(){

		// check calling of routine if error checking is activated
		if(err::check==true){std::cout << "vout::data has been called" << std::endl;}

		// Calculate MPI Timings since last data output
		#ifdef MPICF
		if(vmpi::DetailedMPITiming){

			// Calculate Average times
			MPI_Reduce (&vmpi::TotalComputeTime,&vmpi::AverageComputeTime,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce (&vmpi::TotalWaitTime,&vmpi::AverageWaitTime,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			vmpi::AverageComputeTime/=double(vmpi::num_processors);
			vmpi::AverageWaitTime/=double(vmpi::num_processors);

			// Calculate Maximum times
			MPI_Reduce (&vmpi::TotalComputeTime,&vmpi::MaximumComputeTime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
			MPI_Reduce (&vmpi::TotalWaitTime,&vmpi::MaximumWaitTime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

			// Save times for timing matrix
			vmpi::ComputeTimeArray.push_back(vmpi::TotalComputeTime);
			vmpi::WaitTimeArray.push_back(vmpi::TotalWaitTime);

			// reset until next data output
			vmpi::TotalComputeTime=0.0;
			vmpi::TotalWaitTime=0.0;
		}
		#endif

      // check for open ofstream
      if(!zmag.is_open()){
         // check for checkpoint continue and append data
         if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag) zmag.open("output",std::ofstream::app);
         // otherwise overwrite file
         else{
            zmag.open("output",std::ofstream::trunc);
            // write file header information
            if(vmpi::my_rank==0) write_output_file_header(zmag, file_output_list);
         }
      }

		// Only output 1/output_rate time steps
      if(sim::time%vout::output_rate==0){

		// Output data to output
      if(vmpi::my_rank==0){

      // For gpu acceleration get statistics from device
      if(gpu::acceleration) gpu::stats::get();

		for(unsigned int item=0;item<file_output_list.size();item++){
			switch(file_output_list[item]){
				case 0:
					vout::time(zmag);
					break;
				case 1:
					vout::real_time(zmag);
					break;
				case 2:
					vout::temperature(zmag);
					break;
				case 3:
					vout::Happ(zmag);
					break;
				case 4:
					vout::Hvec(zmag);
					break;
				case 5:
					vout::mvec(zmag);
					break;
				case 6:
					vout::magm(zmag);
					break;
				case 7:
					vout::mean_magm(zmag);
					break;
				case 8:
					vout::mat_mvec(zmag);
					break;
				case 9:
					vout::mat_mean_magm(zmag);
					break;
				case 12:
					vout::mdoth(zmag);
					break;
				case 14:
					vout::systorque(zmag);
					break;
				case 15:
					vout::mean_systorque(zmag);
					break;
				case 16:
					vout::constraint_phi(zmag);
					break;
				case 17:
					vout::constraint_theta(zmag);
					break;
				case 18:
					vout::material_constraint_phi(zmag);
					break;
				case 19:
					vout::material_constraint_theta(zmag);
					break;
				case 20:
					vout::material_mean_systorque(zmag);
					break;
				case 21:
					vout::mean_system_susceptibility(zmag);
					break;
				case 22:
					vout::phonon_temperature(zmag);
					break;
				case 23:
					vout::material_temperature(zmag);
					break;
				case 24:
					vout::material_applied_field_strength(zmag);
					break;
				case 25:
					vout::material_fmr_field_strength(zmag);
					break;
				case 26:
					vout::mat_mdoth(zmag);
					break;
            case 27:
               vout::total_energy(zmag);
               break;
            case 28:
               vout::mean_total_energy(zmag);
               break;
            case 29:
               vout::total_anisotropy_energy(zmag);
               break;
            case 30:
               vout::mean_total_anisotropy_energy(zmag);
               break;
            case 31:
               vout::total_cubic_anisotropy_energy(zmag);
               break;
            case 32:
               vout::mean_total_cubic_anisotropy_energy(zmag);
               break;
            case 33:
               vout::total_surface_anisotropy_energy(zmag);
               break;
            case 34:
               vout::mean_total_surface_anisotropy_energy(zmag);
               break;
            case 35:
               vout::total_exchange_energy(zmag);
               break;
            case 36:
               vout::mean_total_exchange_energy(zmag);
               break;
            case 37:
               vout::total_applied_field_energy(zmag);
               break;
            case 38:
               vout::mean_total_applied_field_energy(zmag);
               break;
            case 39:
               vout::total_magnetostatic_energy(zmag);
               break;
            case 40:
               vout::mean_total_magnetostatic_energy(zmag);
               break;
            case 41:
               vout::total_so_anisotropy_energy(zmag);
               break;
            case 42:
               vout::mean_total_so_anisotropy_energy(zmag);
               break;
            case 43:
               vout::height_mvec(zmag);
               break;
            case 44:
               vout::material_height_mvec(zmag);
               break;
            case 45:
               vout::height_mvec_actual(zmag);
               break;
            case 46:
               vout::material_height_mvec_actual(zmag);
               break;
            case 60:
					vout::MPITimings(zmag);
					break;
			}
		}
		// Carriage return
		if(file_output_list.size()>0) zmag << std::endl;

      } // end of code for rank 0 only
   } // end of if statement for output rate

		// Output data to cout
		if(vmpi::my_rank==0){
      if(sim::time%vout::output_rate==0){ // needs to be altered to separate variable at some point

         for(unsigned int item=0;item<screen_output_list.size();item++){
			switch(screen_output_list[item]){
				case 0:
					vout::time(std::cout);
					break;
				case 1:
					vout::real_time(std::cout);
					break;
				case 2:
					vout::temperature(std::cout);
					break;
				case 3:
					vout::Happ(std::cout);
					break;
				case 4:
					vout::Hvec(std::cout);
					break;
				case 5:
					vout::mvec(std::cout);
					break;
				case 6:
					vout::magm(std::cout);
					break;
				case 7:
					vout::mean_magm(std::cout);
					break;
				case 8:
					vout::mat_mvec(std::cout);
					break;
				case 9:
					vout::mat_mean_magm(std::cout);
					break;
				case 12:
					vout::mdoth(std::cout);
					break;
				case 14:
					vout::systorque(std::cout);
					break;
				case 15:
					vout::mean_systorque(std::cout);
					break;
				case 16:
					vout::constraint_phi(std::cout);
					break;
				case 17:
					vout::constraint_theta(std::cout);
					break;
				case 18:
					vout::material_constraint_phi(std::cout);
					break;
				case 19:
					vout::material_constraint_theta(std::cout);
					break;
				case 20:
					vout::material_mean_systorque(std::cout);
					break;
				case 21:
					vout::mean_system_susceptibility(std::cout);
					break;
				case 22:
					vout::phonon_temperature(std::cout);
					break;
				case 23:
					vout::material_temperature(std::cout);
					break;
				case 24:
					vout::material_applied_field_strength(std::cout);
					break;
				case 25:
					vout::material_fmr_field_strength(std::cout);
					break;
				case 26:
					vout::mat_mdoth(std::cout);
					break;
            case 27:
               vout::total_energy(std::cout);
               break;
            case 28:
               vout::mean_total_energy(std::cout);
               break;
            case 29:
               vout::total_anisotropy_energy(std::cout);
               break;
            case 30:
               vout::mean_total_anisotropy_energy(std::cout);
               break;
            case 31:
               vout::total_cubic_anisotropy_energy(std::cout);
               break;
            case 32:
               vout::mean_total_cubic_anisotropy_energy(std::cout);
               break;
            case 33:
               vout::total_surface_anisotropy_energy(std::cout);
               break;
            case 34:
               vout::mean_total_surface_anisotropy_energy(std::cout);
               break;
            case 35:
               vout::total_exchange_energy(std::cout);
               break;
            case 36:
               vout::mean_total_exchange_energy(std::cout);
               break;
            case 37:
               vout::total_applied_field_energy(std::cout);
               break;
            case 38:
               vout::mean_total_applied_field_energy(std::cout);
               break;
            case 39:
               vout::total_magnetostatic_energy(std::cout);
               break;
            case 40:
               vout::mean_total_magnetostatic_energy(std::cout);
               break;
            case 41:
               vout::total_so_anisotropy_energy(std::cout);
               break;
            case 42:
               vout::mean_total_so_anisotropy_energy(std::cout);
               break;
            case 60:
					vout::MPITimings(std::cout);
					break;
			}
		}

		// Carriage return
		if(screen_output_list.size()>0) std::cout << std::endl;
		}

   } // End of if statement to output data to screen

		if(sim::time%vout::output_grain_rate==0){

		// calculate grain magnetisations
		grains::mag();

		// Output data to zgrain
		if(vmpi::my_rank==0){

			// check for open ofstream
         if(vout::grain_output_list.size() > 0 && !zgrain.is_open()){
            // check for checkpoint continue and append data
            if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag) zgrain.open("grain",std::ofstream::app);
            // otherwise overwrite file
            else zgrain.open("grain",std::ofstream::trunc);
         }

			for(unsigned int item=0;item<vout::grain_output_list.size();item++){
			switch(vout::grain_output_list[item]){
				case 0:
					vout::time(zgrain);
					break;
				case 1:
					vout::real_time(zgrain);
					break;
				case 2:
					vout::temperature(zgrain);
					break;
				case 3:
					vout::Happ(zgrain);
					break;
				case 4:
					vout::Hvec(zgrain);
					break;
				case 10:
					vout::grain_mvec(zgrain);
					break;
				case 11:
					vout::grain_magm(zgrain);
					break;
				case 13:
					vout::grain_mat_mvec(zgrain);
					break;
			   case 22:
					vout::phonon_temperature(zgrain);
               break;
			}
		}

		// Carriage return
		if(vout::grain_output_list.size()>0) zgrain << std::endl;
		}
		}

		vout::config();

      // optionally save checkpoint file
      if(sim::save_checkpoint_flag==true && sim::save_checkpoint_continuous_flag==true && sim::time%sim::save_checkpoint_rate==0) save_checkpoint();

	} // end of data

} // end of namespace vout
