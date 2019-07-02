//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and Rory Pond 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk and rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <sstream>

// Vampire headers
#include "vio.hpp"
#include "errors.hpp"

// vio module headers
#include "internal.hpp"

namespace vin{
	// Function to extract all variables from a string and return a vector
	std::vector<double> doubles_from_string(std::string value){

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
		//zlog << zTs() << "READING" << std::endl;
		// Print informative message to zlog file
		zlog << zTs() << "Opening main input file \"" << filename << "\"." << std::endl;

		std::stringstream inputfile;
		inputfile.str( vin::get_string(filename.c_str(), "input", -1) );

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
	//		std::cout << "\t" << "key:  " << key << std::endl;
//			std::cout << "\t" << "word: " << word << std::endl;
//			std::cout << "\t" << "value:" << value << std::endl;
//		   std::cout << "\t" << "unit: " << unit << std::endl;
			int matchcheck = match(key, word, value, unit, line_counter);
			if(matchcheck==EXIT_FAILURE){
				err::vexit();
			}
			}
		}

		return EXIT_SUCCESS;
	}


}
