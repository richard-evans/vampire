//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and Rory Pond 2016. (Jack Collings 2021) All rights reserved.
//
//   Email: richard.evans@york.ac.uk and rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <sstream>
#include <cstring>
#include <string>
#include <iostream>
// Vampire headers
#include "vio.hpp"
#include "errors.hpp"

// vio module headers
#include "internal.hpp"

namespace vin{
	// Function to extract all variables from a string and return a vector of double
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

	// Function to extract all variables from a string and return a vector of int
	std::vector<int> integers_from_string(std::string value){

		// array for storing variables
		std::vector<int> array(0);

		// set source for ss
		std::istringstream source(value);

		// int variable to store values
		int temp = 0.0;

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

      // Comment and delimiter characters for "input" type files.
      char delim[] = ":=!";
		char delim2[] = "[]";

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
			std::string tmpWord="";
			std::string word="";
			std::string superIndexString = "";
			std::string subIndexString = "";
			std::string value="";
			std::string unit="";

			// get size of string

         // remove everything after comment character
         line = line.substr(0,line.find('#'));

			//for (int i = 0; i < line.length(); ++i){
			//	int start = 0;
			//	if (line.at(i) == ':'){
			//		key = line.substr(0, i + 1 - start);
			//		start = i;
			//	}
			//	else if(line.at(i) == '[' && word == "" ){
			//		word = line.substr(start, i);
			//	}
			//	else if line.at()
			//}

         // convert to c-string style, for tokenisation
         std::vector<char> cstr(line.begin(),line.end());

         // add null terminator, to insure that strtok cannot run over into previously used memory.
         cstr.push_back('\0');

         // tokenise the string, using delimiters from above
         char *token = strtok(&cstr[0],delim); // first call of strtok sets the string to tokenise.
         for (int count = 0; count < 4 && token !=NULL; count++){
            if (count==0){key=token;}       // Format is always the same
            else if(count==1){tmpWord=token;}  // but breaks if EOL found
            else if(count==2){value=token;} // so if unused, keywords will remain as ""
            else if(count==3){unit=token;}
            token = strtok(NULL,delim);
         };

			// Break up tmpWord into the word, superIndex and subIndex (similar to above)
			std::vector<char> wordcstr(tmpWord.begin(), tmpWord.end());
			wordcstr.push_back('\0');
			char *wordToken = strtok(&wordcstr[0],delim2);
			for (int count = 0; count < 3 && wordToken != NULL; ++count){
				if (count == 0){word=wordToken;}
				else if (count == 1){superIndexString=wordToken;}
				else if (count == 2){subIndexString=wordToken;}
				wordToken = strtok(NULL, delim2);
			}

         // Change super/sub indicies from string to int form
			string empty="";
	
			int superIndex = 0;
			int subIndex = 0;

			if(subIndexString != empty){
				subIndex = stoi(subIndexString);
				superIndex = stoi(superIndexString);
			}
			else if(superIndexString != empty){
				superIndex = stoi(superIndexString);
			}

			// Call different overloads depending on whether super and sub indicies are present
			if(key != empty && superIndex == 0 && subIndex == 0){
				//	std::cout << "\t" << "key:  " << key << std::endl;
				//	std::cout << "\t" << "word: " << word << std::endl;
				//	std::cout << "\t" << "value:" << value << std::endl;
				// std::cout << "\t" << "unit: " << unit << std::endl;
				int matchcheck = match(key, word, value, unit, line_counter);
				if(matchcheck==EXIT_FAILURE){
					err::vexit();
				}
			}
			else if (key != empty && superIndex != 0 && subIndex != 0){
				//std::cout << "\t" << "key:				" << key << std::endl;
				//std::cout << "\t" << "word: 			" << word << std::endl;
				//std::cout << "\t" << "value:			" << value << std::endl;
				//std::cout << "\t" << "unit:			" << unit << std::endl;
				//std::cout << "\t" << "superIndex:	" << superIndex << std::endl;
				//std::cout << "\t" << "subIndex:		" << subIndex << std::endl;
				int matchcheck = match(key, word, value, unit, line_counter, superIndex - 1, subIndex - 1);
				if(matchcheck == EXIT_FAILURE){
					err::vexit();
				}
			}
			else if (key != empty && superIndex != 0){
				err::vexit();	// There are no functionalities in vampire which use only a superIndex in the input file.
			}
		}

		return EXIT_SUCCESS;
	}


}
