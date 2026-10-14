//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <sstream>
#include <iostream>

// Vampire headers
#include "errors.hpp"
#include "material.hpp"
#include "montecarlo.hpp"
#include "vio.hpp"

// material module headers
#include "internal.hpp"

namespace mp{

    // temporary array of materials for reading in material data
    std::vector<materials_t> read_material(0);

namespace internal{

    int read_mat_file(std::string const matfile, int const line_number){


        // resize temporary materials array for storage of variables
        mp::read_material.resize(mp::max_materials);
        montecarlo::cmc::cmc_mat.resize(mp::max_materials);

            // Print informative message to zlog file
        zlog << zTs() << "Opening material file \"" << matfile << "\"." << std::endl;

        // Open file read only
        std::stringstream inputfile;
        inputfile.str( vin::get_string(matfile.c_str(), "material", line_number) );
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

            // save a copy of the line before stripping characters in case of error
            std::string original_line = line;

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
            const char* colon=":";  // Word identifier
            const char* eq="=";     // Value identifier
            const char* exc="!";    // Unit identifier
            const char* hash="#";   // Comment identifier
            const char* si="[";     // Index identifier
            const char* ei="]";     // End index identifier

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
                        std::ostringstream message;
                        message << "Causes could be invalid character or outside of range, ie less than 1 or greater than max_materials=" << mp::max_materials << ", exiting";
                        err::fatal(message.str());
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
                                    std::ostringstream message;
                                    message << "Causes could be invalid character or outside of range, ie less than 1 or greater than max_materials=" << mp::max_materials << ", exiting";
                                    err::fatal(message.str());
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
            std::string empty="";
            if(key!=empty){
                //std::cout << key << "[" << super_index << "]:" << word << "[" << sub_index << "]=" << value << " !" << unit << std::endl;
                //std::cout << "\t" << "key:  " << key << std::endl;
                //std::cout << "\t" << "word: " << word << std::endl;
                //std::cout << "\t" << "value:" << value << std::endl;
                //std::cout << "\t" << "unit: " << unit << std::endl;
            int matchcheck = match_material(word, value, unit, line_counter, super_index-1, sub_index-1, original_line, matfile);
                if(matchcheck==EXIT_FAILURE){
                    err::vexit();
                }
            }
        }

        // resize global material array
        mp::material.resize(mp::num_materials);

        // Copy data to global material array
        for(int mat=0;mat<mp::num_materials;mat++){
            mp::material[mat]=mp::read_material[mat];
        }

        // Resize read array to zero
        mp::read_material.resize(0);


        // Close file
        //inputfile.close();

        return EXIT_SUCCESS;

    }

} // end of internal namespace

} // end of mp namespace
