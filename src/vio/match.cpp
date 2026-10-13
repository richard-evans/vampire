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
#include <iostream>

// Vampire headers
#include "anisotropy.hpp"
#include "cells.hpp"
#include "config.hpp"
#include "create.hpp"
#include "dipole.hpp"
#include "environment.hpp"
#include "errors.hpp"
#include "exchange.hpp"
#include "gpu.hpp"
#include "hamr.hpp"
#include "hierarchical.hpp"
#include "ltmp.hpp"
#include "material.hpp"
#include "micromagnetic.hpp"
#include "montecarlo.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "spininitialize.hpp"
#include "spintextures.hpp"
#include "spintorque.hpp"
#include "spintransport.hpp"
#include "spinwaves.hpp" // JRH
#include "stats.hpp"
#include "unitcell.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// vio module headers
#include "internal.hpp"

namespace vin{

   //---------------------------------------------------------------------------
   // Function to match input file keywords to module input parameters. Each
   // module processes its own keywords in src/<module>/interface.cpp
   //---------------------------------------------------------------------------
   int match(string const key, string const word, string const value, string const unit, int const line){

      //-------------------------------------------------------------------
      // Call module input parameters
      //-------------------------------------------------------------------
      if(ltmp::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(anisotropy::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(cells::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(config::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(create::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(dipole::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(gpu::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(exchange::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(hierarchical::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(montecarlo::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(program::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(sim::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(stats::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(st::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(spin_transport::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(unitcell::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(vio::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(vmpi::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(micromagnetic::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(environment::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(hamr::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(sld::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(spininitialize::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
      else if(spinwaves::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS; // JRH spinwaves input parameters
      else if(spintextures::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;

      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      terminaltextcolor(RED);
      std::cerr << "Error - Unknown control statement \'" << key <<":"<< word << "\' on line " << line << " of input file" << std::endl;
      terminaltextcolor(WHITE);
      return EXIT_FAILURE;

   } // end of match function

   int match(string const key, string const word, string const value, string const unit, int const line, int const superIndex, int const subIndex){
      if(unitcell::match_input_parameter(key, word, value, unit, line, superIndex, subIndex)) return EXIT_SUCCESS;

      terminaltextcolor(RED);
      std::cerr << "Error - Unknown control statement \'" << key <<":"<< word << "\' on line " << line << " of input file" << std::endl;
      terminaltextcolor(WHITE);
      return EXIT_FAILURE;
   }

    // temporary array of materials for reading in material data
    std::vector<mp::materials_t> read_material(0);

    int read_mat_file(std::string const matfile, int const line_number){


        // resize temporary materials array for storage of variables
        read_material.resize(mp::max_materials);
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
            int matchcheck = vin::match_material(word, value, unit, line_counter, super_index-1, sub_index-1, original_line, matfile);
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
        //inputfile.close();

        return EXIT_SUCCESS;

    }

   //---------------------------------------------------------------------------
   // Function to match material file keywords to module material parameters
   //---------------------------------------------------------------------------
   int match_material(string const word,
                      string const value,
                      string const unit,
                      int const line,
                      int const super_index,
                      int const sub_index,
                      std::string const line_string,
                      std::string const filename_string)
   {

      //-------------------------------------------------------------------
      // Call module material parameters
      //-------------------------------------------------------------------
      if(vio::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(anisotropy::match_material_parameter(word, value, unit, line, super_index, sub_index, mp::max_materials)) return EXIT_SUCCESS;
      else if(cells::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(config::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(create::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(dipole::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(exchange::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(hierarchical::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(montecarlo::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(program::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(sim::match_material_parameter(word, value, unit, line, super_index)) return EXIT_SUCCESS;
      else if(stats::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(st::match_material(word, value, unit, line, super_index)) return EXIT_SUCCESS;
      else if(spin_transport::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(unitcell::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(micromagnetic::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(environment::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(sld::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(spininitialize::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(spinwaves::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
      else if(spintextures::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;

      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      terminaltextcolor(RED);
      std::cerr << "Error - Unknown control statement '" << line_string << "' on line " << line << " of material file '" << filename_string << "'" << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error - Unknown control statement '" << line_string << " on line " << line << " of material file '" << filename_string << "'" << std::endl;
      return EXIT_FAILURE;

   }

}
