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
//  Revised parser for vampire input, material and ucf files
//
//  (c) R F L Evans 2013
//
//-----------------------------------------------------------------------------
//
// c++ headers
#include <fstream>
#include <iostream>
#include <string>

// vampire headers
#include "errors.hpp"
#include "vio.hpp"

// local vampire headers
#include "vio.ihpp"

namespace vio{
//--------------------------------------------------------------------------
//
//  Function to parse input file for keywords and set simulation parameters.
//  
//  (c) R F L Evans 2013
//
//--------------------------------------------------------------------------
//
void parse_vampire_input_file(std::string const filename){

   using vio::vTs;
   using vio::vlog;
   
   // Initialise input parameter checking
   vio::initialise_dimensions();
   
   // ifstream declaration
   std::ifstream inputfile;

   // Print informative message to vlog file
   vlog << vTs() << "Opening main input file \"" << filename << "\"." << std::endl; 

   // Open file read only
   inputfile.open(filename.c_str());

   // Check for opening
   if(!inputfile.is_open()){
      std::cerr << "Error opening main input file \"" << filename << "\". File does not exist!" << std::endl;
      vlog << vTs() << "Error: Main input file \"" << filename << "\" cannot be opened or does not exist." << std::endl;
      vlog << vTs() << "If file exists then check file permissions to ensure it is readable by the user." << std::endl;
      err::vexit();
   }

   // Print informative message to vlog file
   vlog << zTs() << "Parsing system parameters from main input file." << std::endl;

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
      std::string category="";
      std::string keyword="";
      std::string value="";
      std::string unit="";

      // get size of string
      int linelength = line.length();
      int last=0;
      
      // set character triggers
      const char* colon=":";  // Word identifier
      const char* eq="=";     // Value identifier
      const char* exc="!";    // Unit identifier
      const char* hash="#";   // Comment identifier
      
      // Determine key by looping over characters in line
      for(int i=0;i<linelength;i++){
         char c=line.at(i);
         last=i;
         
         // if character is not ":" or "=" or "!" or "#" interpret as key
         if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
            category.push_back(c);
         }
         else break;
      }
      const int end_key=last;
      
      // Determine the rest
      for(int i=end_key;i<linelength;i++){
         
         char c=line.at(i);

         // period found - interpret as word
         if(c== *colon){
            for(int j=i+1;j<linelength;j++){
               // if character is not special add to value
               char c=line.at(j);
               if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
                  keyword.push_back(c);
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
      }
      string empty="";
      if(category!=empty){
         // Check for match in list of keywords
         int matchcheck = vio::match_keyword(category, keyword, value, unit, line_counter);
         // Check for valid match
         if(matchcheck==EXIT_FAILURE) err::vexit();
      }
   }
   // Close file
   inputfile.close();

   // temporarily exit
   err::vexit();
   
   return;
}

} // end of namespace vio
