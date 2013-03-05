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

namespace vio{
   
   bool vio_help=false;
//--------------------------------------------------------------------------
//
//  Function to match set of category, keyword, value and unit.
//  
//  (c) R F L Evans 2013
//
//--------------------------------------------------------------------------
//
int match_keyword(string const category, string const keyword, string const value, string const unit, int const line){

   // temporary string for comparison
   std::string test;

   // Test for create variables
   //----------------------------------
   //test="create";
   //if(category==test){
      //int frs=vin::match_create(keyword, value, unit, line);
      //return frs;
   //}
   // Test for dimension variables
   //----------------------------------
   //else
   test="dimensions";
   if(category==test){
      int frs=vio::match_dimension(category, keyword, value, unit, line);
      return frs;
   }
   // Test for simulation variables
   //----------------------------------
   /*else
   test="sim";
   if(category==test){
      int frs=vin::match_sim(keyword, value, unit, line);
      return frs;
   }
   // Test for data file output
   //----------------------------------
   else
   test="output";
   if(category==test){
      int frs=vin::match_vout_list(keyword, line, vout::file_output_list);
      return frs;
   }
   // Test for screen output
   //----------------------------------
   else
   test="screen";
   if(category==test){
      int frs=vin::match_vout_list(keyword, line, vout::screen_output_list);
      return frs;
   }
   // Test for grain output
   //----------------------------------
   else
   test="grain";
   if(category==test){
      int frs=vin::match_vout_grain_list(keyword, value, line, vout::grain_output_list);
      return frs;
   }  
   // Test for config output
   //----------------------------------
   else
   test="config";
   if(category==test){
      int frs=vin::match_config(keyword, value, line);
      return frs;
   }  
   // Get material filename
   //-------------------------------------------------------------------
   else
   test="material";
   if(category==test){
      test="file";
      if(keyword==test){
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
            std::cerr << "Error - empty filename in control statement \'material:" << keyword << "\' on line " << line << " of input file" << std::endl;
            return EXIT_FAILURE;
         }
      }
      // Get unit cell filename
      //-------------------------------------------------------------------
      test="unit-cell-file";
      if(keyword==test){
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
            std::cerr << "Error - empty filename in control statement \'material:" << keyword << "\' on line " << line << " of input file" << std::endl;
            return EXIT_FAILURE;
         }
      }
      else{
         std::cerr << "Error - Unknown control statement \'material:" << keyword << "\' on line " << line << " of input file" << std::endl;
         return EXIT_FAILURE;
      }
   }*/
   else
      std::cerr << "Error - Unknown control statement \'" << category <<":"<< keyword << "\' on line " << line << " of input file" << std::endl;
      return EXIT_FAILURE;

} // end of match function


int keyword_help(string const keyword){

   // set help flag to true
   vio_help=true;
   
   string const category="";
   string const value=""; 
   string const unit=""; 
   int const line=0;
   
   // Initialise input parameter database
   vio::initialise_dimensions();
   


   // Test for create variables
   //----------------------------------
   //test="create";
   //if(category==test){
      //int frs=vin::match_create(keyword, value, unit, line);
      //return frs;
   //}
   // Test for dimension variables
   //----------------------------------
   //else
      int frs=vio::match_dimension(category, keyword, value, unit, line);
      if(frs==EXIT_SUCCESS) return frs;
      
   // Test for simulation variables
   //----------------------------------
   /*else
   test="sim";
   if(category==test){
      int frs=vin::match_sim(keyword, value, unit, line);
      return frs;
   }
   // Test for data file output
   //----------------------------------
   else
   test="output";
   if(category==test){
      int frs=vin::match_vout_list(keyword, line, vout::file_output_list);
      return frs;
   }
   // Test for screen output
   //----------------------------------
   else
   test="screen";
   if(category==test){
      int frs=vin::match_vout_list(keyword, line, vout::screen_output_list);
      return frs;
   }
   // Test for grain output
   //----------------------------------
   else
   test="grain";
   if(category==test){
      int frs=vin::match_vout_grain_list(keyword, value, line, vout::grain_output_list);
      return frs;
   }  
   // Test for config output
   //----------------------------------
   else
   test="config";
   if(category==test){
      int frs=vin::match_config(keyword, value, line);
      return frs;
   }  
   // Get material filename
   //-------------------------------------------------------------------
   else
   test="material";
   if(category==test){
      test="file";
      if(keyword==test){
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
            std::cerr << "Error - empty filename in control statement \'material:" << keyword << "\' on line " << line << " of input file" << std::endl;
            return EXIT_FAILURE;
         }
      }
      // Get unit cell filename
      //-------------------------------------------------------------------
      test="unit-cell-file";
      if(keyword==test){
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
            std::cerr << "Error - empty filename in control statement \'material:" << keyword << "\' on line " << line << " of input file" << std::endl;
            return EXIT_FAILURE;
         }
      }
      else{
         std::cerr << "Error - Unknown control statement \'material:" << keyword << "\' on line " << line << " of input file" << std::endl;
         return EXIT_FAILURE;
      }
   }*/
   
      std::cerr << "Keyword \'" << keyword << "\' not found." << std::endl;
      return EXIT_FAILURE;

} // end of match function



} // end of namespace vio