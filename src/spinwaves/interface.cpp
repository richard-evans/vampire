//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sergiu Ruta 2022. All rights reserved.
//
//   Email: sergiu.ruta@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <iostream>
#include <string>

// Vampire headers
#include "material.hpp"
#include "spinwaves.hpp"
#include "errors.hpp"
#include "vio.hpp"

// sw module headers
#include "internal.hpp"

namespace spinwaves{

   int largest = 0;

   //---------------------------------------------------------------------------
   // Function to process input file parameters for sw module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // std::cout << key << std::endl;
      // Check for the generic spinwaves key used for input file and reduction method
      std::string prefix="spinwaves";
      if(key==prefix){

         std::string test="kfile";
         if(word==test){
            std::string kpath_file=value;
            // strip quotes
            kpath_file.erase(std::remove(kpath_file.begin(), kpath_file.end(), '\"'), kpath_file.end());
            test="";
            // if filename not blank set ucf file name
            if(kpath_file!=test){

               // send name to internal variable
               spinwaves::internal::filename=kpath_file;

               // check file exists
               std::ifstream fin(kpath_file);
               if (!fin){
                  terminaltextcolor(RED);
                  std::cerr << "Error - cannot find file \'" << kpath_file << "\' in control statement \'spinwaves:" << word << "\' on line " << line << " of input file" << std::endl;
                  zlog << zTs() << "Error - cannot find file \'" << kpath_file << "\' in control statement \'spinwaves:" << word << "\' on line " << line << " of input file" << std::endl;
                  terminaltextcolor(WHITE);
                  err::vexit();
               }

               return true;
            }
            else{
               terminaltextcolor(RED);
               std::cerr << "Error - empty filename in control statement \'spinwaves:" << word << "\' on line " << line << " of input file" << std::endl;
               terminaltextcolor(WHITE);
               return false;
            }
         }
         test="filetype";
         if(word==test){
            std::string filetype_temp=value;
            filetype_temp.erase(std::remove(filetype_temp.begin(), filetype_temp.end(), '\"'), filetype_temp.end());
            internal::filetype=filetype_temp;

            if (internal::filetype == "path" || internal::filetype== "specific-k"){
               return true;
            }
            else {
               terminaltextcolor(RED);
               std::cerr << "Error - Unexpected filetype for spinwave module. \'spinwaves:" << word << " = " << value << "\' on line " << line << " of input file" << std::endl;
               terminaltextcolor(WHITE);
               return false;
            }
         }
         test="reduction-method";
         if(word==test){
            std::string reduc_string=value;
            reduc_string.erase(std::remove(reduc_string.begin(), reduc_string.end(), '\"'), reduc_string.end());
            internal::reduc_ver=reduc_string;

            if (internal::reduc_ver == "rank0" || internal::reduc_ver == "direct-scatter"){
               return true;
            }
            else {
               terminaltextcolor(RED);
               std::cerr << "Error - Unexpected method for MPI reduction in spinwave model. \'spinwaves:" << word << " = " << value << "\' on line " << line << " of input file" << std::endl;
               terminaltextcolor(WHITE);
               return false;
            }
         }

         test="num-spectrums";
         if(word==test){
            uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
            vin::check_for_valid_int(tt, word, line, prefix, 1, 20,"input","1 - 20");
            internal::nspec = tt;

            // resize spectrum specific arrays
            internal::oss.resize(internal::nspec);
            internal::cm.resize(internal::nspec);
            internal::normk.resize(internal::nspec);
            internal::component.resize(internal::nspec);

            // populate with default values
            for (int j = 0; j < internal::nspec; j++){
               internal::oss[j] = true;
               internal::cm[j] = true;
               internal::normk[j] = true;
               internal::component[j] = "sx";
            }
            return true;
         }
         test="fourier-prefactor";
         if(word==test){

            if (value == "false"){
               internal::prefactor = false;
               return true;
            }
            else if (value == "true"){
               internal::prefactor = true;
               return true;
            }
            else {
               terminaltextcolor(RED);
               std::cerr << "Error - Unknown value in control statement \'spinwaves:" << word << " = " << value << "\' on line " << line << " of input file" << std::endl;
               terminaltextcolor(WHITE);
               return false;
            }
         }
         test="intermediate-structure-factor";
         if(word==test){
            if (value == "false"){
               internal::isf = false;
               return true;
            }
            else if (value == "true"){
               internal::isf = true;
               return true;
            }
            else {
               terminaltextcolor(RED);
               std::cerr << "Error - Unknown value in control statement \'spinwaves:" << word << " = " << value << "\' on line " << line << " of input file" << std::endl;
               terminaltextcolor(WHITE);
               return false;
            }
         }

      }


      // now check if the keys exist for each specific spectrum.
      if (key.find("spinwaves[") != std::string::npos){

         // I've set the maximum number of spectrums = 20. I can't imagine anyone would need that many.
         for (int i = 0 ; i < 20; i++){

            int len=i+1;

            // set the prefix as a substring
            prefix = "spinwaves[" + std::to_string(len) + "]";

            if (key==prefix){

               // array to verify number of spectrums
               internal::super_index_values.push_back(i);

               std::string test="one-sided";
               if(word==test){
                  if (value == "false"){
                     internal::oss[i] = "false";
                     return true;
                  }
                  else if (value == "true"){
                     internal::oss[i] = "false";
                     return true;
                  }
                  else {
                     terminaltextcolor(RED);
                     std::cerr << "Error - Unknown value in control statement \'spinwaves:" << word << " = " << value << "\' on line " << line << " of input file" << std::endl;
                     terminaltextcolor(WHITE);
                     return false;
                  }
               }

               test="material";
               if(word==test){

                  std::vector<int> temp;
                  std::istringstream iss(value);
                  std::string token;

                  // if (value == "all")

                  while (std::getline(iss, token, ',')){
                     //uint64_t tt = std::stoi(token);

                     std::cout << i << " " << std::stoi(token)-1 << std::endl;
                     internal::mat.push_back(std::stoi(token)-1);
                     internal::mat_in_spec.push_back(i);
                  }

                  return true;
               }
               //  ------------------------------------------------------------------


               //  ------------------------------------------------------------------
               // component of magnetisation to use for fourier transform -----------
               //  ------------------------------------------------------------------
               test="complex-magnitude";
               if(word==test){
                  if (value == "false"){
                     internal::cm[i] = false;
                     return true;
                  }
                  else if (value == "true"){
                     internal::cm[i] = true;
                     return true;
                  }
                  else {
                     terminaltextcolor(RED);
                     std::cerr << "Error - Unknown value in control statement \'spinwaves:" << word << " = " << value << "\' on line " << line << " of input file" << std::endl;
                     terminaltextcolor(WHITE);
                     return false;
                  }
               }
               //  ----------------------------------------------------------------

               //  ------------------------------------------------------------------
               // component of magnetisation to use for fourier transform -----------
               //  ------------------------------------------------------------------
               test="normalise-each-k";
               if(word==test){
                  if (value == "false"){
                     internal::normk[i] = false;
                     return true;
                  }
                  else if (value == "true"){
                     internal::normk[i] = true;
                     return true;
                  }
                  else {
                     terminaltextcolor(RED);
                     std::cerr << "Error - Unknown value in control statement \'spinwaves:" << word << " = " << value << "\' on line " << line << " of input file" << std::endl;
                     terminaltextcolor(WHITE);
                     return false;
                  }
               }
               //  ------------------------------------------------------------------

                test="component";
               if(word==test){
                  std::string component_temp=value;
                  component_temp.erase(std::remove(component_temp.begin(), component_temp.end(), '\"'), component_temp.end());
                  internal::component[i]=component_temp;

                  // check input file contains allowed entry.
                  if (component_temp == "sx" || component_temp == "sy" || component_temp == "sz"){
                     return true;
                  }
                  else {
                     terminaltextcolor(RED);
                     std::cerr << "Error - Unexpected method for MPI reduction in spinwave model. \'spinwaves:" << word << " = " << value << "\' on line " << line << " of input file" << std::endl;
                     terminaltextcolor(WHITE);
                     return false;
                  }
               }
            }
         }
      }

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      // add prefix string
      std::string prefix="material:";

      // Check for material id > current array size and if so dynamically expand mp array
      if((unsigned int) super_index + 1 > internal::mp.size() && super_index + 1 < 101) internal::mp.resize(super_index + 1);

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of sw namespace
