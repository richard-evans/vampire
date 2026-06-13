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
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

// Vampire headers
#include "spininitialize.hpp"
#include "errors.hpp"
#include "vio.hpp"

// spininitialize module headers
#include "internal.hpp"

namespace spininitialize{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for spininitialize module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="spin-initialize";
      if(key!=prefix) return false;

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

      //------------------------------------------------------------------------
      std::string test = "initial-spin-direction";
      if(word == test){

         // split comma separated value into trimmed tokens
         std::vector<std::string> tokens = internal::split_csv(value);

         if(tokens.empty()){
            terminaltextcolor(RED);
            std::cerr << "Error on line " << line << " of material file - " << prefix << "[" << super_index+1 << "]:" << word << " has no value." << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }

         internal::mp_t& mat = internal::mp[super_index];

         const std::string& key = tokens[0];

         //---------------------------------------------------------------
         // random spins (infinite temperature)
         //---------------------------------------------------------------
         if(key == "random"){
            mat.texture = internal::random;
            return true;
         }
         //---------------------------------------------------------------
         // domain-wall-x / domain-wall-y / domain-wall-z, centre, width
         //---------------------------------------------------------------
         else if(key == "domain-wall-x" || key == "domain-wall-y" || key == "domain-wall-z"){
            if(tokens.size() != 3){
               terminaltextcolor(RED);
               std::cerr << "Error on line " << line << " of material file - " << prefix << "[" << super_index+1 << "]:" << word << " = " << key << " requires exactly 2 numerical parameters (centre, width)." << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
            mat.texture = internal::domain_wall;
            mat.axis = (key == "domain-wall-x") ? 0 : (key == "domain-wall-y") ? 1 : 2;
            mat.centre[0] = atof(tokens[1].c_str());
            mat.width     = atof(tokens[2].c_str());
            if(mat.width <= 0.0){
               terminaltextcolor(RED);
               std::cerr << "Error on line " << line << " of material file - " << prefix << "[" << super_index+1 << "]:" << word << " domain wall width must be greater than zero." << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
            return true;
         }
         //---------------------------------------------------------------
         // skyrmion, cx, cy, radius [, chirality, polarity]
         //---------------------------------------------------------------
         else if(key == "skyrmion"){
            if(tokens.size() != 4 && tokens.size() != 6){
               terminaltextcolor(RED);
               std::cerr << "Error on line " << line << " of material file - " << prefix << "[" << super_index+1 << "]:" << word << " = skyrmion requires 3 (cx, cy, radius) or 5 (cx, cy, radius, chirality, polarity) numerical parameters." << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
            mat.texture = internal::skyrmion;
            mat.centre[0] = atof(tokens[1].c_str());
            mat.centre[1] = atof(tokens[2].c_str());
            mat.width     = atof(tokens[3].c_str());
            if(mat.width <= 0.0){
               terminaltextcolor(RED);
               std::cerr << "Error on line " << line << " of material file - " << prefix << "[" << super_index+1 << "]:" << word << " skyrmion radius must be greater than zero." << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
            if(tokens.size() == 6){
               mat.chirality = (atof(tokens[4].c_str()) < 0.0) ? -1 : 1;
               mat.polarity  = (atof(tokens[5].c_str()) < 0.0) ? -1 : 1;
            }
            return true;
         }
         //---------------------------------------------------------------
         // spin-spiral-x / spin-spiral-y / spin-spiral-z, wavelength
         //---------------------------------------------------------------
         else if(key == "spin-spiral-x" || key == "spin-spiral-y" || key == "spin-spiral-z"){
            if(tokens.size() != 2){
               terminaltextcolor(RED);
               std::cerr << "Error on line " << line << " of material file - " << prefix << "[" << super_index+1 << "]:" << word << " = " << key << " requires exactly 1 numerical parameter (wavelength)." << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
            mat.texture = internal::spin_spiral;
            mat.axis = (key == "spin-spiral-x") ? 0 : (key == "spin-spiral-y") ? 1 : 2;
            mat.width = atof(tokens[1].c_str());
            if(mat.width <= 0.0){
               terminaltextcolor(RED);
               std::cerr << "Error on line " << line << " of material file - " << prefix << "[" << super_index+1 << "]:" << word << " spin spiral wavelength must be greater than zero." << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
            return true;
         }
         //---------------------------------------------------------------
         // vortex, cx, cy, core-radius [, chirality, polarity]
         //---------------------------------------------------------------
         else if(key == "vortex"){
            if(tokens.size() != 4 && tokens.size() != 6){
               terminaltextcolor(RED);
               std::cerr << "Error on line " << line << " of material file - " << prefix << "[" << super_index+1 << "]:" << word << " = vortex requires 3 (cx, cy, core-radius) or 5 (cx, cy, core-radius, chirality, polarity) numerical parameters." << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
            mat.texture = internal::vortex;
            mat.centre[0] = atof(tokens[1].c_str());
            mat.centre[1] = atof(tokens[2].c_str());
            mat.width     = atof(tokens[3].c_str());
            if(mat.width < 0.0){
               terminaltextcolor(RED);
               std::cerr << "Error on line " << line << " of material file - " << prefix << "[" << super_index+1 << "]:" << word << " vortex core radius must not be negative." << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
            if(tokens.size() == 6){
               mat.chirality = (atof(tokens[4].c_str()) < 0.0) ? -1 : 1;
               mat.polarity  = (atof(tokens[5].c_str()) < 0.0) ? -1 : 1;
            }
            return true;
         }
         //---------------------------------------------------------------
         // vector-field, <filename>
         //---------------------------------------------------------------
         else if(key == "vector-field"){
            if(tokens.size() != 2){
               terminaltextcolor(RED);
               std::cerr << "Error on line " << line << " of material file - " << prefix << "[" << super_index+1 << "]:" << word << " = vector-field requires exactly 1 filename parameter." << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
            mat.texture = internal::vector_field;
            mat.vector_field_id = internal::get_vector_field_id(tokens[1], line);
            return true;
         }
         //---------------------------------------------------------------
         // otherwise expect a unit vector mx,my,mz (or crystallographic [hkl] notation)
         //---------------------------------------------------------------
         else{

            std::vector<double> u = vin::doubles_from_string(value);

            // check for sane input and normalise if necessary
            vin::check_for_valid_unit_vector(u, word, line, prefix, "material");

            mat.texture = internal::uniform_vector;
            mat.initial_spin[0] = u.at(0);
            mat.initial_spin[1] = u.at(1);
            mat.initial_spin[2] = u.at(2);

            return true;
         }

      }

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of spininitialize namespace
