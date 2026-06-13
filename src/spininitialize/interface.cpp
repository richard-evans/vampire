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
      // material[#]:initial-spin-direction = <texture-name>, <parameters...>
      //
      // This single keyword sets the initial spin configuration for material
      // (super_index+1). The value is a comma-separated list: the first token
      // selects the texture (or, for backwards compatibility, is itself the
      // first component of a plain unit vector), and any remaining tokens are
      // the numerical parameters for that texture (e.g. centre/width for a
      // domain wall). Each "else if" branch below handles one texture type.
      //------------------------------------------------------------------------
      std::string test = "initial-spin-direction";
      if(word == test){

         // split comma separated value into trimmed tokens, e.g.
         // "domain-wall-x, 0.5, 0.2" -> {"domain-wall-x", "0.5", "0.2"}
         std::vector<std::string> tokens = internal::split_csv(value);

         if(tokens.empty()){
            terminaltextcolor(RED);
            std::cerr << "Error on line " << line << " of material file - " << prefix << "[" << super_index+1 << "]:" << word << " has no value." << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }

         // mat is the spin-initialisation record for this material; its
         // fields are filled in by whichever branch below matches
         internal::mp_t& mat = internal::mp[super_index];

         // the first token is either a texture name (e.g. "skyrmion") or,
         // for a plain unit vector, the x-component of the vector itself
         const std::string& key = tokens[0];

         //---------------------------------------------------------------
         // random spins (infinite temperature)
         // syntax: initial-spin-direction = random
         //---------------------------------------------------------------
         if(key == "random"){
            mat.texture = internal::random;
            return true;
         }
         //---------------------------------------------------------------
         // domain-wall-x / domain-wall-y / domain-wall-z, centre, width
         // syntax: initial-spin-direction = domain-wall-x, centre, width
         // centre and width are fractions (0-1) of the system size along
         // the chosen axis; see domain_wall_spin() in textures.cpp
         //---------------------------------------------------------------
         else if(key == "domain-wall-x" || key == "domain-wall-y" || key == "domain-wall-z"){
            // exactly 3 tokens expected: the texture name + 2 numbers
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
         // syntax: initial-spin-direction = skyrmion, cx, cy, radius
         //     or: initial-spin-direction = skyrmion, cx, cy, radius, chirality, polarity
         // (cx,cy) and radius are fractions (0-1) of the system size in the
         // x-y plane; chirality and polarity are each +1 or -1 (default +1);
         // see skyrmion_spin() in textures.cpp
         //---------------------------------------------------------------
         else if(key == "skyrmion"){
            // either 4 tokens (name + cx,cy,radius) or 6 tokens
            // (name + cx,cy,radius,chirality,polarity) are accepted
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
         // syntax: initial-spin-direction = spin-spiral-x, wavelength
         // wavelength is a fraction (0-1) of the system size along the
         // chosen axis; see spin_spiral_spin() in textures.cpp
         //---------------------------------------------------------------
         else if(key == "spin-spiral-x" || key == "spin-spiral-y" || key == "spin-spiral-z"){
            // exactly 2 tokens expected: the texture name + 1 number
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
         // syntax: initial-spin-direction = vortex, cx, cy, core-radius
         //     or: initial-spin-direction = vortex, cx, cy, core-radius, chirality, polarity
         // (cx,cy) and core-radius are fractions (0-1) of the system size in
         // the x-y plane; a core-radius of 0 gives a purely in-plane vortex
         // with no out-of-plane core; chirality and polarity are each +1 or
         // -1 (default +1); see vortex_spin() in textures.cpp
         //---------------------------------------------------------------
         else if(key == "vortex"){
            // either 4 tokens (name + cx,cy,core-radius) or 6 tokens
            // (name + cx,cy,core-radius,chirality,polarity) are accepted
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
         // syntax: initial-spin-direction = vector-field, filename
         // filename contains rows "x,y,z,mx,my,mz" (fractional coordinates
         // and a direction vector) which are interpolated for each atom by
         // inverse-distance weighting; see get_vector_field_id() in
         // parse.cpp and vector_field_spin() in vector_field.cpp
         //---------------------------------------------------------------
         else if(key == "vector-field"){
            // exactly 2 tokens expected: the texture name + the filename
            if(tokens.size() != 2){
               terminaltextcolor(RED);
               std::cerr << "Error on line " << line << " of material file - " << prefix << "[" << super_index+1 << "]:" << word << " = vector-field requires exactly 1 filename parameter." << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
            mat.texture = internal::vector_field;
            // load the file (or retrieve it from the cache if another
            // material already references the same filename) and store its
            // index for later use by vector_field_spin()
            mat.vector_field_id = internal::get_vector_field_id(tokens[1], line);
            return true;
         }
         //---------------------------------------------------------------
         // otherwise expect a unit vector mx,my,mz (or crystallographic [hkl] notation)
         // syntax: initial-spin-direction = mx,my,mz
         //     or: initial-spin-direction = [hkl]
         // every atom of this material is then initialised pointing along
         // this single fixed direction; see uniform_vector_spin() in
         // textures.cpp
         //---------------------------------------------------------------
         else{

            // parse the full (un-split) value as a list of doubles, so that
            // crystallographic notation such as "[1-10]" is handled
            // correctly by doubles_from_string
            std::vector<double> u = vin::doubles_from_string(value);

            // check for sane input (3 components) and normalise to a unit
            // vector if necessary
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
