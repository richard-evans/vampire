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
#include <cmath>
#include <fstream>
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

   namespace internal{

      //-----------------------------------------------------------------------------
      // Split a comma-separated string into trimmed tokens
      //-----------------------------------------------------------------------------
      std::vector<std::string> split_csv(const std::string& value){

         std::vector<std::string> tokens;
         std::stringstream ss(value);
         std::string token;

         while(std::getline(ss, token, ',')){
            // trim leading and trailing whitespace
            const std::string whitespace = " \t\r\n";
            size_t start = token.find_first_not_of(whitespace);
            if(start == std::string::npos){
               tokens.push_back("");
               continue;
            }
            size_t end = token.find_last_not_of(whitespace);
            tokens.push_back(token.substr(start, end - start + 1));
         }

         return tokens;

      }

      //-----------------------------------------------------------------------------
      // Normalise a 3-vector in place. Returns false if the vector has zero length,
      // in which case it is left unchanged.
      //-----------------------------------------------------------------------------
      bool normalise(double& sx, double& sy, double& sz){

         const double mod_s = std::sqrt(sx*sx + sy*sy + sz*sz);

         if(mod_s < 1.0e-12) return false;

         sx /= mod_s;
         sy /= mod_s;
         sz /= mod_s;

         return true;

      }

      //-----------------------------------------------------------------------------
      // Load a vector field file containing rows of "x,y,z,mx,my,mz" where x,y,z are
      // fractional coordinates of the system size [0:1] and mx,my,mz is the (not
      // necessarily normalised) spin direction at that point. Files are cached so
      // that multiple materials referencing the same file only load it once.
      //-----------------------------------------------------------------------------
      int get_vector_field_id(const std::string& filename, const int line){

         // check cache for an already-loaded file
         for(size_t i = 0; i < vector_field_filenames.size(); i++){
            if(vector_field_filenames[i] == filename) return int(i);
         }

         // open the file
         std::ifstream ifile(filename.c_str());

         if(!ifile.is_open()){
            terminaltextcolor(RED);
            std::cerr << "Error on line " << line << " of material file - vector field file \'" << filename << "\' could not be opened." << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }

         std::vector<field_point_t> points;

         std::string line_str;
         while(std::getline(ifile, line_str)){

            // skip empty lines and comments
            std::vector<std::string> tokens = split_csv(line_str);
            if(tokens.size() == 1 && tokens[0].empty()) continue;
            if(tokens[0].empty()) continue;
            if(tokens[0][0] == '#') continue;

            if(tokens.size() != 6){
               terminaltextcolor(RED);
               std::cerr << "Error - vector field file \'" << filename << "\' contains an invalid line (expected 6 comma-separated values \'x,y,z,mx,my,mz\'): \'" << line_str << "\'" << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }

            field_point_t point;
            point.x  = atof(tokens[0].c_str());
            point.y  = atof(tokens[1].c_str());
            point.z  = atof(tokens[2].c_str());
            point.mx = atof(tokens[3].c_str());
            point.my = atof(tokens[4].c_str());
            point.mz = atof(tokens[5].c_str());

            // normalise the supplied direction vector
            if(!normalise(point.mx, point.my, point.mz)){
               terminaltextcolor(RED);
               std::cerr << "Error - vector field file \'" << filename << "\' contains a point with a zero-length direction vector: \'" << line_str << "\'" << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }

            points.push_back(point);

         }

         if(points.empty()){
            terminaltextcolor(RED);
            std::cerr << "Error on line " << line << " of material file - vector field file \'" << filename << "\' contains no data points." << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }

         // add to cache and return new id
         vector_field_filenames.push_back(filename);
         vector_field_data.push_back(points);

         return int(vector_field_filenames.size() - 1);

      }

   } // end of internal namespace

} // end of spininitialize namespace
