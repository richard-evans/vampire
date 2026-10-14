//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Daniel Meilak 2021. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include<string>
#include<vector>
#include<fstream>
#include<locale>
#include<cctype>
#include<algorithm>

// vdc headers
#include"vdc.hpp"

namespace vdc {

// forward function declarations
void read(std::vector<input_t> &input_list);
void set(const std::vector<input_t> &input_list);
void checks();
void init_vector_y();


//------------------------------------------------------------------------------
// Input file parsing function
//------------------------------------------------------------------------------
void read_and_set(){

   // vector of vectors to store input key and values
   std::vector<input_t> input_list;

   // read from input file
   read(input_list);

   // set all required parameters
   set(input_list);

   // check input parameter dependencies
   checks();
}

//------------------------------------------------------------------------------
// open input file and read all input parameters and values
//------------------------------------------------------------------------------
void read(std::vector<input_t> &input_list){

   // check if input file exists
   std::ifstream input_file(vdc::input_file);
   if (!input_file.is_open()){
      //std::cerr << "No input file provided. Running default configuration." << std::endl;
      return;
   }
   else { std::cout << "Reading from input file '" << vdc::input_file << "'\n"; }

   // read input into line
   std::string line;

   unsigned int line_number = 1;
   std::string delimiters = " ,(){}[]:=!";

   // temp variables
   std::string temp_val;

   // work through input_file
   while(std::getline(input_file,line)){

      // remove everything after comment character
      line = line.substr(0,line.find('#'));

      // skip empty lines
      if (line.empty()){
         line_number++;
         continue;
      }

      input_t input;
      bool key_set = false;

      // work through line
      for (char &c : line){

         // if uppercase, make lowercase
         c = std::tolower(c);

         // skip delimiters and push last value
         if (delimiters.find(c) != std::string::npos){

            // in case delimiters follow each other
            if (temp_val == ""){ continue; }

            // set key or add to values
            if (key_set){ input.value.push_back(temp_val); }
            else {
               key_set = true;
               input.key = temp_val;
            }

            temp_val.clear();
         }
         // otherwise add char to temp_val
         else { temp_val.push_back(c); }
      }

      // push back last value if not empty
      if (temp_val != ""){
         if (key_set){ input.value.push_back(temp_val); }
         else {
                  key_set = true;
                  input.key = temp_val;
         }
      }

      // add to input_list and empty temp storage
      input.line_number = line_number;
      input_list.push_back(input);
      temp_val.clear();
      line_number++;
   }

   input_file.close();
}

//------------------------------------------------------------------------------
// set all input parameters in input
//------------------------------------------------------------------------------
void set(const std::vector<input_t> &input_list){

   // loop through input lines
   for (const input_t &input : input_list){

      // if key has already been set at command line, print warning and skip
      if (std::find(vdc::cmdl_parameters.begin(), vdc::cmdl_parameters.end(), input.key) != vdc::cmdl_parameters.end() ){
         std::cerr << "Warning - the parameter '" << input.key << "' has been set in the command line\n"
                   << "The same parameter on line " << input.line_number << " of input file '" << vdc::input_file
                   << "' is being ignored\n";
         continue;
      }

      // check key, print error if not found
      if (vdc::key_list.find(input.key) == vdc::key_list.end()){
         std::cerr << "Error - Unknown control statement '" << input.key << "' on line "
                   << input.line_number << " of input file '" << vdc::input_file << "'" << std::endl;
         std::exit(EXIT_FAILURE);
      }

      // set parameters using function wrapper and arguments
      vdc::key_list.at(input.key)(input);
   }
}

//---------------------------------------------------------------------------
// Additional checks on input file parameters
//---------------------------------------------------------------------------
void checks(){

   // check for valid axis initialisations
   if ( z_vector && !x_vector ){

      // check for a z-axis with vector_z[2] = 0
      // (Hence z_vector lies in xy-plane)
      if ( (-0.000001 < vdc::vector_z[2]) && (vdc::vector_z[2] < 0.000001) ){
         // x-axis will lie along {0,0,1}
         vdc::vector_x = {0.0, 0.0, 1.0};

         // find vector_y
         init_vector_y();
      }
      else {
         // find x-axis which lies on plane with normal vector_z
         // (there must exist a vector with coor {1,0,x} normal to z_vector )
         vdc::vector_x = {1.0, 0.0, -1.0*vdc::vector_z[0]/vdc::vector_z[2]};

         // find vector_y
         init_vector_y();
      }
   }
   else if ( !z_vector && x_vector ){

      // x-axis cannot be initialised alone
      std::cerr << "Error - x-axis cannot be initialised alone."
                << "\n" << "To use 1D colour scheme, initialise z-axis instead"
                << std::endl;
      std::exit(EXIT_FAILURE);
   }
   else if ( z_vector && x_vector ){

      // check if input axes are orthogonal
      double zdotx;
      zdotx = vdc::vector_z[0]*vdc::vector_x[0] + vdc::vector_z[1]*vdc::vector_x[1] + vdc::vector_z[2]*vdc::vector_x[2];

      if ( (zdotx > 0.000001) || (zdotx < -0.000001) ){
         std::cerr << "Error - input axes are not orthogonal." << std::endl;
         std::exit(EXIT_FAILURE);
      }

   }
}

//------------------------------------------------------------------------------
// Perform cross product of input vectors vector_x and vector_z to get vector_y
//------------------------------------------------------------------------------
void init_vector_y(){

   vdc::vector_y[0] = vdc::vector_z[1]*vdc::vector_x[2] - vdc::vector_x[1]*vdc::vector_z[2];
   vdc::vector_y[1] = vdc::vector_x[0]*vdc::vector_z[2] - vdc::vector_z[0]*vdc::vector_x[2];
   vdc::vector_y[2] = vdc::vector_z[0]*vdc::vector_x[1] - vdc::vector_x[0]*vdc::vector_z[1];

   return;
}

} // end of namespace vdc
