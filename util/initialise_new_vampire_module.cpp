//------------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//
//    Simple C++ code to initialise placeholder files for new vampire module
//
//    This code creates the following files with the following purposes:
//
//          data.cpp          // Store all data structures local to this module
//          interface.cpp     // Functions for defining the user interface
//          initialise.cpp    // Function to initialise the module
//          internal.hpp      // Header file listing module shared variables
//          <module>.hpp      // Header file for extenrally visible functions
//          makefile          // Module makefile
//
//    This program is invoked using:
//
//          invm --namespace <namespace_name>
//
//    The following options are also accepted
//          --author <author name>
//          --email <email address>
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "stdlib.h"
#include <algorithm>

// Forward declaration of functions
void process_command_line(int argc, char* argv[], std::string& namespace_name, std::string& author, std::string& email, std::string& year);
std::string create_file_header(std::string author, std::string email, std::string year);
bool file_exists(const std::string& file_name);
void create_data(const std::string& file_header, const std::string& namespace_name);
void create_interface(const std::string& file_header, const std::string& namespace_name);
void create_initialise(const std::string& file_header, const std::string& namespace_name);
void create_internal(const std::string& file_header, const std::string& nn);
void create_module(const std::string& file_header, const std::string& nn);
void create_makefile(const std::string& file_header, const std::string& nn);

// global constants
const std::string blank = "";

//------------------------------------------------------------------------------
// Main function to generate module files in current directory
//------------------------------------------------------------------------------
int main(int argc, char* argv[]){

   // Strings for command line options
   std::string namespace_name="";
   std::string author="";
   std::string email="";
   std::string year="2018";

   // determine namespace name, author and email from command line
   process_command_line(argc, argv, namespace_name, author, email, year);

   // create file header
   std::string file_header = create_file_header(author, email, year);

   // Generate data.cpp
   create_data(file_header, namespace_name);

   // Generate interface.cpp
   create_interface(file_header, namespace_name);

   // Generate initialise.cpp
   create_initialise(file_header, namespace_name);

   // Generate internal.hpp
   create_internal(file_header, namespace_name);

   // Generate internal.hpp
   create_module(file_header, namespace_name);

   // Generate makefile
   create_makefile(file_header, namespace_name);

   return EXIT_SUCCESS;

}

//------------------------------------------------------------------------------
// Function to process command line arguments
//------------------------------------------------------------------------------
void process_command_line(int argc, char* argv[],
                          std::string& namespace_name,
                          std::string& author,
                          std::string& email,
                          std::string& year){

   std::cout << "Processing command line arguments" << std::endl;

   for(int arg = 1; arg < argc; arg++){
      std::string sw=argv[arg];
      //---------------------------------------------
      // namespace name
      //---------------------------------------------
      if(sw=="--namespace"){
         // check number of args not exceeded
         if(arg+1 < argc){
            arg++;
            namespace_name = std::string(argv[arg]);
         }
         else{
            std::cerr << "Error - no namespace name specified for \'--namespace\' command line option" << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      //---------------------------------------------
      // author
      //---------------------------------------------
      else if(sw=="--author"){
         // check number of args not exceeded
         if(arg+1 < argc){
            arg++;
            author = std::string(argv[arg]);
         }
         else{
            std::cerr << "Error - no author specified for \'--author\' command line option" << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      //---------------------------------------------
      // email
      //---------------------------------------------
      else if(sw=="--email"){
         // check number of args not exceeded
         if(arg+1 < argc){
            arg++;
            email = std::string(argv[arg]);
         }
         else{
            std::cerr << "Error - no email address specified for \'--email\' command line option" << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      //---------------------------------------------
      // year
      //---------------------------------------------
      else if(sw=="--year"){
         // check number of args not exceeded
         if(arg+1 < argc){
            arg++;
            year = std::string(argv[arg]);
         }
         else{
            std::cerr << "Error - no year specified for \'--year\' command line option" << std::endl;
            exit(EXIT_FAILURE);
         }
      }
      else{
         std::cerr << "Error - unknown command line parameter \'" << sw << "\'" << std::endl;
         exit(EXIT_FAILURE);
      }
   }

   // Check for valid initialisation
   if(namespace_name == blank){
      std::cerr << "Error - no namespace name specified - use the \'--namespace <name> \' command line option" << std::endl;
      exit(EXIT_FAILURE);
   }

   // Remove .,+-={}[]() from name and convert tabs and spaces to underscores
   namespace_name.erase(remove(namespace_name.begin(), namespace_name.end(),'.'), namespace_name.end());
   namespace_name.erase(remove(namespace_name.begin(), namespace_name.end(),','), namespace_name.end());
   namespace_name.erase(remove(namespace_name.begin(), namespace_name.end(),'+'), namespace_name.end());
   namespace_name.erase(remove(namespace_name.begin(), namespace_name.end(),'-'), namespace_name.end());
   namespace_name.erase(remove(namespace_name.begin(), namespace_name.end(),'='), namespace_name.end());
   namespace_name.erase(remove(namespace_name.begin(), namespace_name.end(),'['), namespace_name.end());
   namespace_name.erase(remove(namespace_name.begin(), namespace_name.end(),']'), namespace_name.end());
   namespace_name.erase(remove(namespace_name.begin(), namespace_name.end(),'{'), namespace_name.end());
   namespace_name.erase(remove(namespace_name.begin(), namespace_name.end(),'}'), namespace_name.end());
   namespace_name.erase(remove(namespace_name.begin(), namespace_name.end(),'('), namespace_name.end());
   namespace_name.erase(remove(namespace_name.begin(), namespace_name.end(),')'), namespace_name.end());
   std::replace (namespace_name.begin(), namespace_name.end(),' ', '_');
   std::replace (namespace_name.begin(), namespace_name.end(),'\t', '_');

   std::cout << "The following details have been entered:" << std::endl;
   std::cout << "   Namespace name:  " << namespace_name << std::endl;
   std::cout << "   Author:          " << author << std::endl;
   std::cout << "   Email:           " << email << std::endl;
   std::cout << "   Year:            " << year << std::endl;
   std::cout << "Are these correct (Y/N)? ";
   std::string check;
   check = std::cin.get();
   std::string Yes = "Y";
   std::string yes = "y";
   if(check != Yes && check != yes){
      std::cerr << "Aborting module initialisation" << std::endl;
      exit(EXIT_FAILURE);
   }

   return;

}

//---------------------------------------------------------------------------
// Function to create file header
//---------------------------------------------------------------------------
std::string create_file_header(std::string author, std::string email, std::string year){

   // dec;are temporary string stream
   std::stringstream cfh_ss;
   cfh_ss << "//------------------------------------------------------------------------------\n";
   cfh_ss << "//\n";
   cfh_ss << "//   This file is part of the VAMPIRE open source package under the\n";
   cfh_ss << "//   Free BSD licence (see licence file for details).\n";
   cfh_ss << "//\n";
   if(author != blank){
      cfh_ss << "//   (c) " << author << " " << year << ". All rights reserved.\n";
      cfh_ss << "//\n";
   }
   if(email  != blank){
      cfh_ss << "//   Email: " << email << "\n";
      cfh_ss << "//\n";
   }
   cfh_ss << "//------------------------------------------------------------------------------\n";
   cfh_ss << "//\n";

   return cfh_ss.str();

}

//------------------------------------------------------------------------------
// Function to determine if a file already exists and is accessible
//------------------------------------------------------------------------------
bool file_exists(const std::string& file_name){
    std::ifstream infile(file_name.c_str());
    return infile.good();
}

//------------------------------------------------------------------------------
// Function to capitalize string
//------------------------------------------------------------------------------
std::string capitalize(std::string str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
}

//------------------------------------------------------------------------------
// Function to create data.ccp file
//------------------------------------------------------------------------------
void create_data(const std::string& file_header, const std::string& nn){

   // nn = namespace name

   // Check to see if file exists
   if(file_exists("data.cpp")){
      std::cout << "data.cpp already exists - skipping initialisation." << std::endl;
      return;
   }

   std::cout << "Generating data.cpp" << std::endl;

   // Open file
   std::ofstream ofile;
   ofile.open("data.cpp");

   // Write header
   ofile << file_header << std::endl;
   ofile << "// C++ standard library headers\n" << std::endl;
   ofile << "// Vampire headers" << std::endl;
   // Include module header
   ofile << "#include \"" << nn << ".hpp\"\n" << std::endl;
   ofile << "// " << nn << " module headers" << std::endl;
   ofile << "#include \"internal.hpp\"\n" << std::endl;

   // Namespaces
   ofile << "namespace " << nn << "{\n" << std::endl;
   ofile << "   //------------------------------------------------------------------------------" << std::endl;
   ofile << "   // Externally visible variables" << std::endl;
   ofile << "   //------------------------------------------------------------------------------\n" << std::endl;
   ofile << "   namespace internal{\n" << std::endl;
   ofile << "      //------------------------------------------------------------------------" << std::endl;
   ofile << "      // Shared variables inside " << nn << " module" << std::endl;
   ofile << "      //------------------------------------------------------------------------\n" << std::endl;
   ofile << "      bool enabled; // bool to enable module\n" << std::endl;
   ofile << "      std::vector<internal::mp_t> mp; // array of material properties\n" << std::endl;
   ofile << "   } // end of internal namespace\n" << std::endl;
   ofile << "} // end of " << nn << " namespace\n" << std::endl;

   ofile.close();

   return;

}

void create_interface(const std::string& file_header, const std::string& nn){

   // Check to see if file exists
   if(file_exists("interface.cpp")){
      std::cout << "interface.cpp already exists - skipping initialisation." << std::endl;
      return;
   }

   std::cout << "Generating interface.cpp" << std::endl;

   // Open file
   std::ofstream ofile;
   ofile.open("interface.cpp");

   // Write header
   ofile << file_header << std::endl;
   ofile << "// C++ standard library headers" << std::endl;
   ofile << "#include <string>\n" << std::endl;

   ofile << "// Vampire headers" << std::endl;
   ofile << "#include \"" << nn << ".hpp\"" << std::endl;
   ofile << "#include \"errors.hpp\"" << std::endl;
   ofile << "#include \"vio.hpp\"\n" << std::endl;
   ofile << "// " << nn << " module headers" << std::endl;
   ofile << "#include \"internal.hpp\"\n" << std::endl;

   // Namespaces
   ofile << "namespace " << nn << "{\n" << std::endl;
   ofile << "   //---------------------------------------------------------------------------" << std::endl;
   ofile << "   // Function to process input file parameters for " << nn << " module" << std::endl;
   ofile << "   //---------------------------------------------------------------------------" << std::endl;
   ofile << "   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){\n" << std::endl;
   ofile << "      // Check for valid key, if no match return false" << std::endl;
   ofile << "      std::string prefix=\"" << nn << "\";" << std::endl;
   ofile << "      if(key!=prefix) return false;\n" << std::endl;
   ofile << "      //--------------------------------------------------------------------" << std::endl;
   ofile << "      // Keyword not found" << std::endl;
   ofile << "      //--------------------------------------------------------------------" << std::endl;
   ofile << "      return false;\n" << std::endl;
   ofile << "   }\n" << std::endl;
   ofile << "   //---------------------------------------------------------------------------" << std::endl;
   ofile << "   // Function to process material parameters" << std::endl;
   ofile << "   //---------------------------------------------------------------------------" << std::endl;
   ofile << "   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){\n" << std::endl;
   ofile << "      // add prefix string" << std::endl;
   ofile << "      std::string prefix=\"material:\";\n" << std::endl;
   ofile << "      // Check for material id > current array size and if so dynamically expand mp array" << std::endl;
   ofile << "      if((unsigned int) super_index + 1 > internal::mp.size() && super_index + 1 < 101) internal::mp.resize(super_index + 1);\n" << std::endl;
   ofile << "      //--------------------------------------------------------------------" << std::endl;
   ofile << "      // Keyword not found" << std::endl;
   ofile << "      //--------------------------------------------------------------------" << std::endl;
   ofile << "      return false;\n" << std::endl;
   ofile << "   }\n" << std::endl;
   ofile << "} // end of " << nn << " namespace\n" << std::endl;

   ofile.close();

   return;

}

void create_initialise(const std::string& file_header, const std::string& nn){

   // Check to see if file exists
   if(file_exists("initialize.cpp")){
      std::cout << "initialize.cpp already exists - skipping initialisation." << std::endl;
      return;
   }

   std::cout << "Generating initialize.cpp" << std::endl;

   // Open file
   std::ofstream ofile;
   ofile.open("initialize.cpp");

   // Write header
   ofile << file_header << std::endl;
   ofile << "// C++ standard library headers\n" << std::endl;
   ofile << "// Vampire headers" << std::endl;

   // Include module header
   ofile << "#include \"" << nn << ".hpp\"\n" << std::endl;
   ofile << "// " << nn << " module headers" << std::endl;
   ofile << "#include \"internal.hpp\"\n" << std::endl;

   // Namespaces
   ofile << "namespace " << nn << "{\n" << std::endl;
   ofile << "   //----------------------------------------------------------------------------" << std::endl;
   ofile << "   // Function to initialize " << nn << " module" << std::endl;
   ofile << "   //----------------------------------------------------------------------------" << std::endl;
   ofile << "   void initialize(){\n" << std::endl;
   ofile << "      return;\n" << std::endl;
   ofile << "   }\n" << std::endl;
   ofile << "} // end of " << nn << " namespace\n" << std::endl;

   ofile.close();

   return;

}

void create_internal(const std::string& file_header, const std::string& nn){

   // Check to see if file exists
   if(file_exists("internal.hpp")){
      std::cout << "internal.hpp already exists - skipping initialisation." << std::endl;
      return;
   }

   std::cout << "Generating internal.hpp" << std::endl;

   // create header guard strings
   std::string CAP_NN = capitalize(nn) + "_INTERNAL_H_";
   std::stringstream hguard_top_ss;
   hguard_top_ss << "#ifndef " << CAP_NN << "\n" << "#define " << CAP_NN;
   std::stringstream hguard_bottom_ss;
   hguard_bottom_ss << "#endif //" << CAP_NN;

   // Open file
   std::ofstream ofile;
   ofile.open("internal.hpp");

   // Write header
   ofile << file_header << std::endl;
   // Write header guard
   ofile << hguard_top_ss.str() << std::endl;
   ofile << "//" << std::endl;
   ofile << "//---------------------------------------------------------------------" << std::endl;
   ofile << "// This header file defines shared internal data structures and" << std::endl;
   ofile << "// functions for the " << nn << " module. These functions and" << std::endl;
   ofile << "// variables should not be accessed outside of this module." << std::endl;
   ofile << "//---------------------------------------------------------------------\n" << std::endl;

   ofile << "// C++ standard library headers\n" << std::endl;
   ofile << "// Vampire headers" << std::endl;
   ofile << "#include \"" << nn << ".hpp\"\n" << std::endl;
   ofile << "// " << nn << " module headers" << std::endl;
   ofile << "#include \"internal.hpp\"\n" << std::endl;

   // Namespaces
   ofile << "namespace " << nn << "{\n" << std::endl;
   ofile << "   namespace internal{\n" << std::endl;
   ofile << "      //-------------------------------------------------------------------------" << std::endl;
   ofile << "      // Internal data type definitions" << std::endl;
   ofile << "      //-------------------------------------------------------------------------\n" << std::endl;
   ofile << "      //-----------------------------------------------------------------------------" << std::endl;
   ofile << "      // internal materials class for storing material parameters" << std::endl;
   ofile << "      //-----------------------------------------------------------------------------" << std::endl;
   ofile << "      class mp_t{\n" << std::endl;
   ofile << "          private:\n" << std::endl;
   ofile << "          public:\n" << std::endl;
   ofile << "             //------------------------------" << std::endl;
   ofile << "             // material parameter variables" << std::endl;
   ofile << "             //------------------------------" << std::endl;
   ofile << "             double test;\n" << std::endl;
   ofile << "             // constructor" << std::endl;
   ofile << "             mp_t (const unsigned int max_materials = 100):" << std::endl;
   ofile << "                test(0.0) // constructor initialisation of test variable" << std::endl;
   ofile << "             {" << std::endl;
   ofile << "                // constructor body for initialising more complex data/arrays" << std::endl;
   ofile << "             }; // end of constructor\n" << std::endl;
   ofile << "       }; // end of internal::mp class\n" << std::endl;
   ofile << "      //-------------------------------------------------------------------------" << std::endl;
   ofile << "      // Internal shared variables" << std::endl;
   ofile << "      //-------------------------------------------------------------------------\n" << std::endl;
   ofile << "      extern bool enabled; // bool to enable module\n" << std::endl;
   ofile << "      extern std::vector<internal::mp_t> mp; // array of material properties\n" << std::endl;
   ofile << "      //-------------------------------------------------------------------------" << std::endl;
   ofile << "      // Internal function declarations" << std::endl;
   ofile << "      //-------------------------------------------------------------------------\n" << std::endl;
   ofile << "   } // end of internal namespace\n" << std::endl;
   ofile << "} // end of " << nn << " namespace\n" << std::endl;

   // Write header guard
   ofile << hguard_bottom_ss.str() << std::endl;

   ofile.close();

   return;

}

void create_module(const std::string& file_header, const std::string& nn){

   std::stringstream mfn_ss;
   mfn_ss << nn << ".hpp";
   std::string mfn = mfn_ss.str();

   // Check to see if file exists
   if(file_exists(mfn)){
      std::cout << mfn << " already exists - skipping initialisation." << std::endl;
      return;
   }

   std::cout << "Generating " << mfn  << std::endl;

   // create header guard strings
   std::string CAP_NN = capitalize(nn) + "_H_";
   std::stringstream hguard_top_ss;
   hguard_top_ss << "#ifndef " << CAP_NN << "\n" << "#define " << CAP_NN;
   std::stringstream hguard_bottom_ss;
   hguard_bottom_ss << "#endif //" << CAP_NN;

   // Open file
   std::ofstream ofile;
   ofile.open(mfn);

   // Write header
   ofile << file_header << std::endl;
   // Write header guard
   ofile << hguard_top_ss.str() << "\n" << std::endl;

   ofile << "// C++ standard library headers" << std::endl;
   ofile << "#include <string>\n" << std::endl;
   ofile << "// Vampire headers" << std::endl;
   ofile << "#include \"" << nn << ".hpp\"\n" << std::endl;

   ofile << "//--------------------------------------------------------------------------------" << std::endl;
   ofile << "// Namespace for variables and functions for " << nn << " module" << std::endl;
   ofile << "//--------------------------------------------------------------------------------" << std::endl;
   ofile << "namespace " << nn << "{\n" << std::endl;

   ofile << "   //-----------------------------------------------------------------------------" << std::endl;
   ofile << "   // Function to initialise " << nn << " module" << std::endl;
   ofile << "   //-----------------------------------------------------------------------------" << std::endl;
   ofile << "   void initialize();\n" << std::endl;
   ofile << "   //---------------------------------------------------------------------------" << std::endl;
   ofile << "   // Function to process input file parameters for " << nn << " module" << std::endl;
   ofile << "   //---------------------------------------------------------------------------" << std::endl;
   ofile << "   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);\n" << std::endl;
   ofile << "   //---------------------------------------------------------------------------" << std::endl;
   ofile << "   // Function to process material parameters" << std::endl;
   ofile << "   //---------------------------------------------------------------------------" << std::endl;
   ofile << "   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);\n" << std::endl;

   ofile << "} // end of " << nn << " namespace\n" << std::endl;

   // Write header guard
   ofile << hguard_bottom_ss.str() << std::endl;

   ofile.close();

   return;

}

void create_makefile(const std::string& file_header, const std::string& nn){

   // Check to see if file exists
   if(file_exists("makefile")){
      std::cout << "makefile already exists - skipping initialisation." << std::endl;
      return;
   }

   std::cout << "Generating makefile" << std::endl;

   // Open file
   std::ofstream ofile;
   ofile.open("makefile");

   // Write header
   ofile << "#--------------------------------------------------------------" << std::endl;
   ofile << "#          Makefile for " << nn << " module" << std::endl;
   ofile << "#--------------------------------------------------------------\n" << std::endl;

   ofile << "# List module object filenames" << std::endl;
   ofile << nn << "_objects =\\" << std::endl;
   ofile << "data.o \\" << std::endl;
   ofile << "initialize.o \\" << std::endl;
   ofile << "interface.o\n" << std::endl;

   ofile << "# Append module objects to global tree" << std::endl;
   ofile << "OBJECTS+=$(addprefix obj/" << nn << "/,$(" << nn << "_objects))" << std::endl;

   ofile.close();

      return;

}
