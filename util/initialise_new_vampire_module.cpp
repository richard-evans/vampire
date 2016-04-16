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

// Forward declaration of functions
void process_command_line(int argc, char* argv[], std::string& namespace_name, std::string& author, std::string& email);
std::string create_file_header(std::string author, std::string email);

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

   // determine namespace name, author and email from command line
   process_command_line(argc, argv, namespace_name, author, email);

   // create file headers
   std::string file_header = create_file_header(author, email);

   return EXIT_SUCCESS;

}

//------------------------------------------------------------------------------
// Function to process command line arguments
//------------------------------------------------------------------------------
void process_command_line(int argc, char* argv[],
                          std::string& namespace_name,
                          std::string& author,
                          std::string& email){

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
std::string create_file_header(std::string author, std::string email){

   // dec;are temporary string stream
   std::stringstream cfh_ss;
   cfh_ss << "//------------------------------------------------------------------------------\n";
   cfh_ss << "//\n";
   cfh_ss << "//   This source file is part of the VAMPIRE open source package under the\n";
   cfh_ss << "//   GNU GPL (version 2) licence (see licence file for details).\n";
   cfh_ss << "//\n";
   if(author != blank){
      cfh_ss << "//   (c) " << author << " 2016. All rights reserved.\n";
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
