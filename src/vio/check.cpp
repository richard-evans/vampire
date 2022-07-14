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

// Vampire headers
// Headers
#include "vio.hpp"
#include "errors.hpp"
#include "units.hpp"

// vio module headers
#include "internal.hpp"

namespace vin{

    ///-----------------------------------------------------------------------
    /// Function to check for correct unit type and valid variable range
    ///-----------------------------------------------------------------------
    ///
    void check_for_valid_value(double& value, /// value of variable as in input file
                                        std::string word, /// input file keyword
                                        int line, /// input file line
                                        std::string prefix, /// input file prefix
                                        std::string unit, /// unit specified in input file
                                        std::string unit_type, /// expected unit type
                                        double range_min, /// acceptable minimum value for variable
                                        double range_max, /// acceptable maximum value for variable
                                        std::string input_file_type, ///input file name
                                        std::string range_text) /// customised text
    {

        // Define test unit
        std::string test_unit_type=unit_type;

        // Define integer for unit conversion status
        int convert_status=0;

        // If no unit given, assume internal, otherwise convert to internal units
        if(unit.size() != 0) convert_status = units::convert(unit,value,test_unit_type);

        // Test for valid conversion
        if(convert_status==EXIT_FAILURE){
            terminaltextcolor(RED);
            std::cerr << "Error: Unit \'" << unit << "\' specified on line " << line << " of " << input_file_type << " file is not a valid unit." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error: Unit \'" << unit << "\' specified on line " << line << " of " << input_file_type << " file is not a valid unit." << std::endl;
            err::vexit();
        }

        // Test for change in unit type in case of wrong unit type
        if(unit_type!=test_unit_type){
            terminaltextcolor(RED);
            std::cerr << "Error: Unit \'" << unit << "\' of type \'" << test_unit_type << "\' specified on line " << line << " of " << input_file_type << " is invalid for parameter " << prefix << word << "."<< std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error: Unit \'" << unit << "\' of type \'" << test_unit_type << "\' specified on line " << line << " of " << input_file_type << " is invalid for parameter " << prefix << word << "."<< std::endl;
            err::vexit();
        }

        // Check for valid range
        if((fabs(value)<range_min) || (fabs(value)>range_max)){
            terminaltextcolor(RED);
            std::cerr << "Error: " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error: " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
            err::vexit();
        }

        // Success - input is sane!
        return;

    }

   ///-----------------------------------------------------------------------
   /// Function to check for correct unit type and valid variable range
   ///-----------------------------------------------------------------------
   ///
   void check_for_valid_positive_value(double& value, /// value of variable as in input file
                                     std::string word, /// input file keyword
                                     int line, /// input file line
                                     std::string prefix, /// input file prefix
                                     std::string unit, /// unit specified in input file
                                     std::string unit_type, /// expected unit type
                                     double range_min, /// acceptable minimum value for variable
                                     double range_max, /// acceptable maximum value for variable
                                     std::string input_file_type, ///input file name
                                     std::string range_text) /// customised text
   {

      // check for correct conversion and absolute value
      check_for_valid_value(value, word, line, prefix, unit, unit_type, range_min, range_max, input_file_type, range_text);

      // check for positive value
      if(value < 0.0){
          terminaltextcolor(RED);
          std::cerr << "Error: " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be a positive constant " << range_text << "." << std::endl;
          terminaltextcolor(WHITE);
          zlog << zTs() << "Error: " << prefix << ":"<< word << " on line " << line << " of " << input_file_type << " file must be a positive constant " << range_text << "." << std::endl;
          err::vexit();
      }

      // Success - input is sane!
      return;

   }
    ///
    /// Function to check for valid int variable range
    ///-----------------------------------------------------------------------
    ///
    void check_for_valid_int(  int& value, /// value of variable as in input file
                                        std::string word, /// input file keyword
                                        int line, /// input file line
                                        std::string prefix, /// input file prefix
                                        int range_min, /// acceptable minimum value for variable
                                        int range_max, /// acceptable maximum value for variable
                                        std::string input_file_type, ///input file name
                                        std::string range_text) /// customised text
    {

        // Check for valid range
        if((value<range_min) || (value>range_max)){
        terminaltextcolor(RED);
            std::cerr << "Error: " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
            terminaltextcolor(WHITE);
        zlog << zTs() << "Error: " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
            err::vexit();
        }

        // Success - input is sane!
        return;

    }

    ///
    /// Overloaded function to check for valid uint variable range
    ///-----------------------------------------------------------------------
    ///
    void check_for_valid_int(  unsigned int& value, /// value of variable as in input file
                                        std::string word, /// input file keyword
                                        int line, /// input file line
                                        std::string prefix, /// input file prefix
                                        unsigned int range_min, /// acceptable minimum value for variable
                                        unsigned int range_max, /// acceptable maximum value for variable
                                        std::string input_file_type, ///input file name
                                        std::string range_text) /// customised text
    {

        // Check for valid range
        if((value<range_min) || (value>range_max)){
        terminaltextcolor(RED);
            std::cerr << "Error: " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
            terminaltextcolor(WHITE);
        zlog << zTs() << "Error: " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
            err::vexit();
        }

        // Success - input is sane!
        return;

    }

    //-----------------------------------------------------------------------
    // Overloaded function to check for valid uint64_t variable range
    //-----------------------------------------------------------------------
    void check_for_valid_int(uint64_t& value,             // value of variable as in input file
                             std::string word,            // input file keyword
                             int line,                    // input file line
                             std::string prefix,          // input file prefix
                             uint64_t range_min,          // acceptable minimum value for variable
                             uint64_t range_max,          // acceptable maximum value for variable
                             std::string input_file_type, // input file name
                             std::string range_text)      // customised text
    {

        // Check for valid range
        if( (value < range_min) || (value > range_max) ){
           terminaltextcolor(RED);
           std::cerr << "Error: " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
           terminaltextcolor(WHITE);
           zlog << zTs() << "Error: " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
           err::vexit();
        }

        // Success - input is sane!
        return;

    }

    ///-----------------------------------------------------------------------
    /// Function to check for valid boolean
    ///
    /// (c) R F L Evans 2013
    ///
    /// If input is invalid, then function will output error message and
    /// program will exit from here. Otherwise returns a sanitised bool.
    ///
    ///-----------------------------------------------------------------------
    bool check_for_valid_bool( std::string value, /// variable as in input file
                                        std::string word, /// input file keyword
                                        int line, /// input file line
                                        std::string prefix, /// input file prefix
                                        std::string input_file_type) ///input file name
    {
        // Define string constants
        const std::string t="true";
        const std::string f="false";
        const std::string b="";

        // Check for three possible correct answers
        if(value==t) return true;
        if(value==f) return false;
        if(value==b) return true;

        // Invalid input - print error and exit
        terminaltextcolor(RED);
        std::cerr << "Error: " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be true or false." << std::endl;
        terminaltextcolor(WHITE);
        zlog << zTs() << "Error: " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be true or false." << std::endl;
        err::vexit();

        return false;
    }

    ///
    /// Function to check for correct 3-component vector and ensure length of 1
    ///-------------------------------------------------------------------------
    ///
    void check_for_valid_unit_vector(std::vector<double>& u, /// unit vector
                                        std::string word, /// input file keyword
                                        int line, /// input file line
                                        std::string prefix, /// input file prefix
                                        std::string input_file_type) ///input file name
    {

        // check size
        if(u.size()!=3){
        terminaltextcolor(RED);
            std::cerr << "Error: unit-vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must have three values." << std::endl;
            terminaltextcolor(WHITE);
        zlog << zTs() << "Error: unit-vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must have three values." << std::endl;
            err::vexit();
        }

        // Normalise
        double ULength=sqrt(u.at(0)*u.at(0)+u.at(1)*u.at(1)+u.at(2)*u.at(2));

        // Check for correct length unit vector
        if(ULength < 1.0e-9){
        terminaltextcolor(RED);
            std::cerr << "Error: unit-vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be normalisable (possibly all zero)." << std::endl;
            terminaltextcolor(WHITE);
        zlog << zTs() << "Error: unit-vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be normalisable (possibly all zero)." << std::endl;
            err::vexit();
        }
        u.at(0)/=ULength;
        u.at(1)/=ULength;
        u.at(2)/=ULength;

        // Success - input is sane!
        return;

    }

    ///
    /// Function to check for correct 3-component vector and ensure length of 1
    ///-------------------------------------------------------------------------
    ///
    void check_for_valid_three_vector(std::vector<double>& u, /// unit vector
                               std::string word, /// input file keyword
                               int line, /// input file line
                               std::string prefix, /// input file prefix
                               std::string input_file_type) ///input file name
    {

       // check size
       if(u.size()!=3){
    	  terminaltextcolor(RED);
          std::cerr << "Error: vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must have three values." << std::endl;
          terminaltextcolor(WHITE);
    	  zlog << zTs() << "Error: vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must have three values." << std::endl;
          err::vexit();
       }
       // Check for valid range
       if(fabs(u.at(0)) >1.e10){
    	  terminaltextcolor(RED);
          std::cerr << "Error: first element of vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be between +/- 1e10." << std::endl;
          terminaltextcolor(WHITE);
    	  zlog << zTs() << "Error: first element of vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be between +/- 1e10." << std::endl;
          err::vexit();
       }
       if(fabs(u.at(1)) >1.e10){
    	  terminaltextcolor(RED);
          std::cerr << "Error: second element of vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be between +/- 1e10." << std::endl;
          terminaltextcolor(WHITE);
    	  zlog << zTs() << "Error: second element of vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be between +/- 1e10." << std::endl;
          err::vexit();
       }
       if(fabs(u.at(2)) >1.e10){
    	  terminaltextcolor(RED);
          std::cerr << "Error: third element of vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be between +/- 1e10." << std::endl;
          terminaltextcolor(WHITE);
    	  zlog << zTs() << "Error: third element of vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be between +/- 1e10." << std::endl;
          err::vexit();
       }

       // Success - input is sane!
       return;

    }

    ///
    /// Function to check for correct vector with valid values
    ///-------------------------------------------------------------------------
    ///
    void check_for_valid_vector(std::vector<double>& u, /// unit vector
                                std::string word, /// input file keyword
                                int line, /// input file line
                                std::string prefix, /// input file prefix
                                std::string unit, /// unit specified in input file
                                std::string unit_type, /// expected unit type
                                double range_min, /// acceptable minimum value for variable
                                double range_max, /// acceptable maximum value for variable
                                std::string input_file_type, ///input file name
                                std::string range_text) /// customised text
    {

       //---------------------------------------------------------------------------
       // Check for valid unit
       //---------------------------------------------------------------------------

       for(size_t idx=0; idx<u.size(); idx++){

          double value = u.at(idx);

          // Define test unit
       	std::string test_unit_type=unit_type;

       	// Define integer for unit conversion status
       	int convert_status=0;

       	// If no unit given, assume internal, otherwise convert to internal units
       	if(unit.size() != 0) convert_status = units::convert(unit,value,test_unit_type);

       	// Test for valid conversion
       	if(convert_status==EXIT_FAILURE){
       		terminaltextcolor(RED);
       		std::cerr << "Error: Unit \'" << unit << "\' specified on line " << line << " of " << input_file_type << " file is not a valid unit." << std::endl;
       		terminaltextcolor(WHITE);
       		zlog << zTs() << "Error: Unit \'" << unit << "\' specified on line " << line << " of " << input_file_type << " file is not a valid unit." << std::endl;
       		err::vexit();
       	}

       	// Test for change in unit type in case of wrong unit type
       	if(unit_type!=test_unit_type){
       		terminaltextcolor(RED);
       		std::cerr << "Error: Unit \'" << unit << "\' of type \'" << test_unit_type << "\' specified on line " << line << " of " << input_file_type << " is invalid for parameter " << prefix << word << "."<< std::endl;
       		terminaltextcolor(WHITE);
       		zlog << zTs() << "Error: Unit \'" << unit << "\' of type \'" << test_unit_type << "\' specified on line " << line << " of " << input_file_type << " is invalid for parameter " << prefix << word << "."<< std::endl;
       		err::vexit();
       	}

          // Check for valid range
          if((fabs(value)<range_min) || (fabs(value)>range_max)){
             terminaltextcolor(RED);
             std::cerr << "Error: element " << idx+1 << " of vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
             terminaltextcolor(WHITE);
       	   zlog << zTs() << "Error: element " << idx+1 << " of vector variable " << prefix << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
             err::vexit();
          }

          // save value back to array
          u.at(idx) = value;

       }

       // Success - input is sane!
       return;

    }

    void check_for_valid_vector(std::vector<double>& u, /// unit vector
                                std::string word, /// input file keyword
                                int line, /// input file line
                                std::string prefix, /// input file prefix
                                std::string unit, /// unit specified in input file
                                std::string unit_type, /// expected unit type
                                const std::vector <double>& range_min, /// acceptable minimum value for variable
                                const std::vector <double>& range_max, /// acceptable maximum value for variable
                                std::string input_file_type, ///input file name
                                std::string range_text) /// customised text
    {

       //---------------------------------------------------------------------------
       // Check for valid unit
       //---------------------------------------------------------------------------

       for(size_t idx=0; idx<u.size(); ++idx){

          double value = u.at(idx);
          double minvalue = range_min.at(idx);
          double maxvalue = range_max.at(idx);

          // Define test unit
       	std::string test_unit_type=unit_type;

       	// Define integer for unit conversion status
       	int convert_status=0;

       	// If no unit given, assume internal, otherwise convert to internal units
       	if(unit.size() != 0) convert_status = units::convert(unit,value,test_unit_type);

       	// Test for valid conversion
       	if(convert_status==EXIT_FAILURE){
       		terminaltextcolor(RED);
       		std::cerr << "Error: Unit \'" << unit << "\' specified on line " << line << " of " << input_file_type << " file is not a valid unit." << std::endl;
       		terminaltextcolor(WHITE);
       		zlog << zTs() << "Error: Unit \'" << unit << "\' specified on line " << line << " of " << input_file_type << " file is not a valid unit." << std::endl;
       		err::vexit();
       	}

       	// Test for change in unit type in case of wrong unit type
       	if(unit_type!=test_unit_type){
       		terminaltextcolor(RED);
       		std::cerr << "Error: Unit \'" << unit << "\' of type \'" << test_unit_type << "\' specified on line " << line << " of " << input_file_type << " is invalid for parameter " << prefix << word << "."<< std::endl;
       		terminaltextcolor(WHITE);
       		zlog << zTs() << "Error: Unit \'" << unit << "\' of type \'" << test_unit_type << "\' specified on line " << line << " of " << input_file_type << " is invalid for parameter " << prefix << word << "."<< std::endl;
       		err::vexit();
       	}

          // Check for valid range
          if((value<minvalue) || (value>maxvalue)){
             terminaltextcolor(RED);
             std::cerr << "Error: element " << idx+1 << " of vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
             terminaltextcolor(WHITE);
       	   zlog << zTs() << "Error: element " << idx+1 << " of vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
             err::vexit();
          }

          // save value back to array
          u.at(idx) = value;

       }

       // Success - input is sane!
       return;

    }


    ///
    /// Function to check for vector with valid values for bit sequence:
	 /// only -1, 0, 1 values are accepted
    ///-------------------------------------------------------------------------
    ///
	void check_for_valid_bitsequence(std::vector<int>& u, /// vector of integers
                                std::string word, /// input file keyword
                                int line, /// input file line
                                std::string prefix, /// input file prefix
                                int range_min, /// acceptable minimum value for variable
                                int range_max, /// acceptable maximum value for variable
                                std::string input_file_type, ///input file name
                                std::string range_text) /// customised text
   {

      for(size_t idx=0; idx<u.size(); ++idx){

         int value = u.at(idx);
         int minvalue = range_min;
         int maxvalue = range_max;

         // Check for valid range
         if((value<minvalue) || (value>maxvalue)){
            terminaltextcolor(RED);
            std::cerr << "Error: element " << idx+1 << " of vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error: element " << idx+1 << " of vector variable " << prefix << ":" << word << " on line " << line << " of " << input_file_type << " file must be in the range " << range_text << "." << std::endl;
            err::vexit();
         }
         // save value back to array
         u.at(idx) = value;
      }
      // Success - input is sane!
      return;
   }

}
