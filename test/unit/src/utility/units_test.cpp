//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <cmath>

// include header for test functions
#include "units.hpp"

namespace vmpi{
   int my_rank = 0;
}

namespace ut{

int floaterror(const double value, const double expected_value, const double precision, const std::string function){
   double error = 0.0;

   // check for division by zero
   if(fabs(expected_value) > 0.0) error = value / expected_value;

   // check for both zero
   if(value == expected_value) return 0; // exact floating point comparison

   if(fabs(1.0-error) < precision) return 0;
   else{
      std::cout << "FAIL: Floating point error in test of function " << function << ": value " << value << " should be the same as " << expected_value << " to precision " << precision << std::endl;
      return 1;
   }
}

int stringerror(const std::string value, const std::string expected_value, const std::string function){

   if(value == expected_value) return 0;
   else{
      if(value == "") std::cout << "FAIL: String error in test of function " << function << ": output string is empty but should be the same as " << expected_value << std::endl;
      else            std::cout << "FAIL: String error in test of function " << function << ": string " << value << " should be the same as " << expected_value << std::endl;
      return 1;
   }
}

   namespace utility{

      int convert_test(const std::string test_unit, double test_value, const double expected_value, std::string test_unit_type, const std::string expected_unit_type, const double precision){

         int ec = 0; // error count increment

         // call function to be tested
         units::convert(test_unit, test_value, test_unit_type);

         // check for numerical error
         ec += ut::floaterror(test_value, expected_value, precision, "units::convert");
         // check for string error
         ec += ut::stringerror(test_unit_type, expected_unit_type, "units::convert");

         return ec;
      }

//------------------------------------------------------------------------------
// Function to test utility module functions
//------------------------------------------------------------------------------
int test_units(const bool verbose){

   //extern int convert(std::string input_unit, double& value, std::string& type);
	//extern void convert(std::string input_unit, std::vector<double>& value, std::string& type);

   // specify numerical precision required for all tests
   const double precision = 0.00001;


   int ec = 0; // error counter

   //-----------------------------------------------------------------------------------------
   // Testing units::convert(std::string input_unit, double& value, std::string& type);
   //-----------------------------------------------------------------------------------------
   ec += convert_test("meV", 1.0, 1.602176634e-22, "","energy",precision);
   ec += convert_test("T",   1.0, 1.0,             "","field",precision);
   ec += convert_test("A/m", 1.0, 1.0,             "","magnetisation",precision);
   ec += convert_test("",    1.0, 1.0,             "","none",precision); // test no unit


   return ec;

}

}
}
