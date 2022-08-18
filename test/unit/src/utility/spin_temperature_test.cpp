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
#include <vector>

// include header for test functions
#include "sld.hpp"

#include "constants.hpp"
#include "material.hpp"


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



   namespace utility{
      /*
      compute_spin_temperature(const int start_index, // first atom for exchange interactions to be calculated
                  const int end_index,
                  const std::vector<int>& type_array, // type for atom
                  std::vector<double>& x_spin_array, // coord vectors for atoms
                  std::vector<double>& y_spin_array,
                  std::vector<double>& z_spin_array,
                  std::vector<double>& fields_array_x, //  vectors for fields
                  std::vector<double>& fields_array_y,
                  std::vector<double>& fields_array_z)*/
                  //	int convert(std::string input_unit, double& value, std::string& type)
      int spin_temperature_test(const int test_start_index, const int test_end_index, const std::vector<int>& test_type_array, std::vector<double>& test_x_spin_array, std::vector<double>& test_y_spin_array, std::vector<double>& test_z_spin_array, std::vector<double>& test_fields_array_x, std::vector<double>& test_fields_array_y, std::vector<double>& test_fields_array_z, std::vector<double>& test_mu_s_array, const double expected_value, const double precision){

         int ec = 0; // error count increment

         // call function to be tested
         //convert is from vampire
         double test_value=sld::compute_spin_temperature(test_start_index, test_end_index, test_type_array,  test_x_spin_array, test_y_spin_array, test_z_spin_array,  test_fields_array_x,  test_fields_array_y,  test_fields_array_z, test_mu_s_array);
         // check for numerical error
         ec += ut::floaterror(test_value, expected_value, precision, "sld::compute_spin_temperature");

         return ec;
      }

//------------------------------------------------------------------------------
// Function to test utility module functions
//------------------------------------------------------------------------------
int test_spin_temperature(const bool verbose){

   //extern int convert(std::string input_unit, double& value, std::string& type);
	//extern void convert(std::string input_unit, std::vector<double>& value, std::string& type);

   // specify numerical precision required for all tests
   const double precision = 0.00001;
   const double expected_value=0.0;


   int ec = 0; // error counter

   //-----------------------------------------------------
   // Testing units::convert(std::string input_unit, double& value, std::string& type);
   //-----------------------------------------------------
   //int spin_temperature_test(const std::int test_start_index, const std::int test_start_index, const std::int test_end_index, const std::vector<int>& test_type_array, std::vector<double>& test_x_spin_array, std::vector<double>& test_y_spin_array, std::vector<double>& test_z_spin_array, std::vector<double>& test_fields_array_x, std::vector<double>& test_fields_array_y, std::vector<double>& test_fields_array_z, const double expected_value, const double precision){
   std::vector<double> x_spin_array;
   x_spin_array.resize(2,0);

   std::vector<double> y_spin_array;
   y_spin_array.resize(2,0);

   std::vector<double> z_spin_array;
   z_spin_array.resize(2,1);

   std::vector<double> fields_array_x;
   fields_array_x.resize(2,0);

   std::vector<double> fields_array_y;
   fields_array_y.resize(2,0);

   std::vector<double> fields_array_z;
   fields_array_z.resize(2,100);

   std::vector<int> type_array;
   type_array.resize(2,0);

   std::vector<double> mu_s_array;
   mu_s_array.resize(2,2.22);

   ec += spin_temperature_test(0,2,type_array,x_spin_array,y_spin_array,z_spin_array,fields_array_x,fields_array_y,fields_array_z, mu_s_array, expected_value,precision);

   return ec;

}

}
}
