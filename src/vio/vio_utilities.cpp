//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify 
//  it under the terms of the GNU General Public License as published by 
//  the Free Software Foundation; either version 2 of the License, or 
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but 
//  WITHOUT ANY WARRANTY; without even the implied warranty of 
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License 
//  along with this program; if not, write to the Free Software Foundation, 
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
//  Utilities for file i/o operations
//
// (c) R F L Evans 2013
//
// c++ headers
#include <string>
#include <sstream>
#include <vector>

// vampire headers
#include "./vio.ihpp"

namespace vio{

// Function to extract comma-separated double precision numbers from string
//---------------------------------------------------------------------------
std::vector<double> get_doubles_from_string(std::string value){

   // array for storing variables
   std::vector<double> array(0);
   
   // set source for ss
   std::istringstream source(value);

   // double variable to store values
   double temp = 0.0;
   
   // string to store text
   std::string field;

   // loop over all comma separated values
   while(getline(source,field,',')){
      
      // convert string to ss
      std::stringstream fs(field);

      // read in variable
      fs >> temp;
      
      // push data value back to array
      array.push_back(temp);
      
   }
   
   // return values to calling function
   return array;

}

// Functions converting string to value
// Double
//----------------------------------------------------------------------------
int string_to_value(std::string value_string, double& value){
   value=atof(value_string.c_str());
   return 0;
}

// Int
//----------------------------------------------------------------------------
int string_to_value(std::string value_string, int& value){ 
   value=atoi(value_string.c_str());
   return 0;
}

int string_to_value(std::string value_string, bool& value){
      std::string strue="true";
      std::string sfalse="false";
      if(value_string==strue){
         value=true;
         return 0; // success
      }
      else if(value_string==sfalse){
         value=false;
         return 0; // success
      }
      else return 1; // error 
}

// Vectors
int string_to_value(std::string value_string, vector_t& value){

   std::vector<double> values = get_doubles_from_string(value_string);
   if(values.size()!=3){
      return 2; // error
   }
   else{
      value.x=values.at(0);
      value.y=values.at(1);
      value.z=values.at(2);
      return 0; // success
   }
}
   
// Unit Vectors
int string_to_value(std::string value_string, uvector_t& value){

   std::vector<double> values = get_doubles_from_string(value_string);
   if(values.size()!=3){
      return 2; // error
   }
   else{
      double ulength=sqrt(values.at(0)*values.at(0)+values.at(1)*values.at(1)+values.at(2)*values.at(2));

      // Check for correct length unit vector
      if(ulength < 1.0e-9) return 3; // error

      values.at(0)/=ulength;
      values.at(1)/=ulength;
      values.at(2)/=ulength;
      
      value.x=values.at(0);
      value.y=values.at(1);
      value.z=values.at(2);
      return 0; // success
   }
}

// Tensors
int string_to_value(std::string value_string, tensor_t& value){

   std::vector<double> values = get_doubles_from_string(value_string);
   if(values.size()!=9){
      return 4; // error
   }
   else{
      value.xx=values.at(0);
      value.xy=values.at(1);
      value.xz=values.at(2);
      value.yx=values.at(3);
      value.yy=values.at(4);
      value.yz=values.at(5);
      value.zx=values.at(6);
      value.zy=values.at(7);
      value.zz=values.at(8);
      return 0; // success
   }
}

}  // end of namespace vio
