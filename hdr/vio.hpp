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
#ifndef VIO_H_
#define VIO_H_

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <time.h>
#include <sys/types.h>
#ifdef WIN_COMPILE
  #include <process.h>
  #include <windows.h>
#else
  #include <unistd.h>
#endif

#include "vmpi.hpp"
#include "material.hpp"
#include <iostream>
#include <iomanip>
#include <vector>

// Global Output Streams
extern std::ofstream zinfo;
extern std::ofstream zmag;
extern std::ofstream zgrain;
extern std::ofstream zlog;
extern std::ofstream dp_fields;

enum textcolor {
   RED=31,
   GREEN=32,
   YELLOW=33,
   BLUE=34,
   PURPLE=35,
   WHITE=0
};

void terminaltextcolor(enum textcolor );

// global timestap function
extern std::string zTs();

namespace vin{
   extern int read(std::string const);

   extern void check_for_valid_value(double& value, std::string word, int line, std::string prefix, std::string unit, std::string unit_type,
                                     double range_min, double range_max, std::string input_file_type, std::string range_text);

   extern void check_for_valid_positive_value(double& value, std::string word, int line, std::string prefix, std::string unit, std::string unit_type,
                                              double range_min, double range_max, std::string input_file_type, std::string range_text);

   extern void check_for_valid_int(int& value, std::string word, int line, std::string prefix, int range_min, int range_max,
                                   std::string input_file_type, std::string range_text);

   extern void check_for_valid_int(  unsigned int& value, std::string word, int line, std::string prefix, unsigned int range_min,
                              unsigned int range_max, std::string input_file_type, std::string range_text);

   extern void check_for_valid_int(  uint64_t& value, std::string word, int line, std::string prefix,
                                     uint64_t range_min, uint64_t range_max, std::string input_file_type, std::string range_text);

   extern bool check_for_valid_bool( std::string value, std::string word, int line, std::string prefix, std::string input_file_type);

   extern void check_for_valid_unit_vector(std::vector<double>& u, std::string word, int line, std::string prefix, std::string input_file_type);

   extern void check_for_valid_three_vector(std::vector<double>& u, std::string word, int line, std::string prefix, std::string input_file_type);

   extern void check_for_valid_vector(std::vector<double>& u, std::string word, int line, std::string prefix, std::string unit, std::string unit_type,
                                      double range_min, double range_max, std::string input_file_type, std::string range_text);

   extern void check_for_valid_vector(std::vector<double>& u, std::string word, int line, std::string prefix, std::string unit, std::string unit_type,
                                      const std::vector<double>& range_min, const std::vector<double>& range_max, std::string input_file_type, std::string range_text);

	extern void check_for_valid_bitsequence(std::vector<int>& u, std::string word, int line, std::string prefix, const int range_min, 
													const int range_max, std::string input_file_type, std::string range_text); 

   extern std::vector<double> doubles_from_string(std::string value);
   extern std::vector<int> integers_from_string(std::string value);

   // function to read file on master process and return a std::string of its contents
   extern std::string get_string(std::string const filename, std::string source_file_name, int line);

   // simple functions to extract variables from strings
   extern uint64_t str_to_uint64(std::string input_str);
   extern double str_to_double(std::string input_str);

   extern std::vector<mp::materials_t> read_material;

}

namespace vout{

   extern bool custom_precision; // enable user selectable precision for data output
   extern unsigned int precision; // variable to control output precision (digits)
   extern bool fixed; // fixed precision output
   extern bool header_option; // column headers
	extern std::vector<unsigned int> file_output_list;
	extern std::vector<unsigned int> screen_output_list;

   extern int output_rate;

   extern bool gnuplot_array_format;

	//extern bool output_povray;
	//extern int output_povray_rate;

	//extern bool output_povray_cells;
	//extern int output_povray_cells_rate;

	extern void data();
	extern void zLogTsInit(std::string);
    void output_switch(std::ostream&, unsigned int);
    extern void write_out(std::ostream&,std::vector<unsigned int>&);
	//extern int pov_file();

	void redirect(std::ostream& strm, std::string filename);
	void nullify(std::ostream& strm);

    extern int fw_size;
    extern int fw_size_int;
    extern int max_header;
    
    extern std::string output_file_name;

//class that creates an object which acts like an output
//stream but delivers fixed width output separated by
//tabs
//as a result, outputs that are part of one column
//should be concatenated before output, so as to
//prevent splitting out into multiple columns.
//to use you should use the output stream you
//would normally be using as an argument during construction
//and the width of your columns.
//     e.g.
//     std::ostringstream res;
//     vout::fixed_width_output result(res,vout::fw_size);
//you can then use the <<operator to send output to this
//stream, but formatted.
//     e.g.
//     result << "ID" + std::to_string(mask_id);

class fixed_width_output{
  private:
    int width; // the width of each output
    std::ostringstream& stream_obj; // the initial stream object
  public:
    //constructor, calls constructors of width and stream_obj
    //                                          :-> member initialization list
    fixed_width_output(std::ostringstream& obj, int w): width(w),stream_obj(obj) {};

    template<typename T> // sets up a template for using the operator<<

    // defines a function which returns a pointer to the fixed... object
    // takes one of type T as input.
    fixed_width_output& operator<<(const T& output){
      //sends the formatted output to a stream_obj
      stream_obj <<std::left<<std::setw(width) << output <<"\t";

      //returns the object (dereferenced from the pointer)
      return *this;
    }

    // specialises the function, for when the input is an output stream
    // which is being operated on, such as using <<std::endl;
    fixed_width_output& operator<<(std::ostringstream& (*func)(std::ostringstream&)){
        func(stream_obj);
        return *this;
    };
    std::string str(){
        return stream_obj.str();
    };
};

}

// Checkpoint load/save functions
void load_checkpoint();
void save_checkpoint();

namespace vio{
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

}

#endif /*VIO_H_*/
