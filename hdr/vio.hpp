//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
#ifndef VIO_H_
#define VIO_H_

#include <cstdint>
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

   extern void check_for_valid_int(  std::uint64_t& value, std::string word, int line, std::string prefix,
                                     std::uint64_t range_min, std::uint64_t range_max, std::string input_file_type, std::string range_text);

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
   extern std::uint64_t str_to_uint64(std::string input_str);
   extern int str_to_int(std::string input_str);
   extern double str_to_double(std::string input_str);

}

namespace vout{

   //---------------------------------------------------------------------------
   // Enumerated list of data columns for output: and screen: keywords
   //---------------------------------------------------------------------------
   enum class output_t {
      // simulation time, temperature and fields
      time_steps,
      real_time,
      temperature,
      phonon_temperature,
      material_temperature,
      applied_field_strength,
      applied_field_unit_vector,
      applied_field_alignment,
      material_applied_field_strength,
      material_applied_field_alignment,
      fmr_field_strength,
      material_fmr_field_strength,
      fractional_electric_field_strength,

      // system magnetisation
      magnetisation,
      magnetisation_length,
      mean_magnetisation_length,
      mean_magnetisation,

      // material magnetisation
      material_magnetisation,
      material_mean_magnetisation_length,
      material_mean_magnetisation,

      // height magnetisation
      height_magnetisation,
      height_magnetisation_normalised,
      mean_height_magnetisation,
      mean_height_magnetisation_length,
      material_height_magnetisation,
      material_height_magnetisation_normalised,

      // torque and constraint angles
      total_torque,
      mean_total_torque,
      material_torque,
      material_mean_torque,
      constraint_phi,
      constraint_theta,
      material_constraint_phi,
      material_constraint_theta,

      // susceptibility, specific heat and fluctuations
      mean_susceptibility,
      material_mean_susceptibility,
      mean_specific_heat,
      material_mean_specific_heat,
      binder_cumulant,
      material_binder_cumulant,
      material_standard_deviation,

      // energies
      total_energy,
      mean_total_energy,
      material_total_energy,
      material_mean_total_energy,
      anisotropy_energy,
      mean_anisotropy_energy,
      exchange_energy,
      mean_exchange_energy,
      applied_field_energy,
      mean_applied_field_energy,
      magnetostatic_energy,
      mean_magnetostatic_energy,
      dmi_energy,
      mean_dmi_energy,
      material_dmi_energy,
      material_mean_dmi_energy,

      // spin temperature and spin length
      spin_temperature,
      mean_spin_temperature,
      material_spin_temperature,
      material_mean_spin_temperature,
      mean_spin_length,
      material_mean_spin_length,
      mean_height_spin_length,

      // spin-lattice dynamics
      potential_energy,
      kinetic_energy,
      sld_exchange_energy,
      sld_coupling_energy,
      sld_total_energy,
      mean_coupling_field_vec,
      lattice_temperature,
      mean_lattice_temperature,
      material_lattice_temperature,
      material_mean_lattice_temperature,

      // transport, micromagnetic and domain walls
      resistance,
      current,
      magneto_resistance,
      lfa_ms,
      domain_wall_centre,
      track_pos_x,
      track_pos_z,

      // diagnostics
      mpi_timings
   };

   extern bool custom_precision; // enable user selectable precision for data output
   extern unsigned int precision; // variable to control output precision (digits)
   extern bool fixed; // fixed precision output
   extern bool header_option; // column headers
   extern std::vector<output_t> file_output_list;
   extern std::vector<output_t> screen_output_list;

   extern int output_rate;

   extern bool gnuplot_array_format;

   //extern bool output_povray;
   //extern int output_povray_rate;

   //extern bool output_povray_cells;
   //extern int output_povray_cells_rate;

   extern void data();
   extern void zLogTsInit(std::string);
    extern void write_out(std::ostream&,std::vector<output_t>&);
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

namespace vio{
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);
}

#endif /*VIO_H_*/
