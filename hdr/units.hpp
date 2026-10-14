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
#ifndef UNITS_H_
#define UNITS_H_

#include <string>
#include <vector>

namespace units{

   extern const double pi;

   // Holds all the data for one unit.
   // Keeping these three values in a struct means they always travel together —
   // there is no risk of the name, type, and conversion factor getting out of
   // step with each other (which was easy to do with the old parallel arrays).
   struct unit_t {
      std::string name;        // string the user writes in the input file, e.g. "nm"
      std::string type;        // physical quantity this unit measures, e.g. "length"
      double      conversion;  // multiply the user's value by this to get internal units
   };

   // conversion functions — public API is unchanged; all existing callers still compile
   extern int  convert(std::string input_unit, double& value, std::string& type);
   extern void convert(std::string input_unit, std::vector<double>& value, std::string& type);
   extern int  revert (std::string output_unit, double& value, std::string& type);

}


#endif /* UNITS_H_ */
