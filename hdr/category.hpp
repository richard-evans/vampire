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

#ifndef CATEGORY_H_
#define CATEGORY_H_

#include <vector>

namespace cat{

   class category_t {
   public:

      // list of atoms in each category
      std::vector <int> category_atoms;

      // magnetisation variables
      double magm;
      double mean_magm;
      double mx;
      double my;
      double mz;

      // energy and torque variables
      double torque;
      double mean_torque;
      double energy;

      // constraint variables
      double theta;
      double phi;

      // spin temperature variables
      double spin_temp;
      double mean_spin_temp;
      double SxH2;
      double SH;

      // member functions
      category_t(); /// constructor

   };

   extern std::vector <category_t> category;


} // End of namespace cat

#endif //CATEGORY_H_
