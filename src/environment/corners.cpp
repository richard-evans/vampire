//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2018. All rights reserved.
//
//   Email: sarah.jenkins@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
// C++ standard library headers
#include <string>
#include <cmath>
#include <cstdlib>
#include <iostream>
// C library headers
#include <fenv.h>
#include <signal.h>
#include <math.h>
// Vampire headers
#include "cells.hpp" // needed for cells::cell_id_array but to be removed
#include "dipole.hpp"
#include "vio.hpp"
#include "vutil.hpp"
#include "atoms.hpp"
// Vampire headers
#include "hierarchical.hpp"


// hierarchical module headers
#include "internal.hpp"

using namespace std;

namespace environment{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

    std::vector < std::vector < double> > calculate_corners(double x, double y, double z, double cell_size_x, double cell_size_y, double cell_size_z){

      std::vector < std::vector < double> >corners;
      corners.resize(8);
      for (int i = 0; i < 8; ++i)
        corners[i].resize(3,0.0);

      double min_x = x - cell_size_x/2.0;
      double min_y = y - cell_size_y/2.0;
      double min_z = z - cell_size_z/2.0;

      double max_x = x + cell_size_x/2.0;
      double max_y = y + cell_size_y/2.0;
      double max_z = z + cell_size_z/2.0;

      corners[0][0] = min_x;
      corners[1][0] = min_x;
      corners[2][0] = min_x;
      corners[3][0] = min_x;
      corners[4][0] = max_x;
      corners[5][0] = max_x;
      corners[6][0] = max_x;
      corners[7][0] = max_x;

      corners[0][1] = min_y;
      corners[1][1] = min_y;
      corners[2][1] = max_y;
      corners[3][1] = max_y;
      corners[4][1] = min_y;
      corners[5][1] = min_y;
      corners[6][1] = max_y;
      corners[7][1] = max_y;

      corners[0][2] = min_z;
      corners[1][2] = max_z;
      corners[2][2] = min_z;
      corners[3][2] = max_z;
      corners[4][2] = min_z;
      corners[5][2] = max_z;
      corners[6][2] = min_z;
      corners[7][2] = max_z;

      return corners;

    }
}
}
