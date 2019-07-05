//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <iostream>

// C library headers
#include <fenv.h>
#include <signal.h>

// Vampire headers
#include "cells.hpp" // needed for dp::cell_id_array but to be removed
#include "dipole.hpp"
#include "vio.hpp"
#include "vutil.hpp"
#include "math.h"
#include "errors.hpp"
#include "atoms.hpp"

// dipole module headers
#include "../cells/internal.hpp"
#include "internal.hpp"

// Vampire headers
#include "hierarchical.hpp"


// hierarchical module headers
#include "internal.hpp"

namespace ha = hierarchical::internal;


namespace hierarchical{
  namespace internal{

   //-----------------------------------------------------------------------------
   // Function for calculate magnetisation in cells
   //-----------------------------------------------------------------------------
    int hierarchical_mag(){


    cells::mag();

   //std::cout <<"CELLS:\t" << num_zero_level_cells << "\t" << cells::num_cells <<std::endl;

        for (int cell = 0; cell < cells::num_cells; cell ++){

       ha::mag_array_x[cell] = cells::mag_array_x[cell];
       ha::mag_array_y[cell] = cells::mag_array_y[cell];
       ha::mag_array_z[cell] = cells::mag_array_z[cell];
       //std::cout << cell << '\t' << cells::mag_array_x[cell]*1.0/9.27400915e-24 <<std::endl;
    }

    for (int level = 1; level < ha::num_levels; level ++ ){

      int start = ha::cells_level_start_index[level];
      int end   = ha::cells_level_end_index[level];

        for (int cell = start; cell < end; cell++){

          int start_cell_in_cell = ha::cells_in_cells_start_index[cell];
          int end_cell_in_cell = ha::cells_in_cells_end_index[cell];

          ha::mag_array_x[cell] = 0;
          ha::mag_array_y[cell] = 0;
          ha::mag_array_z[cell] = 0;

          for (int cell_in_cell = start_cell_in_cell; cell_in_cell < end_cell_in_cell; cell_in_cell++){

            ha::mag_array_x[cell] += ha::mag_array_x[cell_in_cell];
            ha::mag_array_y[cell] += ha::mag_array_y[cell_in_cell];
            ha::mag_array_z[cell] += ha::mag_array_z[cell_in_cell];

          }
        }
      }
      return 0;
    }
  }
}
