//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2017. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "dipole.hpp"

// dipole module headers
#include "internal.hpp"

namespace dipole{

   //----------------------------------------------
   // function to finalize dipole solver
   //----------------------------------------------
   void finalize(){

      // if fft solver enabled then deallocate memory
      if(dipole::fft) dipole::internal::finalize_fft_solver();

      return;

   }

} // end of dipole namespace
