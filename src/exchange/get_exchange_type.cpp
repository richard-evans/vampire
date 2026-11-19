//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "exchange.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

   //------------------------------------------------------------------------------
   // Function to return exchange type
   //------------------------------------------------------------------------------
   unsigned int get_exchange_type(){
      return internal::exchange_type;
   }

   //------------------------------------------------------------------------------
   // Function to return whether the exchange set-up needs the vectors i->j of the
   // bilinear neighbour list
   //------------------------------------------------------------------------------
   bool neighbour_vectors_needed(){
      return internal::enable_dmi || internal::enable_kitaev || internal::enable_fourspin;
   }

} // end of exchange namespace
