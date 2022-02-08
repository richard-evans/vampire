//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef EXCHANGE_TYPES_H_
#define EXCHANGE_TYPES_H_

//--------------------------------------------------------------------------------
// Namespace for variables and functions for exchange module
//--------------------------------------------------------------------------------
namespace exchange{

   //---------------------------------------------------------------------------
   // Enumerated list of available exchange types
   //---------------------------------------------------------------------------
   enum exchange_t { isotropic = 0, // isotropic exchange interactions
                     vectorial = 1, // vector exchange Jxx, Jyy, Jzz
                     tensorial = 2 // tensor exchange Jxx, Jxy ... Jzz
   };

}

#endif //EXCHANGE_TYPES_H_
