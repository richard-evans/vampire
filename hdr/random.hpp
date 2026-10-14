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
#ifndef RANDOM_H_
#define RANDOM_H_
#include "mtrand.hpp"
namespace mtrandom
//==========================================================
// Namespace mtrandom
//==========================================================
{
   extern MTRand grnd; /// single sequence of random numbers
   extern double gaussian();
   extern double gaussianc(MTRand&);

   extern int voronoi_seed;
   extern int integration_seed;
}


#endif /*RANDOM_H_*/
