#ifndef RANDOM_H_
#define RANDOM_H_
#include "mtrand.hpp"
namespace mtrandom
//==========================================================
// Namespace mtrandom
//==========================================================
{
	extern MTRand grnd; // single sequence of random numbers
	extern double gaussian();
}


#endif /*RANDOM_H_*/
