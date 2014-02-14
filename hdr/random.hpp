//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify 
//  it under the terms of the GNU General Public License as published by 
//  the Free Software Foundation; either version 2 of the License, or 
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but 
//  WITHOUT ANY WARRANTY; without even the implied warranty of 
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License 
//  along with this program; if not, write to the Free Software Foundation, 
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
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
