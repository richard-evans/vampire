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
#ifndef LLG_H_
#define LLG_H_
/// Header file for LLG namespace
namespace LLG_arrays{
	
//==========================================================
// Namespace to store persistant LLG integration arrays
//==========================================================

	extern std::vector <double> x_euler_array;	
	extern std::vector <double> y_euler_array;	
	extern std::vector <double> z_euler_array;

	extern std::vector <double> x_heun_array;	
	extern std::vector <double> y_heun_array;	
	extern std::vector <double> z_heun_array;

	extern std::vector <double> x_spin_storage_array;	
	extern std::vector <double> y_spin_storage_array;	
	extern std::vector <double> z_spin_storage_array;

	extern std::vector <double> x_initial_spin_array;	
	extern std::vector <double> y_initial_spin_array;	
	extern std::vector <double> z_initial_spin_array;

	extern bool LLG_set;

}
#endif /*LLG_H_*/
