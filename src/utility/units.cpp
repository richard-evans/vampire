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
///
/// @file
/// @brief Contains the units namespace for program parameters and unit conversion. 
///
/// @section notes Implementation Notes
/// None.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    18/01/2010
/// @internal
///	Created:		18/01/2010
///	Revision:	  ---
///=====================================================================================
///
// Headers
#include "errors.hpp"
#include "units.hpp"
#include "vio.hpp"

#include <iostream>
#include <cstdlib>

/// @namespace units
/// @brief Contains program parameters and functions for unit conversion.
/// 
/// @internal
///=====================================================================================
///
namespace units {
	
	const int max_units=43;

	const double pi=3.14;
	//const double bohr_magneton=7.0;
	bool initialised=false;
	// Generic unit names and conversion factors
	std::string unit[max_units];
	std::string type[max_units];
	double conversion[max_units];

/// @brief Initialises system units for conversion.
///
/// @details Sets units::initialised=true. Example usage
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    18/01/2010
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		18/01/2010
///	Revision:	  ---
///=====================================================================================
///
int init(){
		// check calling of routine if error checking is activated
		if(err::check==true){
			std::cout << "units::init has been called " << std::endl;
		}
		using units::unit;
		using units::conversion;
		// Distances
		unit[0]="";				conversion[0]=1.0; 			type[0]="none";		// Default (No conversion)
		unit[1]="A";			conversion[1]=1.0; 			type[1]="length";		// Angstrom (Internal)
		unit[2]="m";			conversion[2]=1.0E10; 		type[2]="length";		// Metres
		unit[3]="cm";			conversion[3]=1.0E8; 		type[3]="length";		// Centimetres
		unit[4]="mm";			conversion[4]=1.0E7;	 		type[4]="length";		// Millimetres
		unit[5]="um";			conversion[5]=1.0E4; 		type[5]="length";		// Micrometres
		unit[6]="nm";			conversion[6]=1.0E1; 		type[6]="length";		// Nanometres
		unit[7]="pm";			conversion[7]=1.0E-3;	 	type[7]="length";		// Picometres
		unit[8]="in";			conversion[8]=2.54E8; 		type[8]="length";		// Imperial Inches
		unit[9]="feet";		conversion[9]=3.12E9; 		type[9]="length";		// Imperial Feet
		// Energies
		unit[10]="J";			conversion[10]=1.0; 					type[10]="energy";		// Joules
		unit[11]="eV";			conversion[11]=1.602176487e-19; 	type[11]="energy";		// Electron Volts
		unit[12]="meV";		conversion[12]=1.602176487e-22; 	type[12]="energy";		// Millielectron Volts
		unit[13]="erg";		conversion[13]=1E-7; 				type[13]="energy";		// Ergs
		unit[14]="Ryd";		conversion[14]=2.17987208E-18;	type[14]="energy";		// Rydbergs
		unit[15]="mRyd";		conversion[15]=2.17987208E-21; 	type[15]="energy";		// MilliRydbergs
		unit[16]="Ht";			conversion[16]=4.35974417E-18;	type[16]="energy";		// Hartrees
		// Spin moment
		unit[17]="J/T";		conversion[17]=1.0; 					type[17]="moment";		// Joules/Tesla
		unit[18]="Amm";		conversion[18]=1.0; 					type[18]="moment";		// Amps metres squared
		unit[19]="erg/G";		conversion[19]=1.0E-3; 				type[19]="moment";		// Erg/Gauss
		unit[20]="abAcmcm";	conversion[20]=1.0E-3; 				type[20]="moment";		// Abampere centimetre squared
		unit[21]="muB";		conversion[21]=9.27400915e-24; 	type[21]="moment";		// Bohr Magnetons
		unit[22]="eV/T";		conversion[22]=1.602176487e-19;	type[22]="moment";		// Electron volts/Tesla
		unit[23]="erg/Oe";	conversion[23]=1.0E-3; 				type[23]="moment";		// Erg/Oersted
		// Magnetisation
		unit[24]="A/m";		conversion[24]=2.17987208E-18;	type[24]="magnetisation";	// Amps/metre
		unit[25]="emu/cc";	conversion[25]=2.17987208E-21; 	type[25]="magnetisation";	// emu/cubic cm
		unit[26]="J/T/AAA";	conversion[26]=1.0;					type[26]="magnetisation";	// Joule/Tesla/Angstrom cubed
		// Anisotropy
		unit[27]="J/atom";	conversion[27]=1.0; 					type[27]="anisotropy";		// Joules/Atom
		unit[28]="J/mmm";		conversion[28]=1.0E-30; 			type[28]="anisotropy";		// Joules/metres cubed
		unit[29]="erg/cc";	conversion[29]=1.0E-31; 			type[29]="anisotropy";		// Erg/cc
		// Field
		unit[30]="T";			conversion[30]=1.0; 					type[30]="field";		// Tesla
		unit[31]="mT";			conversion[31]=1.0E-3;		 		type[31]="field";		// milliTesla
		unit[32]="uT";			conversion[32]=1.0E-6;		 		type[32]="field";		// microTesla
		unit[33]="Oe";			conversion[33]=1.0E-4;				type[33]="field";		// Oersted
		unit[34]="kOe";		conversion[34]=1.0E-1; 				type[34]="field";		// kilo Oersted
		// Time
		unit[35]="s";			conversion[35]=1.0;					type[35]="time"; // seconds
		unit[36]="ms";			conversion[36]=1.0E-3;				type[36]="time"; // milliseconds
		unit[37]="us";			conversion[37]=1.0E-6;				type[37]="time"; // microseconds
		unit[38]="ns";			conversion[38]=1.0E-9;				type[38]="time"; // nanoseconds
		unit[39]="ps";			conversion[39]=1.0E-12;				type[39]="time"; // picoseconds
		unit[40]="fs";			conversion[40]=1.0E-15;				type[40]="time"; // femtoseconds
		unit[41]="as";			conversion[41]=1.0E-18;				type[41]="time"; // attoseconds
		unit[42]="zs";			conversion[42]=1.0E-21;				type[42]="time"; // zeptoseconds

      // temperature C, F, K; angles degrees, rad, mrad;
		// Set initialised flag
		units::initialised=true;
		
		return EXIT_SUCCESS;
	}
	
/// @brief Converts external units to internal units.
///
/// @details Example usage units::convert(string "nm", double& var, string& unit_type)
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    18/01/2010
///
/// @param[in] input_unit name of unit to be converted
/// @param[in] value variable to be converted
/// @param[out] type unit type, eg length, volume 
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		18/01/2010
///	Revision:	  ---
///=====================================================================================
///
	int convert(std::string input_unit, double& value, std::string& type){

		// Populate unit array;
		if(units::initialised==false){
			units::init();
		}
		
		// loop over all possible units
		for(int i=0;i<max_units;i++){
			if(input_unit==units::unit[i]){
				// Convert unit
				value*=conversion[i];
				// Set unit type
				type=units::type[i];
				// return
				return EXIT_SUCCESS;
			}
		}

		return EXIT_FAILURE;

	}

/// @brief Converts array of external units to internal units.
///
/// @details Example usage units::convert(string "nm", std::vector<double>& var, string& unit_type)
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    21/07/2012
///
/// @param[in] input_unit name of unit to be converted
/// @param[in] value vector of variables to be converted
/// @param[out] type unit type, eg length, volume 
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		21/07/2012
///	Revision:	  ---
///=====================================================================================
///
	void convert(std::string input_unit, std::vector<double>& value, std::string& type){

		// Populate unit array;
		if(units::initialised==false){
			units::init();
		}
		
		// loop over all possible units
		for(int i=0;i<max_units;i++){
			if(input_unit==units::unit[i]){
				// Convert unit
				for(unsigned int idx=0; idx < value.size(); idx++) value.at(idx)*=conversion[i];
				// Set unit type
				type=units::type[i];

				return;
			}
		}

		// Error if unit not found
		std::cerr << "Error during unit conversion - unit \'"<< input_unit << "\' not found" << std::endl;
		zlog << zTs() << "Error during unit conversion - unit \'"<< input_unit << "\' not found" << std::endl;
		err::vexit();

	}
	
/// @brief Reverts internal units to external units.
///
/// @details Example usage @code units::revert(string "nm", double& var, string& unit_type) /@code
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    18/01/2010
///
/// @param[in] output_unit name of unit to be converted
/// @param[in] value variable to be converted
/// @param[out] type unit type, eg length, volume 
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		18/01/2010
///	Revision:	  ---
///=====================================================================================
///
	int revert(std::string output_unit, double& value, std::string& type){	

		// Populate unit array;
		if(units::initialised==false){
			units::init();
		}
		
		// loop over all possible units
		for(int i=0;i<max_units;i++){
			if(output_unit==units::unit[i]){
				// Convert unit
				value/=conversion[i];
				// Set unit type
				type=units::type[i];
				// return
				return EXIT_SUCCESS;
			}
		}

		// Error if unit not found
		std::cerr << "Error during unit reversion - unit \'"<< output_unit << "\' not found" << std::endl;
		err::vexit();
		
		return EXIT_SUCCESS;
		
	}

} // End of namespace

