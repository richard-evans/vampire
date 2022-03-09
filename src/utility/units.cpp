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
#include <cmath>
#include <cstdlib>

/// @namespace units
/// @brief Contains program parameters and functions for unit conversion.
///
/// @internal
///=====================================================================================
///
namespace units {

	const int max_units=71;

	const double pi=M_PI;
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
		unit[24]="A/m";		conversion[24]=1e-30;				type[24]="magnetisation";	// Amps/metre
		unit[25]="emu/cc";	conversion[25]=1e-27;			 	type[25]="magnetisation";	// emu/cubic cm
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
      // New
      unit[43]="zJ";			conversion[43]=1.0E-21; 			type[43]="energy";		// zeptoJoules
      unit[44]="Ohm-m";		conversion[44]=1.0; 			      type[44]="resistivity";		// Ohm metres

      unit[45]="V";		   conversion[45]=1.0; 			      type[45]="potential";		// volts
      unit[46]="mV";		   conversion[46]=1.0e-3; 			   type[46]="potential";		// Millivolts
      unit[47]="uV";	    	conversion[47]=1.0e-6; 			   type[47]="potential";		// Microvolts
      unit[48]="nV";		   conversion[48]=1.0e-9; 			   type[48]="potential";		// Nanovolts

		unit[49]="Hz";		   conversion[49]=1.0; 			      type[49]="frequency";		// Hz
      unit[50]="kHz";		conversion[50]=1.0e3; 			   type[50]="frequency";		// kHz
      unit[51]="MHz";	   conversion[51]=1.0e6; 			   type[51]="frequency";		// MHz
      unit[52]="GHz";		conversion[52]=1.0e9; 			   type[52]="frequency";		// GHz
      unit[53]="THz";		conversion[53]=1.0e12; 			   type[53]="frequency";		// THz

		// exchange (Internal unit J/Angstrom)
		unit[54]="J/m";		conversion[54]=1.0E-10; 			type[54]="exchange";			// Joules/metres squared
    	unit[55]="erg/cm";	conversion[55]=1.0E-15; 			type[55]="exchange";			// erg/cm squared
	 	unit[56]="erg/cm2";	conversion[56]=1.0e-23; 			type[56]="mm_energy";		// erg/cm2 squared to J/A2

      // Velocity
      unit[57]="A/s";      conversion[57]=1.0;                 type[57]="velocity";        // Angstrom per second
      unit[58]="m/s";      conversion[58]=1.0E10;              type[58]="velocity";        // Metres per second
      unit[59]="cm/s";     conversion[59]=1.0E8;               type[59]="velocity";        // Centimetres per second
      unit[60]="mm/s";     conversion[60]=1.0E7;               type[60]="velocity";        // Millimetres per second
      unit[61]="um/s";     conversion[61]=1.0E4;               type[61]="velocity";        // Micrometres per second
      unit[62]="nm/s";     conversion[62]=1.0E1;               type[62]="velocity";        // Nanometres per second
      unit[63]="pm/s";     conversion[63]=1.0E-3;              type[63]="velocity";        // Picometres per second
      unit[64]="in/s";     conversion[64]=2.54E8;              type[64]="velocity";        // Imperial Inches per second
      unit[65]="feet/s";   conversion[65]=3.12E9;              type[65]="velocity";        // Imperial Feet per second
      unit[66]="A/ms";     conversion[66]=1.0E3;               type[66]="velocity";        // Angstrom per milliseconds
      unit[67]="A/us";     conversion[67]=1.0E6;               type[67]="velocity";        // Angstrom per microseconds
      unit[68]="A/ns";     conversion[68]=1.0E9;               type[68]="velocity";        // Angstrom per nanoseconds
      unit[69]="A/ps";     conversion[69]=1.0E12;              type[69]="velocity";        // Angstrom per picoseconds
      unit[70]="nm/ns";		conversion[70]=1.0E8;               type[70]="velocity";        // Nanometers per nanosecond

      //

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
		terminaltextcolor(RED);
		std::cerr << "Error during unit conversion - unit \'"<< input_unit << "\' not found" << std::endl;
		terminaltextcolor(WHITE);
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
		terminaltextcolor(RED);
		std::cerr << "Error during unit reversion - unit \'"<< output_unit << "\' not found" << std::endl;
		terminaltextcolor(WHITE);
		err::vexit();

		return EXIT_SUCCESS;

	}

} // End of namespace
