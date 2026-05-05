//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
// (c) R F L Evans 2026. All rights reserved.
//
//------------------------------------------------------------------------------

// Headers
#include "constants.hpp"
#include "errors.hpp"
#include "units.hpp"
#include "vio.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace units {

const double pi = M_PI;

// The unit table: a dynamic list of all recognised units.
// Using std::vector means there is no hard-coded size limit —
// just call add() in init() to register a new unit.
// This vector is an implementation detail; it is not exposed in units.hpp.
static std::vector<unit_t> unit_table;
static bool initialised = false;

// Register one unit with the table.
// Having all three values on one line makes it easy to scan the list
// and spot mistakes (wrong type, missing factor of 10, etc.).
static void add(const std::string& name, double conversion, const std::string& type){
   unit_table.push_back({ name, type, conversion });
}

//-----------------------------------------------------------------------------
// Populate the unit table.
// Each add() call registers one unit: (name, conversion factor, type).
// The conversion factor turns the user-supplied value into internal units.
// Internal units are: Angstrom (length), Joules (energy), Tesla (field),
//                     seconds (time), J/T (moment), J/T/A^3 (magnetisation).
//-----------------------------------------------------------------------------
int init(){

   if(err::check) std::cout << "units::init has been called" << std::endl;

   unit_table.clear();

   // --- Distances (internal unit: Angstrom) ---
   add("",       1.0,    "none");    // no unit specified — pass through unchanged
   add("A",      1.0,    "length");  // Angstrom (internal unit, factor = 1)
   add("m",      1.0E10, "length");  // Metres
   add("cm",     1.0E8,  "length");  // Centimetres
   add("mm",     1.0E7,  "length");  // Millimetres
   add("um",     1.0E4,  "length");  // Micrometres
   add("nm",     1.0E1,  "length");  // Nanometres
   add("pm",     1.0E-3, "length");  // Picometres
   add("in",     2.54E8, "length");  // Imperial inches
   add("feet",   3.12E9, "length");  // Imperial feet

   // --- Energies (internal unit: Joules) ---
   add("J",      1.0,                    "energy");  // Joules (internal unit)
   add("eV",     constants::eV,          "energy");  // Electron volts
   add("meV",    constants::eV * 1.0e-3, "energy");  // Millielectron volts
   add("erg",    1.0E-7,          "energy");  // Ergs
   add("Ryd",    2.17987208E-18,  "energy");  // Rydbergs
   add("mRyd",   2.17987208E-21,  "energy");  // Millirydbergs
   add("Ht",     4.35974417E-18,  "energy");  // Hartrees
   add("zJ",     1.0E-21,         "energy");  // Zeptojoules

   // --- Magnetic moment (internal unit: J/T) ---
   add("J/T",    1.0,             "moment");  // Joules/Tesla (internal unit)
   add("Amm",    1.0,             "moment");  // Amp metres squared
   add("erg/G",  1.0E-3,          "moment");  // Erg/Gauss
   add("abAcmcm",1.0E-3,          "moment");  // Abampere centimetre squared
   add("muB",    constants::muB,   "moment");  // Bohr magnetons
   add("eV/T",   constants::eV,   "moment");  // Electron volts/Tesla
   add("erg/Oe", 1.0E-3,          "moment");  // Erg/Oersted

   // --- Magnetisation (internal unit: J/T/Angstrom^3) ---
   // 1 A/m = 1e-30 J/T/A^3  (1 A/m * 1e-30 m^3/A^3 = 1e-30 J/T/A^3)
   add("A/m",    1.0E-30, "magnetisation");  // Amps/metre
   add("emu/cc", 1.0E-27, "magnetisation");  // emu/cubic cm
   add("J/T/AAA",1.0,     "magnetisation");  // Joules/Tesla/Angstrom^3 (internal unit)

   // --- Anisotropy (internal unit: J/atom) ---
   add("J/atom", 1.0,     "anisotropy");  // Joules/atom (internal unit)
   add("J/mmm",  1.0E-30, "anisotropy");  // Joules/metre^3
   add("erg/cc", 1.0E-31, "anisotropy");  // Erg/cm^3

   // --- Magnetic field (internal unit: Tesla) ---
   add("T",      1.0,     "field");  // Tesla (internal unit)
   add("mT",     1.0E-3,  "field");  // Millitesla
   add("uT",     1.0E-6,  "field");  // Microtesla
   add("Oe",     1.0E-4,  "field");  // Oersted
   add("kOe",    1.0E-1,  "field");  // Kilo-Oersted

   // --- Time (internal unit: seconds) ---
   add("s",      1.0,     "time");  // Seconds (internal unit)
   add("ms",     1.0E-3,  "time");  // Milliseconds
   add("us",     1.0E-6,  "time");  // Microseconds
   add("ns",     1.0E-9,  "time");  // Nanoseconds
   add("ps",     1.0E-12, "time");  // Picoseconds
   add("fs",     1.0E-15, "time");  // Femtoseconds
   add("as",     1.0E-18, "time");  // Attoseconds
   add("zs",     1.0E-21, "time");  // Zeptoseconds

   // --- Resistivity (internal unit: Ohm metres) ---
   add("Ohm-m",  1.0,     "resistivity");  // Ohm metres (internal unit)

   // --- Electric potential (internal unit: Volts) ---
   add("V",      1.0,     "potential");  // Volts (internal unit)
   add("mV",     1.0E-3,  "potential");  // Millivolts
   add("uV",     1.0E-6,  "potential");  // Microvolts
   add("nV",     1.0E-9,  "potential");  // Nanovolts

   // --- Frequency (internal unit: Hz) ---
   add("Hz",     1.0,     "frequency");  // Hertz (internal unit)
   add("kHz",    1.0E3,   "frequency");  // Kilohertz
   add("MHz",    1.0E6,   "frequency");  // Megahertz
   add("GHz",    1.0E9,   "frequency");  // Gigahertz
   add("THz",    1.0E12,  "frequency");  // Terahertz

   // --- Exchange stiffness (internal unit: J/Angstrom) ---
   add("J/m",    1.0E-10, "exchange");  // Joules/metre
   add("erg/cm", 1.0E-15, "exchange");  // Erg/centimetre

   // --- Micromagnetic energy density (internal unit: J/Angstrom^2) ---
   add("erg/cm2",1.0E-23, "mm_energy");  // Erg/centimetre^2

   // --- Velocity (internal unit: Angstrom/second) ---
   add("A/s",    1.0,     "velocity");  // Angstrom/second (internal unit)
   add("m/s",    1.0E10,  "velocity");  // Metres/second
   add("cm/s",   1.0E8,   "velocity");  // Centimetres/second
   add("mm/s",   1.0E7,   "velocity");  // Millimetres/second
   add("um/s",   1.0E4,   "velocity");  // Micrometres/second
   add("nm/s",   1.0E1,   "velocity");  // Nanometres/second
   add("pm/s",   1.0E-3,  "velocity");  // Picometres/second
   add("in/s",   2.54E8,  "velocity");  // Imperial inches/second
   add("feet/s", 3.12E9,  "velocity");  // Imperial feet/second
   add("A/ms",   1.0E3,   "velocity");  // Angstrom/millisecond
   add("A/us",   1.0E6,   "velocity");  // Angstrom/microsecond
   add("A/ns",   1.0E9,   "velocity");  // Angstrom/nanosecond
   add("A/ps",   1.0E12,  "velocity");  // Angstrom/picosecond
   add("nm/ns",  1.0E8,   "velocity");  // Nanometres/nanosecond

   // temperature C, F, K; angles degrees, rad, mrad;

   // Set initialised flag
   initialised = true;

   return EXIT_SUCCESS;

}

//-----------------------------------------------------------------------------
// Convert a value from the named unit into internal units.
// Returns EXIT_SUCCESS if the unit was found, EXIT_FAILURE otherwise.
// The 'type' string is set to the physical quantity (e.g. "length") so
// callers can check that the right kind of unit was supplied.
//-----------------------------------------------------------------------------
int convert(std::string input_unit, double& value, std::string& type){

   if(!initialised) init();

   for(const unit_t& u : unit_table){
      if(input_unit == u.name){
         value *= u.conversion;
         type   = u.type;
         return EXIT_SUCCESS;
      }
   }

   return EXIT_FAILURE;
}

//-----------------------------------------------------------------------------
// Convert a vector of values from the named unit into internal units.
// Calls vexit() on failure because every element needs a valid conversion.
//-----------------------------------------------------------------------------
void convert(std::string input_unit, std::vector<double>& value, std::string& type){

   if(!initialised) init();

   for(const unit_t& u : unit_table){
      if(input_unit == u.name){
         for(double& v : value) v *= u.conversion;
         type = u.type;
         return;
      }
   }

   terminaltextcolor(RED);
   std::cerr << "Error during unit conversion - unit '" << input_unit << "' not found" << std::endl;
   terminaltextcolor(WHITE);
   zlog << zTs() << "Error during unit conversion - unit '" << input_unit << "' not found" << std::endl;
   err::vexit();
}

//-----------------------------------------------------------------------------
// Convert a value from internal units back into the named external unit.
// This is the inverse of convert(): it divides rather than multiplies.
//-----------------------------------------------------------------------------
int revert(std::string output_unit, double& value, std::string& type){

   if(!initialised) init();

   for(const unit_t& u : unit_table){
      if(output_unit == u.name){
         value /= u.conversion;
         type   = u.type;
         return EXIT_SUCCESS;
      }
   }

   terminaltextcolor(RED);
   std::cerr << "Error during unit reversion - unit '" << output_unit << "' not found" << std::endl;
   terminaltextcolor(WHITE);
   err::vexit();

   return EXIT_SUCCESS;
}

} // end of namespace units
