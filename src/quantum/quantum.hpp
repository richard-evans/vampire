//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2025. All rights reserved.
//
//   Email: fried-conrad.weber@uni-potsdam.de
//
//------------------------------------------------------------------------------
//

#ifndef QUANTUM_HPP
#define QUANTUM_HPP

namespace quantum{

   //---------------------------------------------------------------------------
   // Function prototypes
   //---------------------------------------------------------------------------
   void initialize();
   double get_field(int atom, int component, double time_offset);
   void increment_time();
   void llg_step();

   // Function to process input file parameters for quantum module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

   // Material parameter accessors
   double get_A(int material_index);
   double get_Gamma(int material_index);
   double get_omega0(int material_index);

}

#endif
