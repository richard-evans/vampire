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
//  Declarations and matching functions for variables of type dimension
//
//  (c) R F L Evans 2013
//
//-----------------------------------------------------------------------------
//
#include "./vio.ihpp"

namespace vio{

   vampire_iv_t<double> dimensions_a;

   
void initialise_dimensions(){

   //--------------------------------------------------------------
   // Unit cell size
   //--------------------------------------------------------------
   dimensions_a.name="a";
   dimensions_a.category="dimensions";
   dimensions_a.unit_type="length";
   dimensions_a.file_type="input";
   dimensions_a.value=3.54; // Angstroms
   dimensions_a.min_value=0.1;   // 0.1 Angstroms
   dimensions_a.max_value=1.0e7; // 1 millimetre 
   dimensions_a.range_error_message="0.1 Angstroms to 1 Millimetre";
   dimensions_a.list_error_alternatives="";
   dimensions_a.deprecated=false;
   dimensions_a.deprecated_alternative="";
   //|<------------------------------80 characters-------------------------------->|
   dimensions_a.documentation=
     "        Defines the unit cell lattice constant along the $x$ direction.\n"
     "        Default units are Angstroms.";
   dimensions_a.is_initialised=true;
   return;
  
}
   
//int match_dimension(std::string const, std::string const, std::string const, int const);

int match_dimension(std::string const category, std::string const keyword, std::string const value, std::string const unit, int const line){

      if(dimensions_a.test(keyword)){
         cs::unit_cell_size[0]=dimensions_a.check_value(value, keyword, line, unit);
         return EXIT_SUCCESS;
      }
      
         
/*         
         double a=atof(value.c_str());
         string unit_type;
         units::convert(unit,a,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::unit_cell_size[0]=a;
            cs::unit_cell_size[1]=a;
            cs::unit_cell_size[2]=a;
            //std::cout << "ax: " << cs::unit_cell_size[0] << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }

   std::string test="a";
      if(word==test){
         double a=atof(value.c_str());
         string unit_type;
         units::convert(unit,a,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::unit_cell_size[0]=a;
            cs::unit_cell_size[1]=a;
            cs::unit_cell_size[2]=a;
            //std::cout << "ax: " << cs::unit_cell_size[0] << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }
      
            //-------------------------------------------------------------------
      // System dimension variables
      //-------------------------------------------------------------------
      std::string prefix="dimension:";
      std::string none = "none";
      
      std::string test="a";
      if(word==test){
         double a=atof(value.c_str());
         string unit_type;
         units::convert(unit,a,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::unit_cell_size[0]=a;
            cs::unit_cell_size[1]=a;
            cs::unit_cell_size[2]=a;
            //std::cout << "ax: " << cs::unit_cell_size[0] << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }
      else
      test="c";
      if(word==test){
         double c=atof(value.c_str());
         string unit_type;
         units::convert(unit,c,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::unit_cell_size[2]=c;
            //std::cout << "ax: " << cs::unit_cell_size[0] << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }
      else
      //--------------------------------------------------------------------
      test="ax";
      if(word==test){
         double ax=atof(value.c_str());
         string unit_type;
         units::convert(unit,ax,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::unit_cell_size[0]=ax;
            //std::cout << "ax: " << cs::unit_cell_size[0] << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }
      else
      //--------------------------------------------------------------------
      test="ay";
      if(word==test){
         double ay=atof(value.c_str());
         string unit_type;
         units::convert(unit,ay,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::unit_cell_size[1]=ay;
            //std::cout << "ax: " << cs::unit_cell_size[0] << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }
      else
      //--------------------------------------------------------------------
      test="az";
      if(word==test){
         double az=atof(value.c_str());
         string unit_type;
         units::convert(unit,az,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::unit_cell_size[2]=az;
            //std::cout << "az: " << cs::unit_cell_size[0] << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }
      else
      //--------------------------------------------------------------------
      test="d";
      if(word==test){
         double d=atof(value.c_str());
         string unit_type;
         units::convert(unit,d,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::system_dimensions[0]=d;
            cs::system_dimensions[1]=d;
            cs::system_dimensions[2]=d;
            //std::cout << "ax: " << cs::unit_cell_size[0] << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }
      else
      //--------------------------------------------------------------------
      test="dx";
      if(word==test){
         double dx=atof(value.c_str());
         string unit_type;
         units::convert(unit,dx,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::system_dimensions[0]=dx;
            //std::cout << "ax: " << cs::unit_cell_size[0] << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }
      else
      //--------------------------------------------------------------------
      test="dy";
      if(word==test){
         double dy=atof(value.c_str());
         string unit_type;
         units::convert(unit,dy,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::system_dimensions[1]=dy;
            //std::cout << "ax: " << cs::unit_cell_size[0] << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }
      else
      //--------------------------------------------------------------------
      test="dz";
      if(word==test){
         double dz=atof(value.c_str());
         string unit_type;
         units::convert(unit,dz,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::system_dimensions[2]=dz;
            //std::cout << "ax: " << cs::unit_cell_size[0] << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'"<< std::endl;
            err::vexit();
         }
      }
      else
      //--------------------------------------------------------------------
      test="particle-size";
      if(word==test){
         double psize=atof(value.c_str());
         string unit_type;
         units::convert(unit,psize,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::particle_scale=psize;
            //std::cout << "particle_size: " << mp::particle_scale << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }
      else
      //--------------------------------------------------------------------
      test="particle-spacing";
      if(word==test){
         double pspacing=atof(value.c_str());
         string unit_type;
         units::convert(unit,pspacing,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cs::particle_spacing=pspacing;
            //std::cout << "particle_spacing: " << mp::particle_spacing << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }
      else
      //--------------------------------------------------------------------
      test="particle-array-offset-x";
      if(word==test){
         double paox=atof(value.c_str());

         // Test for valid range
         check_for_valid_value(paox, word, line, prefix, unit, "length", 0.0, 1.0e4,"input","0.0 - 1.0 micrometre");

         // If within valid range assign value and return
         cs::particle_array_offset_x=paox;

         return EXIT_SUCCESS;
      }
      else
      //--------------------------------------------------------------------
      test="particle-array-offset-y";
      if(word==test){
         double paoy=atof(value.c_str());

         // Test for valid range
         check_for_valid_value(paoy, word, line, prefix, unit, "length", 0.0, 1.0e4,"input","0.0 - 1.0 micrometre");

         // If within valid range assign value and return
         cs::particle_array_offset_y=paoy;

         return EXIT_SUCCESS;
      }
      else
      //--------------------------------------------------------------------
      test="cell-size";
      if(word==test){
         double cs=atof(value.c_str());
         string unit_type;
         units::convert(unit,cs,unit_type);
         string str="length";
         if(unit_type==str || unit_type==none){
            cells::size=cs;
            //std::cout << "cell size: " << cells::size << std::endl;
            return EXIT_SUCCESS;
         }
         else{
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'dimension:" << word << "\'" << std::endl;
            err::vexit();
         }
      }*/
      //--------------------------------------------------------------------
      else{
         std::cerr << "Error - Unknown control statement "<< category << ":"<< keyword << " on line " << line << " of input file. Exiting." << std::endl;
         return EXIT_FAILURE;
      }
   
}
  

} // end of namespace vio

