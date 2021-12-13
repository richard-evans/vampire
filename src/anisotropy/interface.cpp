//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: sw766@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>
#include <sstream>

// Vampire headers
#include "anisotropy.hpp"
#include "errors.hpp"
#include "units.hpp"
#include "vio.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for anisotropy module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="anisotropy";
      if(key!=prefix) return false;

      //-------------------------------------------------------------------
      std::string test="surface-anisotropy-threshold";
      if(word==test){
          // test for native keyword
          test="native";
          if(value==test){
             internal::native_neel_anisotropy_threshold = true;
             return EXIT_SUCCESS;
          }
          int sat=atoi(value.c_str());
          // Test for valid range
          vin::check_for_valid_int(sat, word, line, prefix, 0, 1000000000,"input","0 - 1,000,000,000");
          internal::neel_anisotropy_threshold = sat;
          return true;
      }
      //-------------------------------------------------------------------
      test="surface-anisotropy-nearest-neighbour-range";
      if(word==test){
          // Test for valid range
          double r = atof(value.c_str());
          vin::check_for_valid_value(r, word, line, prefix, unit, "length", 0.0, 1.0e9,"input","0.0 - 1,000,000,000");
          internal::nearest_neighbour_distance = r;
          return true;
      }
      //-------------------------------------------------------------------
      test="enable-bulk-neel-anisotropy";
      if(word==test){
          // Enable large threshold to force calculation of Neel anisotropy for all atoms
          internal::neel_anisotropy_threshold = 1000000000;
          return true;
      }
      //-------------------------------------------------------------------
      test="neel-anisotropy-exponential-factor";
      if(word==test){
          // Enable range dependent Neel anisotropy Lij(r) = exp(-F(r-r0)/r0)
          // F should be steepness of function indicating rate of decay with r
          double F = atof(value.c_str());
          vin::check_for_valid_value(F, word, line, prefix, unit, "none", 0.01, 100.0,"input","0.01 - 100");
          internal::neel_exponential_factor = F;
          internal::neel_range_dependent = true;
          return true;
      }
      //-------------------------------------------------------------------
      test="neel-anisotropy-exponential-range";
      if(word==test){
          // Enable range dependent Neel anisotropy Lij(r) = exp(-F(r-r0)/r0)
          // r should be approximately nearest neighbour range ~ 2.5 angstroms
          double r = atof(value.c_str());
          vin::check_for_valid_value(r, word, line, prefix, unit, "length", 0.0001, 1000.0,"input","0.0001 - 1,000");
          internal::neel_exponential_range = r;
          internal::neel_range_dependent = true;
          return true;
      }
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index, const int max_materials){

      // add prefix string
      std::string prefix="material:";

      // Check for empty material parameter array and resize to avoid segmentation fault
      if(internal::mp.size() == 0){
         internal::mp.resize(max_materials);
      }

      //------------------------------------------------------------
      // Check for material properties
      //------------------------------------------------------------
      //Minimal orthogonality TBD
      std::string test  = "second-order-uniaxial-anisotropy-constant"; // new form (preferred)
      std::string test2 = "uniaxial-anisotropy-constant"; // legacy form (deprecated but probably never obsoleted)
      if( (word == test) || (word == test2) ){
         double ku2 = atof(value.c_str());
         vin::check_for_valid_value(ku2, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e-17 J/atom");
         internal::mp[super_index].ku2 = ku2;
         internal::enable_uniaxial_second_order = true; // Switch on second order tensor calculation for all spins
         return true;
      }

    //  Triaxial anisotropy in second and fourth order
      test="second-order-triaxial-anisotropy-vector";
      if(word == test){
         std::vector<double> u(3);
         // read values from string
         u = vin::doubles_from_string(value);
         //std::cout << u[0] << '\t' << u[1] << "\t" << u[2] << std::endl;
         // check for sane input and normalise if necessary
         vin::check_for_valid_vector(u, word, line, prefix, unit, "anisotropy", -1e-10, 1e-10,"material"," < +/- 1.0e-10");
         // Copy sanitised unit vector to material
         internal::ku_triaxial_vector_x[super_index] = u[0];
         internal::ku_triaxial_vector_y[super_index] = u[1];
         internal::ku_triaxial_vector_z[super_index] = u[2];
         internal::enable_triaxial_anisotropy = true;
         return true;

      }
      test="fourth-order-triaxial-anisotropy-vector";
      if(word == test){
         std::vector<double> u(3);
         // read values from string
         u = vin::doubles_from_string(value);
         //std::cout << u[0] << '\t' << u[1] << "\t" << u[2] << std::endl;
         // check for sane input and normalise if necessary
         vin::check_for_valid_vector(u, word, line, prefix, unit, "anisotropy", -1e-10, 1e-10,"material"," < +/- 1.0e-10");
         // Copy sanitised unit vector to material
         internal::ku4_triaxial_vector_x[super_index] = u[0];
         internal::ku4_triaxial_vector_y[super_index] = u[1];
         internal::ku4_triaxial_vector_z[super_index] = u[2];
         internal::enable_triaxial_fourth_order = true;
         return true;

      }

      test="second-order-triaxial-basis-vector-1";
      if(word == test){
         std::vector<double> u(3);
         // read values from string
         u = vin::doubles_from_string(value);
         //std::cout << u[0] << '\t' << u[1] << "\t" << u[2] << std::endl;
         // check for sane input and normalise if necessary
         vin::check_for_valid_unit_vector(u, word, line, prefix, "material");
         // Copy sanitised unit vector to material
         internal::ku_triaxial_basis1x[super_index] = u[0];
         internal::ku_triaxial_basis1y[super_index] = u[1];
         internal::ku_triaxial_basis1z[super_index] = u[2];
         internal::triaxial_second_order_fixed_basis[super_index] = false;
      //   std::cout << u[0] << '\t' << u[1] << '\t' << u[2] << "\t" << super_index <<  std::endl;
      //   std::cin.get();
         return true;

      }

      test="second-order-triaxial-basis-vector-2";
      if(word == test){
         std::vector<double> u(3);
         // read values from string
         u = vin::doubles_from_string(value);
         //std::cout << u[0] << '\t' << u[1] << "\t" << u[2] << std::endl;
         // check for sane input and normalise if necessary
        vin::check_for_valid_unit_vector(u, word, line, prefix, "material");
         // Copy sanitised unit vector to material
         internal::ku_triaxial_basis2x[super_index] = u[0];
         internal::ku_triaxial_basis2y[super_index] = u[1];
         internal::ku_triaxial_basis2z[super_index] = u[2];
         internal::triaxial_second_order_fixed_basis[super_index] = false;
         return true;

      }

      test="second-order-triaxial-basis-vector-3";
      if(word == test){
         std::vector<double> u(3);
         // read values from string
         u = vin::doubles_from_string(value);
         //std::cout << u[0] << '\t' << u[1] << "\t" << u[2] << std::endl;
         // check for sane input and normalise if necessary
         vin::check_for_valid_unit_vector(u, word, line, prefix, "material");
         // Copy sanitised unit vector to material
         internal::ku_triaxial_basis3x[super_index] = u[0];
         internal::ku_triaxial_basis3y[super_index] = u[1];
         internal::ku_triaxial_basis3z[super_index] = u[2];
         internal::triaxial_second_order_fixed_basis[super_index] = false;
         return true;

      }


      test="fourth-order-triaxial-basis-vector-1";
      if(word == test){
         std::vector<double> u(3);
         // read values from string
         u = vin::doubles_from_string(value);
         //std::cout << u[0] << '\t' << u[1] << "\t" << u[2] << std::endl;
         // check for sane input and normalise if necessary
         vin::check_for_valid_unit_vector(u, word, line, prefix, "material");
         // Copy sanitised unit vector to material
         internal::ku4_triaxial_basis1x[super_index] = u[0];
         internal::ku4_triaxial_basis1y[super_index] = u[1];
         internal::ku4_triaxial_basis1z[super_index] = u[2];
         internal::triaxial_fourth_order_fixed_basis[super_index] = false;
      //   std::cout << u[0] << '\t' << u[1] << '\t' << u[2] << "\t" << super_index <<  std::endl;
      //   std::cin.get();
         return true;

      }

      test="fourth-order-triaxial-basis-vector-2";
      if(word == test){
         std::vector<double> u(3);
         // read values from string
         u = vin::doubles_from_string(value);
         //std::cout << u[0] << '\t' << u[1] << "\t" << u[2] << std::endl;
         // check for sane input and normalise if necessary
         vin::check_for_valid_unit_vector(u, word, line, prefix, "material");
         // Copy sanitised unit vector to material
         internal::ku4_triaxial_basis2x[super_index] = u[0];
         internal::ku4_triaxial_basis2y[super_index] = u[1];
         internal::ku4_triaxial_basis2z[super_index] = u[2];
         internal::triaxial_fourth_order_fixed_basis[super_index] = false;
         return true;

      }

      test="fourth-order-triaxial-basis-vector-3";
      if(word == test){
         std::vector<double> u(3);
         // read values from string
         u = vin::doubles_from_string(value);
         //std::cout << u[0] << '\t' << u[1] << "\t" << u[2] << std::endl;
         // check for sane input and normalise if necessary
         vin::check_for_valid_unit_vector(u, word, line, prefix, "material");
         // Copy sanitised unit vector to material
         internal::ku4_triaxial_basis3x[super_index] = u[0];
         internal::ku4_triaxial_basis3y[super_index] = u[1];
         internal::ku4_triaxial_basis3z[super_index] = u[2];
         internal::triaxial_fourth_order_fixed_basis[super_index] = false;
         return true;

      }

      //------------------------------------------------------------
      //Minimal orthogonality TBD
      test = "fourth-order-uniaxial-anisotropy-constant";
      if( word == test ){
         double ku4 = atof(value.c_str());
         vin::check_for_valid_value(ku4, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e-17 J/atom");
         internal::mp[super_index].ku4 = ku4;
         internal::enable_uniaxial_fourth_order = true; // Switch on second order tensor calculation for all spins (from spherical harmonics)
         return true;
      }
      //------------------------------------------------------------
      //Implementation of biaxial fourth-order anisotropy (simple version)
      test = "fourth-order-biaxial-anisotropy-constant";
      if( word == test ){
         double ku4 = atof(value.c_str());
         vin::check_for_valid_value(ku4, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e-17 J/atom");
         internal::mp[super_index].ku4 = ku4;
         internal::enable_biaxial_fourth_order_simple = true;
         return true;
      }
      //------------------------------------------------------------      
      //Minimal orthogonality
      test = "sixth-order-uniaxial-anisotropy-constant";
      if( word == test ){
         double ku6 = atof(value.c_str());
         vin::check_for_valid_value(ku6, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e-17 J/atom");
         internal::mp[super_index].ku6 = ku6;
         internal::enable_uniaxial_sixth_order = true; // Switch on second order tensor calculation for all spins (from spherical harmonics)
         return true;
      }
      //------------------------------------------------------------
      //Minimal orthogonality
      test = "fourth-order-cubic-anisotropy-constant"; // new form (preferred)
      test2  = "cubic-anisotropy-constant"; // legacy form (deprecated but probably never obsoleted)
      if( (word == test) || (word == test2) ){
         double kc4 = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(kc4, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e-17 J/atom");
         internal::mp[super_index].kc4 = kc4;
         if(internal::enable_cubic_fourth_order_rotation == false) internal::enable_cubic_fourth_order = true; // Switch on second order tensor calculation for all spins (from spherical harmonics)
         return true;
      }
      //------------------------------------------------------------
      //Minimal orthogonality
      test = "sixth-order-cubic-anisotropy-constant";
      if( word == test ){
         double kc6 = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(kc6, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e-17 J/atom");
         internal::mp[super_index].kc6 = kc6;
         internal::enable_cubic_sixth_order = true; // Switch on second order tensor calculation for all spins (from spherical harmonics)
         return true;
      }
      //------------------------------------------------------------
      test="neel-anisotropy-constant"; // new and preferred form
      test2="surface-anisotropy-constant"; // deprecated but never obsolete form
      if( (word == test) || (word == test2) ){
         double kij = atof(value.c_str());
         vin::check_for_valid_value(kij, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e17");
         internal::mp[super_index].kij[sub_index] = kij;
         internal::enable_neel_anisotropy = true;
         return true;
      }
      //------------------------------------------------------------
      test = "lattice-anisotropy-constant";
      if( word == test ){
         double kl = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(kl, word, line, prefix, unit, "energy", -1.0e-17, 1.0e17,"material","-1e17 - 1e17 J/atom");
         internal::mp[super_index].k_lattice = kl;
         internal::enable_lattice_anisotropy = true;
         return true;
      }
      //------------------------------------------------------------
      test="lattice-anisotropy-file";
      if(word==test){

         // Open lattice file
         std::stringstream latt_file;
         latt_file.str( vin::get_string(value.c_str(), "material", line) );

         // specify number of points to be read
         int num_pts=0;

         // Read in number of temperature points
         latt_file >> num_pts;

         // Check for valid number of points
         if(num_pts<=1){
            std::cerr << "Error in lattice-anisotropy-file " << value.c_str() << " on line " << line << " of material file. The first number must be an integer greater than 1. Exiting." << std::endl;
            zlog << zTs() << "Error in lattice-anisotropy-file " << value.c_str() << " on line " << line << " of material file. The first number must be an integer greater than 1. Exiting." << std::endl;
            return false;
         }

         // Loop over all lines in file
         for(int c=0;c<num_pts;c++){

            // temporary variables
            double T;
            double k;

            // Read in points and add them to material
            latt_file >> T >> k;

            // Check for premature end of file
            if(latt_file.eof()){
               std::cerr << "Error in lattice anisotropy-file " << value.c_str() << " on line " << line << " of material file. End of file reached before reading all values. Exiting" << std::endl;
               zlog << zTs() << "Error in lattice anisotropy-file " << value.c_str() << " on line " << line << " of material file. End of file reached before reading all values. Exiting" << std::endl;
               return false;
            }
            internal::mp[super_index].lattice_anisotropy.add_point(T,k);
         }

         return true;

      }
      //------------------------------------------------------------
      test = "uniaxial-anisotropy-direction";
      if(word == test){
         // set up test comparisons
         test="random";
         test2="random-grain";
         // test for random anisotropy directions
         if( value == test ){
            internal::mp[super_index].random_anisotropy = true;
            internal::mp[super_index].random_grain_anisotropy = false;
            internal::enable_random_anisotropy = true; // enable initialisation of random anisotropy
         }
         // test for random grain anisotropy
         else if( value == test2 ){
            internal::mp[super_index].random_anisotropy = false;
            internal::mp[super_index].random_grain_anisotropy = true;
            internal::enable_random_anisotropy = true;
         }
         // set easy axis unit vector for material
         else{
            // temporary storage container
            std::vector<double> u(3);
            // read values from string
            u = vin::doubles_from_string(value);
            // check for sane input and normalise if necessary
            vin::check_for_valid_unit_vector(u, word, line, prefix, "material");
            // Copy sanitised unit vector to material
            internal::mp[super_index].ku_vector = u;
         }
         return true;
      }
      //--------------------------------------
      // Direction 1
      //--------------------------------------
      test = "cubic-anisotropy-direction-1";
      if(word == test){
         // temporary storage container
         std::vector<double> u(3);
         // read values from string
         u = vin::doubles_from_string(value);
         // check for sane input and normalise if necessary
         vin::check_for_valid_unit_vector(u, word, line, prefix, "material");
         // Copy sanitised unit vector to material
         internal::mp[super_index].kc_vector1 = u;
         // enable rotated anisotropy and disable normal anisotropy
         internal::enable_cubic_fourth_order_rotation = true;
         internal::enable_cubic_fourth_order = false;
         return true;
      }
      //--------------------------------------
      // Direction 2
      //--------------------------------------
      test = "cubic-anisotropy-direction-2";
      if(word == test){
         // temporary storage container
         std::vector<double> u(3);
         // read values from string
         u = vin::doubles_from_string(value);
         // check for sane input and normalise if necessary
         vin::check_for_valid_unit_vector(u, word, line, prefix, "material");
         // Copy sanitised unit vector to material
         internal::mp[super_index].kc_vector2 = u;
         // enable rotated anisotropy and disable normal anisotropy
         internal::enable_cubic_fourth_order_rotation = true;
         internal::enable_cubic_fourth_order = false;
         return true;
      }
      //------------------------------------------------------------
      // Fourth order rotational anisotropy
      //------------------------------------------------------------
      test = "fourth-order-rotational-anisotropy-constant";
      if(word == test){
         double k4r = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(k4r, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e-17 J/atom");
         internal::mp[super_index].k4r = k4r;
         internal::enable_fourth_order_rotational = true;
         return true;
      }
      //------------------------------------------------------------
      /*test="uniaxial-anisotropy-tensor";
      if(word==test){
         std::vector<double> t;
         // read values from string
         t = vin::doubles_from_string(value);
         // check size is nine elements
         if( t.size() != 9 ){
            terminaltextcolor(RED);
            std::cerr << "Error in input file - material[" << super_index << "]:uniaxial-anisotropy-tensor must have nine values." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error in input file - material[" << super_index << "]:uniaxial-anisotropy-tensor must have nine values." << std::endl;
            return false;
         }
         string unit_type = "energy";
         // if no unit given, assume internal
         if( unit.size() != 0 ){
            units::convert(unit, t, unit_type);
         }
         string str = "energy";
         if( unit_type == str ){
            // Copy anisotropy vector to material
            internal::mp[super_index].ku_tensor = t;
            internal::enable_second_order_tensor = true; // Switch on second order tensor calculation for all spins
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'material:" << word << "\'"<< std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }*/
      //------------------------------------------------------------
      /*test="cubic-anisotropy-tensor";
      if(word==test){
         std::vector<double> t;
         // read values from string
         t = vin::doubles_from_string(value);
         // check size is nine elements
         if( t.size() != 9 ){
            terminaltextcolor(RED);
            std::cerr << "Error in input file - material[" << super_index << "]:uniaxial-anisotropy-tensor must have nine values." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error in input file - material[" << super_index << "]:uniaxial-anisotropy-tensor must have nine values." << std::endl;
            return false;
         }
         string unit_type = "energy";
         // if no unit given, assume internal
         if( unit.size() != 0 ){
            units::convert(unit, t, unit_type);
         }
         string str = "energy";
         if( unit_type == str ){
            // Copy anisotropy vector to material
            internal::mp[super_index].kc_tensor = t;
            internal::enable_fourth_order_tensor = true; // Switch on fourth order tensor calculation for all spins
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - unit type \'" << unit_type << "\' is invalid for parameter \'material:" << word << "\'"<< std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }*/
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of anisotropy namespace
