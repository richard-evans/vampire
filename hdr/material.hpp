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

#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <iostream>
#include <string>
#include <vector>

class zkval_t{
   public:
   double K;

   // constructor
   zkval_t():
      K(0.0)
   {
   };
};

class zkten_t{
   public:
   double K[3][3];

   // constructor
   zkten_t()
   {
      K[0][0]=0.0;
      K[0][1]=0.0;
      K[0][2]=0.0;

      K[1][0]=0.0;
      K[1][1]=0.0;
      K[1][2]=0.0;

      K[2][0]=0.0;
      K[2][1]=0.0;
      K[2][2]=0.0;
   };
};

namespace mp
{

   //----------------------------------
   // Material Container
   //----------------------------------

   const int max_materials=100;
   extern int num_materials;

   class materials_t {
      public:
      // input parameters
      std::string name;
      std::string element;


      bool micromagnetic_enabled;
      std::vector<double> SAF;
      bool enable_SAF;
      std::vector < double > pinning_field_unit_vector;
      double pinning_field_strength;

      std::vector< double> EF_MM;
      std::vector<bool> override_atomsitic;

      double alpha;
                double alpha_eq;
      double mu_s_SI;
       double mass;
      double magnetisation;

      double gamma_rel;

      double core_shell_size;
      double interface_roughness;
      double density;
      double intermixing[max_materials];
      double cutoff;

      bool alloy_master;
      int alloy_class;
      double alloy[max_materials];

      bool continuous;  ///< Specifies if a material is continuous (overrides granularity in the layer)
      bool moment_flag; ///< Specifies whether moment is set explicitly or from magnetisation

      double one_oneplusalpha_sq;
      double alpha_oneplusalpha_sq;
      double H_th_sigma;
//modified SLD M Strungaru
                double H_th_sigma_eq;

      bool constrained; /// specifies primary or alternate integrator

      double temperature; /// Kelvin
      double maximum_temperature; /// Kelvin
      double minimum_temperature; /// Kelvin
      bool couple_to_phonon_temperature; ///true/false
      double applied_field_strength; /// Tesla
      std::vector<double> applied_field_unit_vector; /// unit vector for material uniaxial anisotropy
      double fmr_field_strength; // Tesla
      double fmr_field_frequency; // Hz
      std::vector<double> fmr_field_unit_vector; /// unit vector for material uniaxial anisotropy
      bool fill; /// flag to determine if material fills voided space
      double temperature_rescaling_alpha; // temperature rescaling exponent
      double temperature_rescaling_Tc; // temperaure rescaling Tc
      int non_magnetic;

      materials_t();
      int print();
   };



   extern std::vector <materials_t> material;

   // temporary array of materials for reading in material data
   extern std::vector <materials_t> read_material;

   extern double dt_SI;
   extern double dt;
   extern double half_dt;
   extern double gamma_SI;

   // Unrolled material parameters for speed
   extern std::vector <double> mu_s_array;
   extern std::vector <zkval_t> MaterialScalarAnisotropyArray;
   extern std::vector <zkten_t> MaterialTensorAnisotropyArray;
   extern std::vector <double> material_second_order_anisotropy_constant_array;
   extern std::vector <double> material_sixth_order_anisotropy_constant_array;
   extern std::vector <double> material_spherical_harmonic_constants_array;
   extern std::vector <double> MaterialCubicAnisotropyArray;

   // Functions
   extern int initialise(std::string);
   extern int print_mat();
   extern int default_system();
   extern int single_spin_system();
   extern int set_derived_parameters();

   // input and material file parameter matching functions
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

}

#endif // MATERIAL_H_
