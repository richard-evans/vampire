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

//--------------------------------------------------------------------
// Class to contain parameters for lattice anisotropy calculation.
//
// Tabulated values are read from a file and added point-wise to
// the class. During variable initialisation interpolating functions
// are determined to calculate k(T)
//
class lattice_anis_t{

   private:

      unsigned int Tmax; // maximum array value in tabulated function
      double k_Tmax; // value of anisotropy at k_Tmax (used for all T > Tmax)

      std::vector<unsigned int> T; // input temperature points
      std::vector<double> k; // input lattice anisotropy points
      std::vector<double> m; // calculated m value
      std::vector<double> c; // calculated c value

   public:

      void add_point(double, double);
      void set_interpolation_table();
      double get_lattice_anisotropy_constant(double);
      void output_interpolated_function(int);

};

namespace mp
{
using std::string;

	//----------------------------------
	// Material Container
	//----------------------------------

	const int max_materials=100;
	extern int num_materials;

	class materials_t {
		public:
		// input parameters
		string name;
		string element;
		
		double alpha;
		double mu_s_SI;
		double magnetisation;
		double Ku1_SI; /// SI uniaxial anisotropy constant
		double Ku2_SI;
      double Ku3_SI;
      double Klatt_SI;
		std::vector<double> KuVec_SI; /// SI anisotropy tensor
		double Ku; /// normalised uniaxial anisotropy constant
      double Ku2; /// normalised uniaxial anisotropy constant
      double Ku3; /// normalised uniaxial anisotropy constant
      double Klatt; /// normalised lattice anisotropy
      std::vector<double> KuVec; /// normalised anisotropy tensor
		std::vector<double> UniaxialAnisotropyUnitVector; /// unit vector for material uniaxial anisotropy
		double Kc1_SI;
		double Kc2_SI;
		double Kc;
		double Ks_SI;
		double Ks;
		
		double gamma_rel;
		double Jij_matrix_SI[max_materials];

		double initial_spin[3];
		bool random_spins;
		
		double min;
		double max;
		int geometry; ///< 0 (geometry disabled, 1+ geometry enabled with 1+ points
		double geometry_coords[100][2];
		double core_shell_size;
		double interface_roughness;
		double density;
		double intermixing[max_materials];
		double cutoff;

		bool alloy_master;
		int alloy_class;
		double alloy[max_materials];
				
		bool continuous;	///< Specifies if a material is continuous (overrides granularity in the layer)
		bool moment_flag;	///< Specifies whether moment is set explicitly or from magnetisation
		bool anis_flag;	///< Specifies whether anisotropy is set explicitly or as energy density
		
		double one_oneplusalpha_sq;
		double alpha_oneplusalpha_sq;
		double Jij_matrix[max_materials];
		double H_th_sigma;
		bool constrained; /// specifies primary or alternate integrator

		double temperature; /// Kelvin
		bool couple_to_phonon_temperature; ///true/false
		double applied_field_strength; /// Tesla
		std::vector<double> applied_field_unit_vector; /// unit vector for material uniaxial anisotropy
		double fmr_field_strength; // Tesla
		double fmr_field_frequency; // Hz
		std::vector<double> fmr_field_unit_vector; /// unit vector for material uniaxial anisotropy
		bool fill; /// flag to determine of material fills voided space
     lattice_anis_t lattice_anisotropy; // class containing lattice anisotropy data
		
		materials_t();
		int print();
	};



	extern std::vector <materials_t> material;

	extern double dt_SI;
	extern double dt;
	extern double half_dt;
	extern double gamma_SI;

	// Unrolled material parameters for speed
	extern std::vector <double> MaterialMuSSIArray;
	extern std::vector <zkval_t> MaterialScalarAnisotropyArray;
	extern std::vector <zkten_t> MaterialTensorAnisotropyArray;
   extern std::vector <double> material_second_order_anisotropy_constant_array;
   extern std::vector <double> material_sixth_order_anisotropy_constant_array;
	extern std::vector <double> MaterialCubicAnisotropyArray;

	
	// Functions
	extern int initialise(std::string);
	extern int print_mat();
	extern int default_system();	
	extern int single_spin_system();
	extern int set_derived_parameters();
	

}

/// Alias deprecated material_parameters to mp namespace
namespace material_parameters=mp;
