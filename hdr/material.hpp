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
		double Ku1_SI; // SI uniaxial anisotropy constant
		double Ku2_SI;
		std::vector<double> KuVec_SI; // SI anisotropy tensor
		double Ku; // normalised uniaxial anisotropy constant
		std::vector<double> KuVec; // normalised anisotropy tensor
		std::vector<double> UniaxialAnisotropyUnitVector; // unit vector for material uniaxial anisotropy
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
		bool constrained; // specifies primary or alternate integrator
		
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
	extern std::vector <double> MaterialCubicAnisotropyArray;

	
	// Functions
	extern int initialise(std::string);
	extern int print_mat();
	extern int default_system();	
	extern int single_spin_system();
	extern int set_derived_parameters();
	

}

// Alias deprecated material_parameters to mp namespace
namespace material_parameters=mp;
