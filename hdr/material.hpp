#include <iostream>
#include <string>
#include <vector>

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
		double Ku1_SI;
		double Ku2_SI;
		std::vector<double> KuVec_SI;
		double Ku;
		std::vector<double> KuVec;
		double Kc1_SI;
		double Kc2_SI;
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


	//extern materials_t material;
	//Integration parameters
	//extern double alpha;
	extern double dt_SI;
	extern double dt;
	extern double half_dt;
	extern double gamma_SI;
	//extern double Jij_SI;
	//extern double mu_s_SI;
	//extern double Ku_SI;
	
	//extern double one_oneplusalpha_sq;
	//extern double alpha_oneplusalpha_sq;
	//extern double Jij;
	//extern double Ku;
	//extern double H_th_sigma;

	
	// Functions
	extern int initialise(std::string);
	extern int print_mat();
	extern int default_system();	
	extern int single_spin_system();
	extern int set_derived_parameters();
	

}

// Alias deprecated material_parameters to mp namespace
namespace material_parameters=mp;
