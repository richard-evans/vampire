#include <iostream>
#include <string>
#include <valarray>

namespace mp
{
using std::string;
using std::valarray;
	//----------------------------------
	// Material Container
	//----------------------------------

	const int max_materials=100;
	extern int num_materials;

	class materials_t {
		public:
		// input parameters
		string name;
		string hamiltonian_type;
		string element;
		string crystal_structure;
		
		double alpha;
		double mu_s_SI;
		double magnetisation;
		double Ku1_SI;
		double Ku2_SI;
		double Kc1_SI;
		double Kc2_SI;
		double Ks_SI;
		
		double gamma_rel;
		double Jij_matrix_SI[max_materials];

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
		//int ianis_flag;
		// derived parameters
		int num_nearest_neighbours;
		int hamiltonian_num_neighbours;
		
		double one_oneplusalpha_sq;
		double alpha_oneplusalpha_sq;
		double Jij_matrix[max_materials];
		double Ku;
		double H_th_sigma;
		
		materials_t();
		int print();
	};



	extern valarray <materials_t> material;


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
	//----------------------------------
	//Input System Parameters
	//----------------------------------
	extern int particle_creation_parity;
	extern double system_dimensions[3];
	extern double particle_scale;
	extern double particle_spacing;
	
	
	//System Parameters
	extern double lattice_constant[3];
	extern double lattice_space_conversion[3];
	extern string crystal_structure;
	extern bool single_spin;
	extern string hamiltonian_type;
	//extern string atomic_element[4];
	
	
	
	//----------------------------------
	//Derived System Parameters
	//----------------------------------
	//extern int num_nearest_neighbours;
	extern int hamiltonian_num_neighbours;
	extern int int_system_dimensions[3];
	
	//----------------------------------
	// System creation flags
	//----------------------------------
	
	extern int system_creation_flags[10];
	
	// Functions
	extern int initialise(std::string);
	extern int print_mat();
	extern int default_system();	
	extern int single_spin_system();
	extern int set_derived_parameters();
	

}

// Alias deprecated material_parameters to mp namespace
namespace material_parameters=mp;
