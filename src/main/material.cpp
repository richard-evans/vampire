#include "material.hpp"

namespace mp{
	// Constructor
materials_t::materials_t ():
	name(""),
	hamiltonian_type("generic"),
	element("   "),
	crystal_structure(""),
	alpha(1.0),
	mu_s_SI(9.27400915e-24),
	magnetisation(0.0),
	Ku1_SI(0.0),
	Ku2_SI(0.0),
	Kc1_SI(0.0),
	Kc2_SI(0.0),
	Ks_SI(0.0),
	gamma_rel(1.0),
	random_spins(false),
	min(0.0),
	max(1.0),
	geometry(0),
	core_shell_size(0.0),
	interface_roughness(0.0),
	density(1.0),
	cutoff(0.0),
	alloy_master(false),
	alloy_class(0),
	continuous(false),
	moment_flag(true),
	anis_flag(true),
	num_nearest_neighbours(0),
	hamiltonian_num_neighbours(0),
	one_oneplusalpha_sq(0.5),
	alpha_oneplusalpha_sq(0.5),
	Ku(0.0),
	H_th_sigma(0.0)
	
	{

	//std::cout << "constructor " << anis_flag << "\t" << ianis_flag << std::endl;	
	// derived parameters
	for(int i=0;i<100;i++){
		geometry_coords[i][0]=0.0;
		geometry_coords[i][1]=0.0;
	}	
	// array variables
	for(int i=0;i<mp::max_materials;i++){
		Jij_matrix_SI[i]=0.0;
		Jij_matrix[i]=0.0;
		intermixing[i]=0.0;
		alloy[i]=0.0;
	}
	initial_spin[0]=0.0;
	initial_spin[0]=0.0;
	initial_spin[0]=1.0;

}

int materials_t::print(){

	std::cout << "----------------------------------------------------------------" << std::endl;
	std::cout << " Material " << name << std::endl;
	std::cout << "----------------------------------------------------------------" << std::endl;
	std::cout << "alpha          = " << alpha << std::endl;
	for(int j=0;j<num_materials;j++){
		//std::cout << " Jij_matrix_SI = " << material[i].Jij_matrix_SI[j] << "\t" << j << std::endl;
		std::cout << " Jij_matrix_SI[" << j << "] = " << Jij_matrix_SI[j] << std::endl;
	}
	std::cout << "mu_s_SI        = " << mu_s_SI << std::endl;
	std::cout << "Ku1_SI          = " << Ku1_SI << std::endl;
	std::cout << "gamma_rel      = " << gamma_rel << std::endl;
	
	return 0;
	
}	

}
