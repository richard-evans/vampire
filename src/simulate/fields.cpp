//====================================================================================================
//
//       				                    	Fields
//
//  			 		Subroutines to calculate fields for the hamiltonian
//	 
//									Version 1.0 R Evans 20/10/2008
//
//==================================================================================================== 
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "demag.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

#include <algorithm>
#include <iostream>

//========================
//function prototypes
//========================

int calculate_exchange_fields(const int,const int);
int calculate_uniaxial_anis_fields(const int,const int);
int calculate_cubic_anis_fields(const int,const int);
int calculate_applied_fields(const int,const int);
int calculate_thermal_fields(const int,const int);
int calculate_dipolar_fields(const int,const int);
int demag_field_update();

int calculate_spin_fields(const int start_index,const int end_index){
	//======================================================
	// 		Subroutine to calculate spin dependent fields
	//
	//			Version 1.0 R Evans 20/10/2008
	//======================================================

	//const int num_atoms = atoms::num_atoms;

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_spin_fields has been called" << std::endl;}
	
	// Initialise Total Spin Fields to zero
	//fill (atoms::x_total_spin_field_array.begin(),atoms::x_total_spin_field_array.end(),0.0);
	fill (atoms::x_total_spin_field_array.begin()+start_index,atoms::x_total_spin_field_array.begin()+end_index,0.0);
        fill (atoms::y_total_spin_field_array.begin()+start_index,atoms::y_total_spin_field_array.begin()+end_index,0.0);
        fill (atoms::z_total_spin_field_array.begin()+start_index,atoms::z_total_spin_field_array.begin()+end_index,0.0);
	//for(int atom=start_index;atom<end_index;atom++){
	//	atoms::x_total_spin_field_array[atom] = 0.0;
	//atoms::y_total_spin_field_array[atom] = 0.0;
	//	atoms::z_total_spin_field_array[atom] = 0.0;
	//}
	
	//const int num_atoms = atoms::num_atoms;
	//std::vector<double> spin_array(3*num_atoms);
        //for(int i=0;i<num_atoms;i++){
        //  spin_array[3*i+0]=atoms::x_spin_array[i];
        //  spin_array[3*i+1]=atoms::y_spin_array[i];
        //  spin_array[3*i+2]=atoms::z_spin_array[i];
        //}
	// Exchange Fields
	if(sim::hamiltonian_simulation_flags[0]==1) calculate_exchange_fields(start_index,end_index);
	
	// Anisotropy Fields
	if(sim::hamiltonian_simulation_flags[1]==1) calculate_uniaxial_anis_fields(start_index,end_index);
	if(sim::hamiltonian_simulation_flags[1]==2) calculate_cubic_anis_fields(start_index,end_index);
	//if(sim::hamiltonian_simulation_flags[1]==3) calculate_local_anis_fields();
	
	// Spin Dependent Extra Fields
	//if(sim::hamiltonian_simulation_flags[4]==1) calculate_??_fields();
	
	return 0;
}

int calculate_external_fields(const int start_index,const int end_index){
	//======================================================
	// 		Subroutine to calculate external fields
	//
	//			Version 1.0 R Evans 20/10/2008
	//======================================================
	//const int num_atoms = atoms::num_atoms;

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_external_fields has been called" << std::endl;}

	// Initialise Total External Fields to zero
	//for(int atom=start_index;atom<end_index;atom++){
	//	atoms::x_total_external_field_array[atom] = 0.0;
	//	atoms::y_total_external_field_array[atom] = 0.0;
	//	atoms::z_total_external_field_array[atom] = 0.0;	
	//}
        fill (atoms::x_total_external_field_array.begin()+start_index,atoms::x_total_external_field_array.begin()+end_index,0.0);
        fill (atoms::y_total_external_field_array.begin()+start_index,atoms::y_total_external_field_array.begin()+end_index,0.0);
        fill (atoms::z_total_external_field_array.begin()+start_index,atoms::z_total_external_field_array.begin()+end_index,0.0);
	
	// Thermal Fields
        if(sim::hamiltonian_simulation_flags[3]==1) calculate_thermal_fields(start_index,end_index);

	// Applied Fields
	if(sim::hamiltonian_simulation_flags[2]==1) calculate_applied_fields(start_index,end_index);
	
	// Thermal Fields
	//if(sim::hamiltonian_simulation_flags[3]==1) calculate_thermal_fields(start_index,end_index);

	// Dipolar Fields
	if(sim::hamiltonian_simulation_flags[4]==1) calculate_dipolar_fields(start_index,end_index);
	
	return 0;
}

int calculate_exchange_fields(const int start_index,const int end_index){
	//======================================================
	// 		Subroutine to calculate exchange fields
	//
	//			Version 1.0 R Evans 20/10/2008
	//======================================================

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_exchange_fields has been called" << std::endl;}

	const int prank=1;
	//const int num_atoms = atoms::num_atoms;

	//std::vector<double> spin_array(3*num_atoms);
	//for(int i=0;i<num_atoms;i++){
	//	spin_array[3*i+0]=atoms::x_spin_array[i];
	//	spin_array[3*i+1]=atoms::y_spin_array[i];
	//	spin_array[3*i+2]=atoms::z_spin_array[i];
	//}
		//if(vmpi::my_rank==prank){
		//std::cout << "--------------------------------------------------------------------------------------" << std::endl;
		//}
	for(int atom=start_index;atom<end_index;atom++){
		const int imaterial=atoms::type_array[atom];
		//double tot[3]={0.0,0.0,0.0};
		for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){
			
			const int natom = atoms::neighbour_list_array[nn];
			const int jmaterial=atoms::type_array[natom]; // 1D array?
			const double Jij = material_parameters::material[imaterial].Jij_matrix[jmaterial];
			atoms::x_total_spin_field_array[atom] -= Jij*atoms::x_spin_array[natom];
			atoms::y_total_spin_field_array[atom] -= Jij*atoms::y_spin_array[natom];
			atoms::z_total_spin_field_array[atom] -= Jij*atoms::z_spin_array[natom];
			//tot[0]-=Jij*spin_array[3*natom+0];
			//tot[1]-=Jij*spin_array[3*natom+1];
			//tot[2]-=Jij*spin_array[3*natom+2];
			//if(vmpi::my_rank==prank){
			//std::cout << "\t" << atom << " " << natom << " " << atoms::x_spin_array[natom] << " " << atoms::y_spin_array[natom] << " " << atoms::z_spin_array[natom] << std::endl;
			//}
			}
		//atoms::x_total_spin_field_array[atom]+=tot[0];
		//atoms::y_total_spin_field_array[atom]+=tot[1];
		//atoms::z_total_spin_field_array[atom]+=tot[2];
		if(vmpi::my_rank==prank){
		//std::cout << atom << "\texchange fields\t" << atoms::x_total_spin_field_array[atom] << "\t";
		//std::cout << atoms::y_total_spin_field_array[atom] << "\t";
		//std::cout << atoms::z_total_spin_field_array[atom] << std::endl;
		//std::cout << "\t=================================================================" << std::endl;
		//std::cout << atom << "\texchange fields\t" << tot[0] << "\t";
		//std::cout << tot[1] << "\t";
		//std::cout << tot[2] << std::endl;
		}
		//std::cin.get();
	}
	//system("sleep 2");
	//exit(0);
	
	return EXIT_SUCCESS;
	}

int calculate_uniaxial_anis_fields(const int start_index,const int end_index){
	//======================================================
	// 	Subroutine to calculate uniaxial anisotropy fields
	//
	//			Version 1.0 R Evans 20/10/2008
	//======================================================

	//const int num_atoms = atoms::num_atoms;

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_uniaxial_anis_fields has been called" << std::endl;}

	for(int atom=start_index;atom<end_index;atom++){
		const int imaterial=atoms::type_array[atom];
		atoms::z_total_spin_field_array[atom] -= 2.0*mp::material[imaterial].Ku*atoms::z_spin_array[atom];
		//std::cout << atom << "\tanisotropy fields\t" << mp::material[imaterial].Ku << "\t";
		//std::cout << atoms::x_total_spin_field_array[atom] << "\t";
		//std::cout << atoms::y_total_spin_field_array[atom] << "\t";
		//std::cout << atoms::z_total_spin_field_array[atom] << std::endl;
	}
	return 0;
}

int calculate_cubic_anis_fields(const int start_index,const int end_index){
	//======================================================
	// 	Subroutine to calculate cubic anisotropy fields
	//
	//			Version 1.0 R Evans 20/10/2008
	//======================================================
	//const int num_atoms = atom_variables::num_atoms;
	//for(int atom=0;atom<num_atoms;atom++){
	//	atom_fields::total_spin_field_array[atom][2] += 2.0*material_parameters::Ku*atom_variables::spin_array[atom][2];
	//}
  for(int i=start_index;i<end_index;i++){
  }
	return 0;
}

int calculate_applied_fields(const int start_index,const int end_index){
	//======================================================
	// 	Subroutine to calculate applied fields
	//
	//			Version 1.0 R Evans 20/10/2008
	//======================================================
	//const int num_atoms = atoms::num_atoms;

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_applied_fields has been called" << std::endl;}

	for(int atom=start_index;atom<end_index;atom++){
		atoms::x_total_external_field_array[atom] += sim::H_vec[0]*sim::H_applied;
		atoms::y_total_external_field_array[atom] += sim::H_vec[1]*sim::H_applied;
		atoms::z_total_external_field_array[atom] += sim::H_vec[2]*sim::H_applied;

		//std::cout << atom << "\tapplied fields\t" << sim::H_vec[0]*sim::H_applied << "\t";
		//std::cout << sim::H_vec[1]*sim::H_applied << "\t";
		//std::cout << sim::H_vec[2]*sim::H_applied << std::endl;
	}
	return 0;
}

int calculate_thermal_fields(const int start_index,const int end_index){
	//======================================================
	// 		Subroutine to calculate thermal fields
	//
	//			Version 1.0 R Evans 20/10/2009
	//======================================================

	//const int num_atoms = atoms::num_atoms;
	const double sqrt_T=sqrt(sim::temperature);
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_thermal_fields has been called" << std::endl;}

        //const int num_atoms = atoms::num_atoms;

        //std::vector<double> tfield_array(3*num_atoms);
	generate (atoms::x_total_external_field_array.begin()+start_index,atoms::x_total_external_field_array.begin()+end_index, mtrandom::gaussian);
        generate (atoms::y_total_external_field_array.begin()+start_index,atoms::y_total_external_field_array.begin()+end_index, mtrandom::gaussian);
        generate (atoms::z_total_external_field_array.begin()+start_index,atoms::z_total_external_field_array.begin()+end_index, mtrandom::gaussian);
	//generate (tfield_array.begin()+3*start_index,tfield_array.begin()+3*end_index, mtrandom::gaussian);

	for(int atom=start_index;atom<end_index;atom++){
		const int imaterial=atoms::type_array[atom];
		const double H_th_sigma = sqrt_T*material_parameters::material[imaterial].H_th_sigma;
		atoms::x_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
		atoms::y_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
		atoms::z_total_external_field_array[atom] *= H_th_sigma; //*mtrandom::gaussian();
		//for(int i=0;i<3;i++){
		//  tfield_array[3*atom+i]*=H_th_sigma;
		//}
		//atoms::x_total_external_field_array[atom]+=tfield_array[3*atom+0];
                //atoms::y_total_external_field_array[atom]+=tfield_array[3*atom+1];
		//atoms::z_total_external_field_array[atom]+=tfield_array[3*atom+2];
 
	}

	//for(int atom=start_index;atom<end_index;atom++){
	//  atoms::x_total_external_field_array[atom]+=tfield_array[3*atom+0];
	//  atoms::y_total_external_field_array[atom]+=tfield_array[3*atom+1];
	//  atoms::z_total_external_field_array[atom]+=tfield_array[3*atom+2];
	//}

	return EXIT_SUCCESS;
}

int calculate_dipolar_fields(const int start_index,const int end_index){
	//======================================================
	// 		Subroutine to calculate dipolar fields
	//
	//			Version 1.0 R Evans 02/11/2009
	//======================================================

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_dipolar_fields has been called" << std::endl;}
	
	// Check for update of dipolar fields
	if(demag::update_counter%demag::update_rate==0){
		demag_field_update();
		demag::update_counter=0;
	}

	// Add dipolar fields
	for(int atom=start_index;atom<end_index;atom++){
		atoms::x_total_external_field_array[atom] += atoms::x_dipolar_field_array[atom];
		atoms::y_total_external_field_array[atom] += atoms::y_dipolar_field_array[atom];
		atoms::z_total_external_field_array[atom] += atoms::z_dipolar_field_array[atom];
	}
	// Update counter
	demag::update_counter++;

	return 0;
}
