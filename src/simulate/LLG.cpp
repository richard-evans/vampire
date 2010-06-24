//====================================================================================================
//
//       				                    	LLG
//
//  			 Subroutine to simulate an atomistic system with LLG integration scheme
//	 
//									Version 1.0 R Evans 02/10/2008
//
//==================================================================================================== 
/// \file LLG.cpp
/// Contains LLG namespace and serial version of the integrator
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "vmpi.hpp"

#include <iostream>
int calculate_spin_fields(const int,const int);
int calculate_external_fields(const int,const int);
int LLG_serial_heun(const int);
int LLG_relax_serial_heun(const int);
int LLG_mpi(const int);

/// \namespace LLG_arrays
/// \brief Defines array variables for storage of temporary data during integration
/// 
/// Arrays are initialised to zero. set_LLG is called after system generation to set size of arrays to the number of atoms
namespace LLG_arrays{
	using std::valarray;
//==========================================================
// Namespace to store persistant LLG integration arrays
//==========================================================

	valarray <double> x_euler_array;
	valarray <double> y_euler_array;	
	valarray <double> z_euler_array;

	valarray <double> x_heun_array;	
	valarray <double> y_heun_array;	
	valarray <double> z_heun_array;

	valarray <double> x_spin_storage_array;	
	valarray <double> y_spin_storage_array;	
	valarray <double> z_spin_storage_array;

	valarray <double> x_initial_spin_array;	
	valarray <double> y_initial_spin_array;	
	valarray <double> z_initial_spin_array;

	bool LLG_set=false; ///< Flag to define state of LLG arrays (initialised/uninitialised)

}
/// Initialises LLG Arrays in namespace LLG_arrays to be of size num_atoms and initialised to 0.0. Also sets LLG_set=true.
int set_LLG(){
	//======================================================
	// Subroutine to allocate LLG variables
	//======================================================

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "set_LLG has been called" << std::endl;}

	if(LLG_arrays::LLG_set==true){
		std::cerr << "Warning - LLG arrays have been reinitialised" << std::endl;
		std::cerr << "\t Check for additional calls to set_LLG()" << std::endl;
	}

	LLG_arrays::x_spin_storage_array.resize(atoms::num_atoms,0.0);
	LLG_arrays::y_spin_storage_array.resize(atoms::num_atoms,0.0);
	LLG_arrays::z_spin_storage_array.resize(atoms::num_atoms,0.0);

	LLG_arrays::x_initial_spin_array.resize(atoms::num_atoms,0.0);
	LLG_arrays::y_initial_spin_array.resize(atoms::num_atoms,0.0);
	LLG_arrays::z_initial_spin_array.resize(atoms::num_atoms,0.0);

	LLG_arrays::x_euler_array.resize(atoms::num_atoms,0.0);
	LLG_arrays::y_euler_array.resize(atoms::num_atoms,0.0);
	LLG_arrays::z_euler_array.resize(atoms::num_atoms,0.0);

	LLG_arrays::x_heun_array.resize(atoms::num_atoms,0.0);
	LLG_arrays::y_heun_array.resize(atoms::num_atoms,0.0);
	LLG_arrays::z_heun_array.resize(atoms::num_atoms,0.0);

	LLG_arrays::LLG_set=true;

	/*

	if(LLG_arrays::spin_storage_array==NULL){
		try{LLG_arrays::spin_storage_array=new double*[atom_variables::num_atoms];
    		for(int i=0; i<atom_variables::num_atoms ; i++)LLG_arrays::spin_storage_array[i]=new double[3];}
  		catch(...){cerr << "error allocating spin_storage_array" << std::endl;exit(1);} 
	}
	else{cerr << "WARNING - Attempted reallocation of LLG_arrays::spin_storage_array " << LLG_arrays::spin_storage_array << std::endl;}

	if(LLG_arrays::initial_spin_array==NULL){
		try{LLG_arrays::initial_spin_array=new double*[atom_variables::num_atoms];
			for(int i=0; i<atom_variables::num_atoms ; i++)LLG_arrays::initial_spin_array[i]=new double[3];}
		catch(...){cerr << "error allocating initial_spin_array" << std::endl;exit(1);} 
	}
	else{cerr << "WARNING - Attempted reallocation of LLG_arrays::initial_spin_array " << LLG_arrays::initial_spin_array << std::endl;}

	if(LLG_arrays::euler_array==NULL){
  		try{LLG_arrays::euler_array=new double*[atom_variables::num_atoms];
    		for(int i=0; i<atom_variables::num_atoms ; i++)LLG_arrays::euler_array[i]=new double[3];}
  		catch(...){cerr << "error allocating euler_array" << std::endl;exit(1);} 
	}
	else{cerr << "WARNING - Attempted reallocation of LLG_arrays::euler_array " << LLG_arrays::euler_array << std::endl;}

	if(LLG_arrays::heun_array==NULL){
  		try{LLG_arrays::heun_array=new double*[atom_variables::num_atoms];
   	 	for(int i=0; i<atom_variables::num_atoms ; i++)LLG_arrays::heun_array[i]=new double[3];}
  		catch(...){cerr << "error allocating heun_array" << std::endl;exit(1);} 
	}
	else{cerr << "WARNING - Attempted reallocation of LLG_arrays::heun_array " << LLG_arrays::heun_array << std::endl;}
	*/
  	return 0;
}

namespace sim{
/// Master LLG Function - dispatches code path to desired LLG routine
/// \f$ \frac{\partial S}{\partial t} \f$
int LLG(const int num_steps){

   //----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "LLG has been called" << std::endl;}

	#ifdef MPICF
		LLG_mpi(num_steps);
	#else
		LLG_serial_heun(num_steps);
	#endif
	
	return 0;
}

}

/// Performs serial Heun integration of the Landau-Lifshitz-Gilbert Equation of motion
int LLG_serial_heun(const int num_steps){

	using namespace LLG_arrays;
	
	//----------------------------------------
	// Local variables for system integration
	//----------------------------------------
	const int num_atoms = atoms::num_atoms;
	double xyz[3];		// Local Delta Spin Components
	double S_new[3];	// New Local Spin Moment
	double mod_S;		// magnitude of spin moment 

	for(int t=0;t<num_steps;t++){

		//----------------------------------------
		// Store initial spin positions
		//----------------------------------------
		
		for(int atom=0;atom<num_atoms;atom++){
			x_initial_spin_array[atom] = atoms::x_spin_array[atom];
			y_initial_spin_array[atom] = atoms::y_spin_array[atom];
			z_initial_spin_array[atom] = atoms::z_spin_array[atom];
		}
			
		//----------------------------------------
		// Set field arrays to zero
		//----------------------------------------
		//----------------------------------------
		// Calculate fields
		//----------------------------------------	
		
		calculate_spin_fields(0,num_atoms);
		calculate_external_fields(0,num_atoms);
		//----------------------------------------
		// Calculate Euler Step
		//----------------------------------------	
		
		for(int atom=0;atom<num_atoms;atom++){

			const int imaterial=atoms::type_array[atom];
			const double one_oneplusalpha_sq = material_parameters::material[imaterial].one_oneplusalpha_sq;
			const double alpha_oneplusalpha_sq = material_parameters::material[imaterial].alpha_oneplusalpha_sq;

			// Store local spin in Sand local field in H
			const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
			const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
										atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
										atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};

			// Calculate Delta S
			xyz[0]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
			xyz[1]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
			xyz[2]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

			// Store dS in euler array
			x_euler_array[atom]=xyz[0];
			y_euler_array[atom]=xyz[1];
			z_euler_array[atom]=xyz[2];

			// Calculate Euler Step
			S_new[0]=S[0]+xyz[0]*material_parameters::dt;
			S_new[1]=S[1]+xyz[1]*material_parameters::dt;
			S_new[2]=S[2]+xyz[2]*material_parameters::dt;
			
			// Normalise Spin Length
			mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);
			
			S_new[0]=S_new[0]*mod_S;
			S_new[1]=S_new[1]*mod_S;
			S_new[2]=S_new[2]*mod_S;

			//Writing of Spin Values to Storage Array
			x_spin_storage_array[atom]=S_new[0];
			y_spin_storage_array[atom]=S_new[1];
			z_spin_storage_array[atom]=S_new[2];		
		}
		
		//----------------------------------------
		// Copy new spins to spin array
		//----------------------------------------
		for(int atom=0;atom<num_atoms;atom++){
			atoms::x_spin_array[atom]=x_spin_storage_array[atom];
			atoms::y_spin_array[atom]=y_spin_storage_array[atom];
			atoms::z_spin_array[atom]=z_spin_storage_array[atom];
		}
		
		//----------------------------------------
		// Recalculate spin dependent fields
		//----------------------------------------
		calculate_spin_fields(0,num_atoms);
		
		//----------------------------------------
		// Calculate Heun Gradients
		//----------------------------------------	
		
		for(int atom=0;atom<num_atoms;atom++){

			const int imaterial=atoms::type_array[atom];;
			const double one_oneplusalpha_sq = material_parameters::material[imaterial].one_oneplusalpha_sq;
			const double alpha_oneplusalpha_sq = material_parameters::material[imaterial].alpha_oneplusalpha_sq;

			// Store local spin in Sand local field in H
			const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
			const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
										atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
										atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};

			// Calculate Delta S
			xyz[0]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
			xyz[1]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
			xyz[2]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

			// Store dS in heun array
			x_heun_array[atom]=xyz[0];
			y_heun_array[atom]=xyz[1];
			z_heun_array[atom]=xyz[2];
		}

		//----------------------------------------
		// Calculate Heun Step
		//----------------------------------------	

		for(int atom=0;atom<num_atoms;atom++){
			S_new[0]=x_initial_spin_array[atom]+material_parameters::half_dt*(x_euler_array[atom]+x_heun_array[atom]);
			S_new[1]=y_initial_spin_array[atom]+material_parameters::half_dt*(y_euler_array[atom]+y_heun_array[atom]);
			S_new[2]=z_initial_spin_array[atom]+material_parameters::half_dt*(z_euler_array[atom]+z_heun_array[atom]);
			
			// Normalise Spin Length
			mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);
			
			S_new[0]=S_new[0]*mod_S;
			S_new[1]=S_new[1]*mod_S;
			S_new[2]=S_new[2]*mod_S;

			//----------------------------------------
			// Copy new spins to spin array
			//----------------------------------------
			atoms::x_spin_array[atom]=S_new[0];
			atoms::y_spin_array[atom]=S_new[1];
			atoms::z_spin_array[atom]=S_new[2];
		}

	}
	
	return 0;
	}

namespace sim{
/// Master LLG_relax Function - dispatches code path to desired LLG routine
/// \f$ \frac{\partial S}{\partial t} \f$
int LLG_relax(const int num_steps){

   //----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "LLG_relax has been called" << std::endl;}

	#ifdef MPICF
		std::cerr << "MPI version of LLG_relax is yet to be implemented, exiting" << std::endl;
		exit(EXIT_FAILURE);
	#else
		LLG_relax_serial_heun(num_steps);
	#endif
	
	return 0;
}
}
/// Performs serial Heun integration of the relaxational part of the Landau-Lifshitz-Gilbert Equation of motion
int LLG_relax_serial_heun(const int num_steps){

	using namespace LLG_arrays;
	
	//----------------------------------------
	// Local variables for system integration
	//----------------------------------------
	const int num_atoms = atoms::num_atoms;
	double xyz[3];		// Local Delta Spin Components
	double S_new[3];	// New Local Spin Moment
	double mod_S;		// magnitude of spin moment 

	for(int t=0;t<num_steps;t++){

		//----------------------------------------
		// Store initial spin positions
		//----------------------------------------
		
		for(int atom=0;atom<num_atoms;atom++){
			x_initial_spin_array[atom] = atoms::x_spin_array[atom];
			y_initial_spin_array[atom] = atoms::y_spin_array[atom];
			z_initial_spin_array[atom] = atoms::z_spin_array[atom];
		}
			
		//----------------------------------------
		// Calculate fields
		//----------------------------------------	
		
		calculate_spin_fields(0,num_atoms);
		calculate_external_fields(0,num_atoms);

		//----------------------------------------
		// Calculate Euler Step
		//----------------------------------------	
		
		for(int atom=0;atom<num_atoms;atom++){

			const int imaterial=atoms::type_array[atom];
			const double alpha_oneplusalpha_sq = material_parameters::material[imaterial].alpha_oneplusalpha_sq;

			// Store local spin in Sand local field in H
			const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
			const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
										atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
										atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};

			// Calculate Delta S
			xyz[0]= (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
			xyz[1]= (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
			xyz[2]= (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

			// Store dS in euler array
			x_euler_array[atom]=xyz[0];
			y_euler_array[atom]=xyz[1];
			z_euler_array[atom]=xyz[2];

			// Calculate Euler Step
			S_new[0]=S[0]+xyz[0]*material_parameters::dt;
			S_new[1]=S[1]+xyz[1]*material_parameters::dt;
			S_new[2]=S[2]+xyz[2]*material_parameters::dt;
			
			// Normalise Spin Length
			mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);
			
			S_new[0]=S_new[0]*mod_S;
			S_new[1]=S_new[1]*mod_S;
			S_new[2]=S_new[2]*mod_S;

			//Writing of Spin Values to Storage Array
			x_spin_storage_array[atom]=S_new[0];
			y_spin_storage_array[atom]=S_new[1];
			z_spin_storage_array[atom]=S_new[2];		
		}
		
		//----------------------------------------
		// Copy new spins to spin array
		//----------------------------------------
		for(int atom=0;atom<num_atoms;atom++){
			atoms::x_spin_array[atom]=x_spin_storage_array[atom];
			atoms::y_spin_array[atom]=y_spin_storage_array[atom];
			atoms::z_spin_array[atom]=z_spin_storage_array[atom];
		}
		
			calculate_spin_fields(0,num_atoms);

		//----------------------------------------
		// Calculate Heun Gradients
		//----------------------------------------	
		
		for(int atom=0;atom<num_atoms;atom++){

			const int imaterial=atoms::type_array[atom];;
			const double alpha_oneplusalpha_sq = material_parameters::material[imaterial].alpha_oneplusalpha_sq;

			// Store local spin in Sand local field in H
			const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
			const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
										atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
										atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};

			// Calculate Delta S
			xyz[0]=(alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
			xyz[1]=(alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
			xyz[2]=(alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

			// Store dS in heun array
			x_heun_array[atom]=xyz[0];
			y_heun_array[atom]=xyz[1];
			z_heun_array[atom]=xyz[2];
		}

		//----------------------------------------
		// Calculate Heun Step
		//----------------------------------------	

		for(int atom=0;atom<num_atoms;atom++){
			S_new[0]=x_initial_spin_array[atom]+material_parameters::half_dt*(x_euler_array[atom]+x_heun_array[atom]);
			S_new[1]=y_initial_spin_array[atom]+material_parameters::half_dt*(y_euler_array[atom]+y_heun_array[atom]);
			S_new[2]=z_initial_spin_array[atom]+material_parameters::half_dt*(z_euler_array[atom]+z_heun_array[atom]);
			
			// Normalise Spin Length
			mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);
			
			S_new[0]=S_new[0]*mod_S;
			S_new[1]=S_new[1]*mod_S;
			S_new[2]=S_new[2]*mod_S;

			//----------------------------------------
			// Copy new spins to spin array
			//----------------------------------------
			atoms::x_spin_array[atom]=S_new[0];
			atoms::y_spin_array[atom]=S_new[1];
			atoms::z_spin_array[atom]=S_new[2];
		}

	}
	
	return 0;
}
