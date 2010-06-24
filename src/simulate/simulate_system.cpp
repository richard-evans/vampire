//====================================================================================================
//
//       				                    simulate_system
//
//  			 Subroutine to simulate an atomistic system with predefined integration scheme
//				 simulation time, temperature etc
//	 
//								Version 1.0 R Evans 02/10/2008
//
//==================================================================================================== 
#include "atoms.hpp"
#include "program.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include <iostream>

	//int LLG(const int);
	int set_LLG();
	int set_demag();
	//int curie_temperature(bool);
	//int program_hamr_run();

	//==========================================================
// Namespace simulation variables
//==========================================================
namespace sim{
	std::ofstream mag_file;
	int time;
	int total_time=1;
	int loop_time=1;
	int partial_time=1;
	int equilibration_time=0;
	
	double Tmax=300;
	double temperature;
	double H_applied=0.0;
	double H_vec[3]={0.0,0.0,1.0};
	double Hmin=-1.0; // T
	double Hmax=+1.0; // T
	double Hinc= 0.1; // T	

	int system_simulation_flags;
	int hamiltonian_simulation_flags[10];

	// derived variables
	
int initialise(){

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "initialise_system has been called" << std::endl;}

   //##########################################################
	//if system simulation is not initialised then
	//if(simulation_variables::time==0){
	//##########################################################
	//simulation_variables::total_time = 200000;
	//simulation_variables::loop_time = 5000;
	//simulation_variables::partial_time=0;
	//simulation_variables::temperature = 0.0;
	//material_parameters::H_th_sigma = sqrt(2.0*material_parameters::alpha*1.3806503e-23*simulation_variables::temperature/(material_parameters::mu_s_SI*material_parameters::gamma_SI*material_parameters::dt_SI));
	
	//for(int atom=0;atom<=atom_variables::num_atoms-1;atom++){
	//	atom_variables::spin_array[atom][0] = 0.0;
	//	atom_variables::spin_array[atom][1] = 1.0/sqrt(2.0);
	//	atom_variables::spin_array[atom][2] = -1.0/sqrt(2.0);
	//}

    for(int atom=0;atom<=atoms::num_atoms-1;atom++){
		atoms::x_spin_array[atom] = 0.0;
		atoms::y_spin_array[atom] = 0.1;
		atoms::z_spin_array[atom] = 0.9;
	}
	//std::cout.setf(std::ios::fixed,std::ios::floatfield);
	
  	sim::mag_file.open ("M_vs_T.txt");
      
	

	set_LLG();
	
	//######################################
	//}
	// close if not initialised
	//######################################
	return 0;

}

} // Namespace sim
//int mpi_init_halo_swap();
//int mpi_complete_halo_swap();
//int output_pov_file();

int simulate_system(){
    //program::curie_temperature(true);
	//program::hamr_run();
	//program::static_hysteresis();
	//program::two_temperature_pulse();
	program::bmark();
	//program::LLB_Boltzmann();
	//program::hysteresis();
	return 0;
}
  /*
//const int num_atoms = mpi_comms::num_core_atoms+mpi_comms::num_boundary_atoms;
const int num_atoms = atoms::num_atoms;


//std::cout << " core/boundary " << mpi_comms::num_core_atoms << "\t" << mpi_comms::num_boundary_atoms << std::endl;
statistics::inv_num_atoms = 1.0/double(num_atoms);
	//----------------------------------------
	// function prototypes
	//----------------------------------------

	// reset material types
	//for(int atom =0;atom<num_atoms;atom++){
		//if(mpi_create_variables::mpi_atom_comm_class_array[atom]==1){
		//	atoms::type_array[atom]=1;
		//}
	//}

	//int LLG();
	//int set_LLG();

	//int cs_create_crystal_structure(string,int[],int,int**,string*);
	//int cs_create_system_type(int[],int&,int**,string*);
	//int cs_create_neighbourlist(int[],int,int**,string*,int**);
	//int cs_set_atom_vars(int,int**,string*,int**);
	//----------------------------------------
	// Local variables for system simulation
	//----------------------------------------
	//double** total_spin_field_array;		// Total spin dependent fields
	//double** total_external_field_array;	// Total external fields
	//double temperature=0.0;

	//int atom;
	//ofstream mag_file;
	//int time;
	//int total_time;
	//int loop_time;
	//int partial_time;
	
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "simulate_system has been called" << std::endl;}

	//=============================================================
	//      Setup LLG arrays
	//=============================================================

	//set_LLG();
	initialise_system();
	//set_demag();

	//=============================================================
	//      Perform system simulation for 1 timestep
	//=============================================================
	for(int atom =0;atom<num_atoms;atom++){
		//if(mpi_generic::my_rank%2!=0){
		if(atoms::type_array[atom]==0){

			//atoms::x_spin_array[atom]=2.0*mtrandom::grnd()-1.0;
			//atoms::y_spin_array[atom]=2.0*mtrandom::grnd()-1.0;
			//atoms::z_spin_array[atom]=2.0*mtrandom::grnd()-1.0;
			atoms::x_spin_array[atom]=0.0;			
			atoms::y_spin_array[atom]=0.0;			
			atoms::z_spin_array[atom]=1.0;
		}
		//else{
		//	atoms::x_spin_array[atom]=0.0;			
		//	atoms::y_spin_array[atom]=0.0;			
		//	atoms::z_spin_array[atom]=1.0;
		//}
	}

	output_pov_file();
		
	for(simulation_variables::time=0;simulation_variables::time<=simulation_variables::total_time;simulation_variables::time++){
		//----------------------------------------
		// Alter system variables (temp, H etc)
		//----------------------------------------
		if(simulation_variables::partial_time==0){
			if(simulation_variables::time!=0){
				//std::cout << simulation_variables::time << "\t" << simulation_variables::temperature << "\t" << statistics::total_mag_m/statistics::data_counter << std::endl;
				simulation_variables::mag_file << simulation_variables::time << "\t" << simulation_variables::temperature << "\t" << statistics::total_mag_m/statistics::data_counter << std::endl;
			}
			simulation_variables::temperature = simulation_variables::temperature + 25.0;
			material_parameters::H_th_sigma = sqrt(2.0*material_parameters::alpha*1.38e-23*simulation_variables::temperature/(material_parameters::mu_s_SI*material_parameters::gamma_SI*material_parameters::dt_SI));
			
	statistics::total_mag_m=0.0;
	statistics::data_counter=0.0;

			//std::cout << temperature << std::endl;
		}
		//double reduced_time = double(simulation_variables::time)*0.1; // - 20000.0;
		
		//simulation_variables::temperature = 300.0+700.0*exp(-reduced_time*reduced_time/10000.0);
		//material_parameters::H_th_sigma = sqrt(2.0*material_parameters::alpha*1.38e-23*simulation_variables::temperature/(material_parameters::mu_s_SI*material_parameters::gamma_SI*material_parameters::dt_SI));

		//if(simulation_variables::time>100000){
		//	simulation_variables::H_applied=double(simulation_variables::time-100000)*(-1.0e-6);
		//}
		//else{
			simulation_variables::H_applied=0.0;
		//}

		#ifdef MPICF
		LLG_mpi();
		#else
		LLG();
		#endif
		//system("sleep 2");

		// Calculate mag_m, mag
		statistics::mag[0]=0.0;
		statistics::mag[1]=0.0;
		statistics::mag[2]=0.0;
		statistics::mag_m =0.0;
		for(int atom=0;atom<num_atoms;atom++){
			statistics::mag[0] = statistics::mag[0] + atoms::x_spin_array[atom];
			statistics::mag[1] = statistics::mag[1] + atoms::y_spin_array[atom];
			statistics::mag[2] = statistics::mag[2] + atoms::z_spin_array[atom];
		}
		
		statistics::mag[0] = statistics::mag[0]*statistics::inv_num_atoms;
		statistics::mag[1] = statistics::mag[1]*statistics::inv_num_atoms;
		statistics::mag[2] = statistics::mag[2]*statistics::inv_num_atoms;

		statistics::mag_m = sqrt(statistics::mag[0]*statistics::mag[0]+statistics::mag[1]*statistics::mag[1] + statistics::mag[2]*statistics::mag[2]);
				
		//if(simulation_variables::partial_time>1000){
		//	statistics::total_mag_m+=statistics::mag_m;
		//	statistics::data_counter+=1.0;
		//}
		
		//statistics::mag_m = sqrt(pow(statistics::mag[0],2) + pow(statistics::mag[1],2) +pow(statistics::mag[2],2));
		if((simulation_variables::time % -1)==0){
		simulation_variables::mag_file << simulation_variables::time << "\t" << 
								statistics::mag[0] << "\t" << 
								statistics::mag[1] << "\t" << 
								statistics::mag[2] << "\t" << 
								statistics::mag_m  << "\t" << 
								simulation_variables::temperature<< std::endl;
		}
		if((simulation_variables::time % 1000)==0){
			output_pov_file();
		}
		if((simulation_variables::time % 10)==0){
	
		#ifdef MPICF
			MPI::COMM_WORLD.Allreduce(&statistics::mag[0],&statistics::mag[0],1,MPI_DOUBLE,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(&statistics::mag[1],&statistics::mag[1],1,MPI_DOUBLE,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(&statistics::mag[2],&statistics::mag[2],1,MPI_DOUBLE,MPI_SUM);
			MPI::COMM_WORLD.Allreduce(&statistics::mag_m,&statistics::mag_m,1,MPI_DOUBLE,MPI_SUM);
		#endif
		if(mpi_generic::my_rank==0){
		std::cout << simulation_variables::time << "\t" << 
							statistics::mag[0]/double(mpi_generic::num_processors) << "\t" << 
							statistics::mag[1]/double(mpi_generic::num_processors) << "\t" << 
							statistics::mag[2]/double(mpi_generic::num_processors) << "\t" << 
							simulation_variables::H_applied << "\t" <<
							//simulation_variables::temperature << "\t" <<
							statistics::mag_m/double(mpi_generic::num_processors)  << std::endl;
		vmag << simulation_variables::time << "\t" << 
							statistics::mag[0]/double(mpi_generic::num_processors) << "\t" << 
							statistics::mag[1]/double(mpi_generic::num_processors) << "\t" << 
							statistics::mag[2]/double(mpi_generic::num_processors) << "\t" << 
							simulation_variables::H_applied << "\t" <<
							//simulation_variables::temperature << "\t" <<
							statistics::mag_m/double(mpi_generic::num_processors)  << std::endl;
		}
		}
		// increment partial timer and reset if end of loop is reached
		//simulation_variables::partial_time++;
		//if(simulation_variables::partial_time==simulation_variables::loop_time)simulation_variables::partial_time=0;
	}

	if(1==0){
		ofstream spin_file;
		spin_file.open ("spins.dat");
		spin_file << atoms::num_atoms << std::endl;
		spin_file << "" << std::endl;
	  	
	  	for(int atom=0; atom<atoms::num_atoms; atom++){
	  		spin_file << material_parameters::material[atoms::type_array[atom]].element << "\t" << 
	  					atoms::x_coord_array[atom]*material_parameters::lattice_space_conversion[0] << "\t" << 
	  					atoms::y_coord_array[atom]*material_parameters::lattice_space_conversion[1] << "\t" << 
	  					atoms::z_coord_array[atom]*material_parameters::lattice_space_conversion[2] << "\t" <<
	  					atoms::x_spin_array[atom] << "\t" << 
	  					atoms::y_spin_array[atom] << "\t" << 
	  					atoms::z_spin_array[atom] << "\t" << std::endl;
	  	}
	
		spin_file.close();
	}

	//simulation_variables::time++;
	//if(simulation_variables::time>=simulation_variables::total_time) exit(0);
	//#################################
	// if simulation ended
	//#################################
	//mag_file.close();
  	//std::cout << "System Simulated Successfully" << std::endl;
  	*/
 	

