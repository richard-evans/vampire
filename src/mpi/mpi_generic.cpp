//====================================================================================
//
//														mpi_generic
//
//									Functions for basic mpi functionality
//
//										Version 1.0 R Evans 07/08/2009
//
//====================================================================================
//
//		Locally allocated variables: 	
//
//=====================================================================================

#include "public.hpp"
//#include <mpi.h>
#include "vmpi.hpp"
#include <iostream>
//using namespace std;
//==========================================================
// Namespace mpi_generic
//==========================================================
/*namespace mpi_generic{
	int mpi_mode=0;
	int my_rank=0;
	int num_processors=1;
	char hostname[20];
	double start_time;
	double end_time;
}
//==========================================================
// Namespace mpi_create_variables
//==========================================================
namespace mpi_create_variables{
	int int_mpi_offset[3];
	int int_mpi_halo_min[3];
	int int_mpi_halo_max[3];
	int int_mpi_core_min[3];
	int int_mpi_core_max[3];
	int mpi_interaction_range=0;
	double mpi_offset[3];
	double mpi_full_system_dimensions[3];
	double mpi_system_min[3];	
	double mpi_system_max[3];
	string decomposition_type;
	// Defines local(0), boundary(1), halo(2) atoms
	valarray<int> mpi_atom_comm_class_array(0,1);
	valarray<int> mpi_atom_location_array(0,1);
	valarray<int> mpi_atom_global_coord_array(0,1);

	valarray<int> mpi_global_atom_number_array(0,1);

	bool mpi_comms_identify=false;	// identifies halo and boundary atoms in xyz output

}
//==========================================================
// Namespace mpi_stats
//==========================================================
namespace mpi_stats{
}
//==========================================================
// Namespace mpi_comms
//==========================================================
namespace mpi_comms{

	int num_core_atoms=0;
	int num_boundary_atoms=0;
	int num_halo_atoms=0;

	int num_boundary_swaps=0;
	int num_halo_swaps=0;
	int num_messages=0;
	
	#ifdef MPICF
		valarray<MPI::Request> requests(1);
		valarray<MPI::Status> stati(1);
	#endif

	valarray<int> send_atom_translation_array(0,1);
	valarray<int> send_start_index_array(0,1);
	valarray<int> send_num_array(0,1);
	valarray<double> send_spin_data_array(0.0,1);

	valarray<int> recv_atom_translation_array(0,1);
	valarray<int> recv_start_index_array(0,1);
	valarray<int> recv_num_array(0,1);
	valarray<double> recv_spin_data_array(0.0,1);

	}*/
#ifdef MPICF
int initialise_mpi(){
	//====================================================================================
	//
	//												initialise_mpi
	//
	//									Startup MPI and Output host information
	//
	//										Version 1.0 R Evans 16/07/2009
	//
	//====================================================================================
	//
	//		Locally allocated variables: 	
	//
	//====================================================================================
	// C++ MPI Binding Reference:
	// https://computing.llnl.gov/tutorials/mpi/
	//====================================================================================

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(error_checking::error_check==true){std::cout << "initialise_mpi has been called" << std::endl;}

	int resultlen;

	//-----------------------
	// Initialise MPI
	//-----------------------
	MPI::Init();

	//------------------------------------
	// Get number of processors and rank
	//------------------------------------

	//mpi_generic::my_rank = MPI::COMM_WORLD.Get_rank();
	//mpi_generic::num_processors = MPI::COMM_WORLD.Get_size();
	//MPI::Get_processor_name(mpi_generic::hostname, resultlen);

	vmpi::my_rank = MPI::COMM_WORLD.Get_rank();
	vmpi::num_processors = MPI::COMM_WORLD.Get_size();
	MPI::Get_processor_name(vmpi::hostname, resultlen);

	//------------------------------------
	// Start MPI Timer
	//------------------------------------
	vmpi::start_time=MPI_Wtime();
	//mpi_generic::start_time=MPI_Wtime();

	return 0;	
}

int mpi_hosts(){
	//====================================================================================
	//
	//													mpi_hosts
	//
	//										Print MPI hostnames to screen
	//
	//										Version 1.0 R Evans 07/08/2009
	//
	//====================================================================================
	//
	//		Locally allocated variables: 	
	//
	//====================================================================================

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(error_checking::error_check==true){std::cout << "mpi_hosts has been called" << std::endl;}

	//-------------------------
	// Wait for all processors
	//-------------------------
   MPI::COMM_WORLD.Barrier();

	for(int i=0;i<vmpi::num_processors;i++){
		if(vmpi::my_rank==i){
			//--------------------------------------------
			// Output rank, num_procs, hostname to screen
			//--------------------------------------------
			std::cout << "Processor " << vmpi::my_rank+1 << " of " << vmpi::num_processors;
			std::cout << " online on host " << vmpi::hostname << std::endl;
		}
	   MPI::COMM_WORLD.Barrier();
	}

	//-------------------------
	// Wait for all processors
	//-------------------------
   MPI::COMM_WORLD.Barrier();

	return 0;	
}

int finalise_mpi(){
	//====================================================================================
	//
	//												finalise_mpi
	//
	//									Finalise MPI and output wall time
	//
	//										Version 1.0 R Evans 07/08/2009
	//
	//====================================================================================
	//
	//		Locally allocated variables: 	
	//
	//====================================================================================

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(error_checking::error_check==true){std::cout << "finalise_mpi has been called" << std::endl;}

	//-------------------------
	// Wait for all processors
	//-------------------------
   MPI::COMM_WORLD.Barrier();

	//-------------------------------------
	// Stop MPI Timer and output to screen
	//-------------------------------------
	vmpi::end_time=MPI_Wtime();
	if(vmpi::my_rank==0){
		std::cout << "MPI Wall Time: " << vmpi::end_time-vmpi::start_time << std::endl;
	}

	//-----------------------
	// Finalise MPI
	//-----------------------
	MPI::Finalize();

	return 0;	
}

#endif 

