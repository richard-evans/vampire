//====================================================================================
//
//														mpi_generic
//
//									Functions for basic mpi functionality
//
//										Version 1.0 R Evans 07/08/2009
//
//====================================================================================

#include "errors.hpp"
#include "vmpi.hpp"
#include <iostream>

#ifdef MPICF
namespace vmpi{

int initialise(){
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
	if(err::check==true){std::cout << "initialise_mpi has been called" << std::endl;}

	int resultlen;

	// Initialise MPI
	MPI::Init();

	// Get number of processors and rank
	vmpi::my_rank = MPI::COMM_WORLD.Get_rank();
	vmpi::num_processors = MPI::COMM_WORLD.Get_size();
	MPI::Get_processor_name(vmpi::hostname, resultlen);

	// Start MPI Timer
	vmpi::start_time=MPI_Wtime();

	return EXIT_SUCCESS;	
}

int hosts(){
	//====================================================================================
	//
	//													mpi_hosts
	//
	//										Print MPI hostnames to screen
	//
	//										Version 1.0 R Evans 07/08/2009
	//
	//====================================================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "mpi_hosts has been called" << std::endl;}

	// Wait for all processors
	//MPI::COMM_WORLD.Barrier();

   //for(int i=0;i<vmpi::num_processors;i++){
   //	if(vmpi::my_rank==i){
	if(vmpi::num_processors<=512){
			// Output rank, num_procs, hostname to screen
			std::cout << "Processor " << vmpi::my_rank+1 << " of " << vmpi::num_processors;
			std::cout << " online on host " << vmpi::hostname << std::endl;
			//	}
			//MPI::COMM_WORLD.Barrier();
	   //}
	}
	// Wait for all processors
			//MPI::COMM_WORLD.Barrier();

	return EXIT_SUCCESS;	
}

int finalise(){
	//====================================================================================
	//
	//												finalise_mpi
	//
	//									Finalise MPI and output wall time
	//
	//										Version 1.0 R Evans 07/08/2009
	//
	//====================================================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "finalise_mpi has been called" << std::endl;}

	// Wait for all processors
   MPI::COMM_WORLD.Barrier();

	// Stop MPI Timer and output to screen
	vmpi::end_time=MPI_Wtime();
	if(vmpi::my_rank==0){
		std::cout << "MPI Simulation Time: " << vmpi::end_time-vmpi::start_time << std::endl;
	}

	// Finalise MPI
	MPI::Finalize();

	return EXIT_SUCCESS;	
}

} // end of namespace vmpi
#endif 

