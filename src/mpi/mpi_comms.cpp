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
//====================================================================================
//
//														mpi_comms
//
//						Functions for performing halo swap during parallel execution
//
//										Version 0.1 R Evans 14/09/2009
//
//====================================================================================
//
//		Locally allocated variables:
//
//=====================================================================================
#include "atoms.hpp"
#include "errors.hpp"
#include "vmpi.hpp"
#include <iostream>

namespace vmpi{

void mpi_init_halo_swap(){
	//====================================================================================
	//
	///												mpi_init_halo_swap
	///
	///									Initiates halo swap for spin data
	///
	///										Version 1.0 R Evans 16/09/2009
	//
	//====================================================================================
	//
	//		Locally allocated variables:
	//
	//====================================================================================

	//using namespace mpi_comms;
#ifdef MPICF
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){
		std::cout << "mpi_init_halo_swap has been called" << "\t";
		std::cout << vmpi::my_rank << std::endl;
	}

	//----------------------------------------------------------
	// Pack spins for sending
	//----------------------------------------------------------

	//std::cout << "rank: " << vmpi::my_rank << " num_halo_swaps: " << vmpi::send_atom_translation_array.size();
	//std::cout << "\t" << "num_boundary_swaps: " << vmpi::recv_atom_translation_array.size() << std::endl;

	for(unsigned int i=0;i<vmpi::send_atom_translation_array.size();i++){
		int atom = vmpi::send_atom_translation_array[i];
		vmpi::send_spin_data_array[3*i+0] = atoms::x_spin_array[atom];
		vmpi::send_spin_data_array[3*i+1] = atoms::y_spin_array[atom];
		vmpi::send_spin_data_array[3*i+2] = atoms::z_spin_array[atom];
		//std::cout << vmpi::my_rank << "\t" << i << "\t";
		//std::cout << atoms::x_spin_array[vmpi::send_atom_translation_array[i]] << "\t";
		//std::cout << atoms::y_spin_array[vmpi::send_atom_translation_array[i]] << "\t";
		//std::cout << atoms::z_spin_array[vmpi::send_atom_translation_array[i]] << std::endl;
	}

	//----------------------------------------------------------
	// Send spin data array and post receives
	//----------------------------------------------------------

	vmpi::requests.resize(0);
	MPI_Request req;

	for (int p=0;p<vmpi::num_processors;p++){
		if(vmpi::send_num_array[p]!=0){
			int num_pts = 3*vmpi::send_num_array[p];
			int si = 3*vmpi::send_start_index_array[p];
			vmpi::requests.push_back(req);
			MPI_Isend(&vmpi::send_spin_data_array[si],num_pts,MPI_DOUBLE,p,48, MPI_COMM_WORLD, &vmpi::requests.back());
		}
		if(vmpi::recv_num_array[p]!=0){
			int num_pts = 3*vmpi::recv_num_array[p];
			int si = 3*vmpi::recv_start_index_array[p];
			vmpi::requests.push_back(req);
			MPI_Irecv(&vmpi::recv_spin_data_array[si],num_pts,MPI_DOUBLE,p,48, MPI_COMM_WORLD, &vmpi::requests.back());
		}
	}

   #endif
	//----------------------------------------------------------
	// Return
	//----------------------------------------------------------

	return;

}

void mpi_complete_halo_swap(){
	//====================================================================================
	///
	///												mpi_complete_halo_swap
	///
	///									Completes halo swap for spin data
	///
	///										Version 1.0 R Evans 16/09/2009
	//
	//====================================================================================

   #ifdef MPICF

	// check calling of routine if error checking is activated
	if(err::check==true){
		std::cout << "mpi_complete_halo_swap has been called" << "\t";
		std::cout << vmpi::my_rank << std::endl;
	}

	// Swap timers compute -> wait
	vmpi::TotalComputeTime+=vmpi::SwapTimer(vmpi::ComputeTime, vmpi::WaitTime);

	// Wait for all comms to complete
	vmpi::stati.resize(vmpi::requests.size());
	MPI_Waitall(vmpi::requests.size(),&vmpi::requests[0],&vmpi::stati[0]);

	// Swap timers wait -> compute
	vmpi::TotalWaitTime+=vmpi::SwapTimer(vmpi::WaitTime, vmpi::ComputeTime);

	// Unpack received spins
	for(unsigned int i=0;i<vmpi::recv_atom_translation_array.size();i++){
		atoms::x_spin_array[vmpi::recv_atom_translation_array[i]] = vmpi::recv_spin_data_array[3*i+0];
		atoms::y_spin_array[vmpi::recv_atom_translation_array[i]] = vmpi::recv_spin_data_array[3*i+1];
		atoms::z_spin_array[vmpi::recv_atom_translation_array[i]] = vmpi::recv_spin_data_array[3*i+2];
		//std::cout << vmpi::my_rank << "\t" << i << "\t" << vmpi::recv_atom_translation_array[i] << "\t";
		//std::cout << atoms::x_spin_array[vmpi::recv_atom_translation_array[i]] << "\t";
		//std::cout << atoms::y_spin_array[vmpi::recv_atom_translation_array[i]] << "\t";
		//std::cout << atoms::z_spin_array[vmpi::recv_atom_translation_array[i]] << std::endl;
	}

   #endif

	return;

}

} // end of namespace vmpi
