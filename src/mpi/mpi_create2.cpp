

#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vmpi.hpp"
#include <iostream>
#include <list>
#include <vector>
#include <fstream>

namespace vmpi{
	int mpi_mode;
	int my_rank;
	int num_processors;
	int num_core_atoms;
	int num_bdry_atoms;
	int num_halo_atoms;

	char hostname[20];
	double start_time;
	double end_time;
	double min_dimensions[3]; ///< Minimum coordinates of system on local cpu
	double max_dimensions[3]; ///< Maximum coordinates of system on local cpu
	
	std::vector<int> send_atom_translation_array;
	std::vector<int> send_start_index_array;
	std::vector<int> send_num_array;
	std::vector<double> send_spin_data_array;

	std::vector<int> recv_atom_translation_array;
	std::vector<int> recv_start_index_array;
	std::vector<int> recv_num_array;
	std::vector<double> recv_spin_data_array;
	#ifdef MPICF
	std::vector<MPI::Request> requests(0);
	std::vector<MPI::Status> stati(0);
	#endif
}
	
#ifdef MPICF
	
	
	
	
	namespace vmpi{
	
	int geometric_decomposition(int num_cpus, double system_dimensions[3]){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "vmpi::geometric_decomposition has been called" << std::endl;}

	// set local variables
	int x=num_cpus;
	int nx,ny,nz;	// Number of cpus in x,y,z
	std::vector<int> factor_array; // to store the factors of each given n_cpu
	factor_array.reserve(50);
	int counter_factor=0; // to count the number of factors
	int n1=1; // store n solutions temporary
	int n2=1;
	int n3=1;
	double lx = system_dimensions[0];
	double ly = system_dimensions[1];
	double lz = system_dimensions[2];
	
	double surface_volumn=0.0;
	double compare_sv=10000000.0; // set a very large number for comparing each surface_volumn to find the minimum
	
	// Check for zero cpu's
	if(num_cpus==0){
		std::cerr << "Error - zero cpu's for mpi decomposition, check initialisation of mpi variables" << std::endl;
		err::vexit();
	}
	
	//---------------------------------------------------
	// Determine number of cpu's in x,y,z
	//---------------------------------------------------
	
	// find all the factors of given n_cpu
	for (int i=1;i<x+1;i++){
		if ((x%i)==0){
		  factor_array.push_back(i);
			//cout << i << "\t"<< counter_factor << "\t" << factor_array[counter_factor] << endl;
			counter_factor++;
		}
	}
	
	// set the remaining elements of the array as 1 if there are no other factors
	for (int i=counter_factor+1;i<factor_array.size();i++){
		factor_array[i]=1;
	}

	//cout << counter_factor << endl;
	if (counter_factor==2 && num_cpus>10) {
		std::cerr << num_cpus << "\t" << "cpus cannot be decomposed efficiently, exiting" << std::endl;
		err::vexit();
	}
	else {
		for (int i=0;i<counter_factor;i++){
			for (int j=0;j<counter_factor;j++){
				for (int k=0;k<counter_factor;k++){
					n1=factor_array[i];
					n2=factor_array[j];
					n3=factor_array[k];
					if (n1*n2*n3==x){
						surface_volumn = 2.0*(double(n1)/lx+double(n2)/ly+double(n3)/lz);
						if (surface_volumn < compare_sv) {
							compare_sv=surface_volumn;
							nx=n1;
							ny=n2;
							nz=n3;
						}
					}
				}
			}
		}
		if(vmpi::my_rank==0){
			std::cout << "System decomposed into" << "\t" << nx << " x " << ny << " x "<< nz << " CPUs for parallel execution" << std::endl;
		}
	}
	//---------------------------------------------------
	// Calculate local mpi_dimensions assuming box
	//---------------------------------------------------
	
	int my_rank = vmpi::my_rank;
	double x1,x2,y1,y2,z1,z2; // start and end of each sub-cube
	
	double dx=lx/double(nx); 
	double dy=ly/double(ny); 
	double dz=lz/double(nz);

	// calculate each rank on x, y, z directions respectively
	int my_rank_x= int(my_rank%((nx)*(ny))/(ny));
	int my_rank_y=(my_rank%((nx)*(ny)))%(ny);
	int my_rank_z= int(my_rank/((nx)*(ny)));
		
	//cout <<my_rank_x << "\t" << my_rank_y << "\t" << my_rank_z << endl;

	// x1, x2 stand for start_position and end_position on x axis respectively. Similar as y1,y2,z1,z2.
	x1=double(my_rank_x)*dx;
	x2=double(x1)+dx;
	y1=double(my_rank_y)*dy;
	y2=double(y1)+dy;
	z1=double(my_rank_z)*dz;
	z2=double(z1)+dz;
	
	//std::cout << my_rank << "\t" << x1 << "\t" << x2 << "\t" << y1 << "\t" << y2 << "\t" << z1 << "\t" << z2 << std::endl;
	
	// set namespaced variables
	vmpi::min_dimensions[0]=x1;
	vmpi::min_dimensions[1]=y1;
	vmpi::min_dimensions[2]=z1;
	vmpi::max_dimensions[0]=x2;
	vmpi::max_dimensions[1]=y2;
	vmpi::max_dimensions[2]=z2;
	
	return EXIT_SUCCESS;

}
	
	int crystal_xyz(std::vector<cs::catom_t> & catom_array){
	//====================================================================================
	//
	//												mpi_crystal_xyz
	//
	//					Reduce Coordinate data to head node and output xyz file
	//
	//										Version 1.0 R Evans 10/08/2009
	//
	//====================================================================================
	//
	//		Locally allocated variables: num_atoms_array
	//											  coord_data_array
	//
	//====================================================================================

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){
		std::cout << "vmpi::crystal_xyz has been called " << vmpi::my_rank << std::endl;
	}

	const int num_processors=vmpi::num_processors;
	const int my_rank = vmpi::my_rank;
	const int num_atoms = catom_array.size();
	
	//int* num_atoms_array;
	std::vector<int> num_atoms_array(num_processors);
	//double** coord_data_array;

	// Identify Atoms subject to MPI comms
	/*if(mpi_comms_identify==true){
		for(int atom=0; atom<cs_num_atoms; atom++){
			if(mpi_create_variables::mpi_atom_comm_class_array[atom]==1) atom_type_array[atom]="Li ";
			else if(mpi_create_variables::mpi_atom_comm_class_array[atom]==2) atom_type_array[atom]="H  ";
			else atom_type_array[atom]="Ag ";
		}
	}*/

	//--------------------------------------------------------------------------
	// Wait for root process
	//--------------------------------------------------------------------------
	MPI::COMM_WORLD.Barrier();

	//--------------------------------------------------------------------------
	// Find number of atoms on each node
	//--------------------------------------------------------------------------
	if(my_rank==0){
		for(int p=1;p<num_processors;p++){
			MPI::COMM_WORLD.Recv(&num_atoms_array[p],1,MPI_INT,p,34);
			//std::cout << p << "\t" << num_atoms_array[p] << std::endl;
		}
	}
	else{
		MPI::COMM_WORLD.Send(&num_atoms,1,MPI_INT,0,34);
	}

	//--------------------------------------------------------------------------
	// Send/Receive data from all nodes and output
	//--------------------------------------------------------------------------

	if(my_rank==0){
		int total_num_atoms=num_atoms;
		for(int p=1;p<num_processors;p++){
			total_num_atoms+=num_atoms_array[p];
		}

		std::cout << "Total atoms(all cpu's): " << total_num_atoms << std::endl;
	
		std::ofstream xyz_file;
	  	xyz_file.open ("crystal.xyz");
	  	xyz_file << total_num_atoms + 80<< std::endl;
	  	xyz_file << "" << std::endl;

		// Output axes
		for (int i=0;i<100;i+=5){
			xyz_file << "O\t" << float(i) << "\t" << 0.0 << "\t" << 0.0 << std::endl;
			xyz_file << "O\t" << 0.0 << "\t" << float(i) << "\t" << 0.0 << std::endl;
			xyz_file << "O\t" << material_parameters::system_dimensions[0] << "\t" << material_parameters::system_dimensions[1]-float(i) << "\t" << 0.0 << std::endl;
			xyz_file << "O\t" << material_parameters::system_dimensions[0]-float(i) << "\t" << 	material_parameters::system_dimensions[1] << "\t" << 0.0 << std::endl;
		}
	  	for(int atom=0; atom<num_atoms; atom++){
				//if(mpi_create_variables::mpi_comms_identify==true){
				if(catom_array[atom].mpi_type==1) xyz_file << "Li\t";
				else if(catom_array[atom].mpi_type==2) xyz_file << "Fe\t";
				else if(catom_array[atom].mpi_type==3) xyz_file << "Na\t";
				else xyz_file << mp::material[catom_array[atom].material].element << "\t";
				//}
				//else xyz_file << material_parameters::material[atoms::type_array[atom]].element << "\t";
	  		//xyz_file << mp::material[catom_array[atom].material].element << "\t" << 
	  		xyz_file << catom_array[atom].x << "\t" << 
	  					catom_array[atom].y << "\t" << 
	  					catom_array[atom].z << "\t" << std::endl;
	  		
			//xyz_file << (coord_array[atom][0]+int_mpi_offset[0])*material_parameters::lattice_space_conversion[0] << "\t" << 
  			//				(coord_array[atom][1]+int_mpi_offset[1])*material_parameters::lattice_space_conversion[1] << "\t" << 
  			//				(coord_array[atom][2]+int_mpi_offset[2])*material_parameters::lattice_space_conversion[2] << std::endl;
			//std::cout << atom_type_array[atom] << "\t";
			//std::cout << coord_array[atom][0] << "\t";
			//std::cout << coord_array[atom][1] << "\t";
			//std::cout << coord_array[atom][2] << std::endl;
	  	}

		//string atomt_array[6]={"Ag","H","Co","O","Cl","Li"};

		for(int p=1;p<num_processors;p++){
			std::vector<double> mpi_data_array(3*num_atoms_array[p]);
			//vector<char> mpi_char_array(3*num_atoms_array[p]);
			std::vector<int> mpi_char_array(num_atoms_array[p]);
			std::vector<int> mpi_type_array(num_atoms_array[p]);
			//vector<int> mpi_comms_array(num_atoms_array[p]);
			// Get data from processors
			//std::cout << "Receiving data from rank " << p << std::endl;
			//std::cout << "\t" << "Number of data points expected: " << 3*num_atoms_array[p] << std::endl;
			MPI::COMM_WORLD.Recv(&mpi_data_array[0],3*num_atoms_array[p],MPI_DOUBLE,p,35);
			MPI::COMM_WORLD.Recv(&mpi_char_array[0],num_atoms_array[p],MPI_INT,p,36);
			MPI::COMM_WORLD.Recv(&mpi_type_array[0],num_atoms_array[p],MPI_INT,p,37);
			//MPI::COMM_WORLD.Recv(&mpi_comms_array[0],num_atoms_array[p],MPI_INT,p,37);

			//void MPI::Comm::Recv(void* buf, int count, const MPI::Datatype& datatype,
         //            int source, int tag) const;


			//std::cout << "Receiving data from rank " << p << " completed" << std::endl;
			//for(int i=0;i<100;i++){
			for(int i=0;i<num_atoms_array[p];i++){
				//std::cout << mpi_data_array[3*i+1] << std::endl;
				//xyz_file << atomt_array[p%6] << "\t";
				//xyz_file << mpi_char_array[3*i] << mpi_char_array[3*i+1] << mpi_char_array[3*i+2] << "\t";
				//if(mpi_create_variables::mpi_comms_identify==true){
				//	if(mpi_comms_array[i]==1) xyz_file << "Li\t";
				//	else if(mpi_comms_array[i]==2) xyz_file << "Fe\t";
				//	else xyz_file << material_parameters::material[mpi_char_array[i]].element << "\t";
				//}
				//else 
				if(mpi_type_array[i]==1) xyz_file << "Li\t";
				else if(mpi_type_array[i]==2) xyz_file << "Fe\t";
				else if(mpi_type_array[i]==3) xyz_file << "Na\t";
				else xyz_file << mp::material[mpi_char_array[i]].element << "\t";

				xyz_file << mpi_data_array[3*i];
				xyz_file << "\t" << mpi_data_array[3*i+1]; // + double(p)*mp::system_dimensions[0];
				xyz_file << "\t" << mpi_data_array[3*i+2]; //+material_parameters::system_dimensions[2]*p;
				xyz_file << std::endl;
			}
		}
	}
	else{
		std::vector<double> mpi_data_array(3*num_atoms);
		std::vector<int> mpi_char_array(num_atoms);
		std::vector<int> mpi_type_array(num_atoms);
		//vector<int> mpi_comms_array(num_atoms);
		// Pack data for sending
		for(int i=0;i<num_atoms;i++){
			// Convert atom_type_array to cstr
			//const char * tchar = atom_type_array[i].c_str();
			//char buf[3];
			//for(int c=0;c<3;c++){
			//	buf[c]=tchar[c];
			//}
			mpi_data_array[3*i+0]=catom_array[i].x;
			mpi_data_array[3*i+1]=catom_array[i].y;
			mpi_data_array[3*i+2]=catom_array[i].z;
			
			
			//for(int j=0;j<3;j++){
			//	mpi_data_array[3*i+j]=(coord_array[i][j]+int_mpi_offset[j])*material_parameters::lattice_space_conversion[j];
			//	mpi_char_array[3*i+j]=buf[j];
			//}
			mpi_char_array[i]=catom_array[i].material;
			mpi_type_array[i]=catom_array[i].mpi_type;
			//mpi_comms_array[i]=mpi_create_variables::mpi_atom_comm_class_array[i];
		}
		MPI::COMM_WORLD.Send(&mpi_data_array[0],3*num_atoms,MPI_DOUBLE,0,35);
		MPI::COMM_WORLD.Send(&mpi_char_array[0],num_atoms,MPI_INT,0,36);
		MPI::COMM_WORLD.Send(&mpi_type_array[0],num_atoms,MPI_INT,0,37);
		//MPI::COMM_WORLD.Send(&mpi_comms_array[0],num_atoms,MPI_INT,0,37);
	}

	return EXIT_SUCCESS;
}


int copy_halo_atoms(std::vector<cs::catom_t> & catom_array){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "vmpi::copy_boundary_atoms has been called" << std::endl;}
	
	// Record initial number of atoms
	//const int num_local_atoms=catom_array.size();

	// Populate atoms with correct cpuid
	for(unsigned int atom=0; atom < catom_array.size(); atom++){
		catom_array[atom].mpi_cpuid = vmpi::my_rank;
	}

	// Array to store all interaction ranges
	std::vector<double> cpu_range_array(6*vmpi::num_processors,0.0); // Linear Memory for MPI comms
	
	// Determine range+interaction range of all CPU's
	double max_interaction_range=0.0;
	for(int mat=0;mat<mp::num_materials;mat++){
		if(mp::material[mat].cutoff > max_interaction_range) max_interaction_range = mp::material[mat].cutoff;
	}
	
	// Populate local ranges
	cpu_range_array[6*vmpi::my_rank+0]=vmpi::min_dimensions[0] - max_interaction_range*mp::lattice_constant[0];
	cpu_range_array[6*vmpi::my_rank+1]=vmpi::min_dimensions[1] - max_interaction_range*mp::lattice_constant[1];
	cpu_range_array[6*vmpi::my_rank+2]=vmpi::min_dimensions[2] - max_interaction_range*mp::lattice_constant[2];
	cpu_range_array[6*vmpi::my_rank+3]=vmpi::max_dimensions[0] + max_interaction_range*mp::lattice_constant[0];
	cpu_range_array[6*vmpi::my_rank+4]=vmpi::max_dimensions[1] + max_interaction_range*mp::lattice_constant[1];
	cpu_range_array[6*vmpi::my_rank+5]=vmpi::max_dimensions[2] + max_interaction_range*mp::lattice_constant[2];
	
	// Reduce data on all CPUs
	//MPI::COMM_WORLD.Allreduce(&cpu_range_array[0], &cpu_range_array[0],6*vmpi::num_processors, MPI_DOUBLE,MPI_SUM);
	MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &cpu_range_array[0],6*vmpi::num_processors, MPI_DOUBLE,MPI_SUM);
	//std::cout << vmpi::my_rank << "\t";
	//for(int i=0;i<6*vmpi::num_processors;i++){
	//	std::cout << cpu_range_array[i] << "\t";
	//}
	//std::cout << std::endl;
	
	// Determine number of atoms on local CPU needed by other CPUs
	std::vector<int> num_send_atoms(vmpi::num_processors,0);
	std::vector<int> num_recv_atoms(vmpi::num_processors,0);
	
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		for(int cpu=0;cpu<vmpi::num_processors;cpu++){
			if(cpu!=vmpi::my_rank){
				if(	((catom_array[atom].x >= cpu_range_array[6*cpu+0]) && (catom_array[atom].x <= cpu_range_array[6*cpu+3])) &&
						((catom_array[atom].y >= cpu_range_array[6*cpu+1]) && (catom_array[atom].y <= cpu_range_array[6*cpu+4])) &&
						((catom_array[atom].z >= cpu_range_array[6*cpu+2]) && (catom_array[atom].z <= cpu_range_array[6*cpu+5]))) {
					num_send_atoms[cpu]++;
				}
			}
		}
	}
	
	std::vector<MPI::Request> requests(0);
	std::vector<MPI::Status> stati(0);

	// Send/receive number of boundary/halo atoms
	for(int cpu=0;cpu<vmpi::num_processors;cpu++){
		if(cpu!=vmpi::my_rank){
		    // Bsend requires an explicit buffer, use asynchronous comms instead!
		    //MPI::COMM_WORLD.Bsend(&num_send_atoms[cpu],1,MPI_INT,cpu,35);
		    //MPI::COMM_WORLD.Recv(&num_recv_atoms[cpu],1,MPI_INT,cpu,35);
			requests.push_back(MPI::COMM_WORLD.Isend(&num_send_atoms[cpu],1,MPI_INT,cpu,35));
			requests.push_back(MPI::COMM_WORLD.Irecv(&num_recv_atoms[cpu],1,MPI_INT,cpu,35));
		}
	}

	stati.resize(requests.size());
	MPI::Request::Waitall(requests.size(),&requests[0],&stati[0]);
		
	//std::cout << vmpi::my_rank << "\t";
	//for(int i=0;i<vmpi::num_processors;i++){
	//	std::cout << " CPU: " << i << " S " << num_send_atoms[i] << " R " << num_recv_atoms[i] << "\t";
	//}
	//std::cout << std::endl;
	
	// Calculate total number of boundary and halo atoms on local CPU
	int num_halo_atoms=0;
	int num_bdry_atoms=0;
	for(int cpu=0;cpu<vmpi::num_processors;cpu++){
		num_halo_atoms += num_recv_atoms[cpu];
		num_bdry_atoms += num_send_atoms[cpu];
	}
	
	// Reserve catom array to accomodate halo atoms for neighbourlist calculation
	catom_array.reserve(catom_array.size()+num_halo_atoms);
	
	// Arrays for sending/receiving data
	std::vector<double> send_coord_array(3*num_bdry_atoms,0.0);
	std::vector<double> recv_coord_array(3*num_halo_atoms,0.0);
	std::vector<int> send_material_array(num_bdry_atoms,0);
	std::vector<int> recv_material_array(num_halo_atoms,0);
	std::vector<int> send_cpuid_array(num_bdry_atoms,0);
	std::vector<int> recv_cpuid_array(num_halo_atoms,0);
	std::vector<int> send_mpi_atom_num_array(num_bdry_atoms,0);
	std::vector<int> recv_mpi_atom_num_array(num_halo_atoms,0);

	// Pack up data for sending
	int counter=0;	// array index

	for(int cpu=0;cpu<vmpi::num_processors;cpu++){
		if(cpu!=vmpi::my_rank){
			for(unsigned int atom=0;atom<catom_array.size();atom++){
				if(	((catom_array[atom].x >= cpu_range_array[6*cpu+0]) && (catom_array[atom].x <= cpu_range_array[6*cpu+3])) &&
						((catom_array[atom].y >= cpu_range_array[6*cpu+1]) && (catom_array[atom].y <= cpu_range_array[6*cpu+4])) &&
						((catom_array[atom].z >= cpu_range_array[6*cpu+2]) && (catom_array[atom].z <= cpu_range_array[6*cpu+5]))) {
					send_coord_array[3*counter+0] = catom_array[atom].x;
					send_coord_array[3*counter+1] = catom_array[atom].y;
					send_coord_array[3*counter+2] = catom_array[atom].z;
					send_material_array[counter]  = catom_array[atom].material;
					send_cpuid_array[counter]  	= vmpi::my_rank; //catom_array[atom].mpi_cpu;
					send_mpi_atom_num_array[counter] = atom;
					counter++;
				}
			}
		}
	}
	
	int send_index=0;
	int recv_index=0;
	
	// Wait for everyone to get here
	//MPI::COMM_WORLD.Barrier();
	//std::cout << vmpi::my_rank << " Exchanging Data ";
	//for(int i=0;i<vmpi::num_processors;i++){
	//	std::cout << " CPU: " << i << " S " << num_send_atoms[i] << " R " << num_recv_atoms[i] << "\t";
	//}
	//std::cout << std::endl;
	
	// Exchange boundary/halo data
	
	//std::vector<MPI::Request> requests(0);
	//std::vector<MPI::Status> stati(0);
	
	for(int cpu=0;cpu<vmpi::num_processors;cpu++){
		if(num_send_atoms[cpu]>0){
			//MPI::COMM_WORLD.Ssend(&send_coord_array[3*send_index],3*num_send_atoms[cpu],MPI_DOUBLE,cpu,50);
			//MPI::COMM_WORLD.Ssend(&send_material_array[send_index],num_send_atoms[cpu],MPI_INT,cpu,51);
			//MPI::COMM_WORLD.Ssend(&send_cpuid_array[send_index],num_send_atoms[cpu],MPI_INT,cpu,52);
			//MPI::COMM_WORLD.Ssend(&send_mpi_atom_num_array[send_index],num_send_atoms[cpu],MPI_INT,cpu,53);
			
			requests.push_back(MPI::COMM_WORLD.Isend(&send_coord_array[3*send_index],3*num_send_atoms[cpu],MPI_DOUBLE,cpu,50));
			requests.push_back(MPI::COMM_WORLD.Isend(&send_material_array[send_index],num_send_atoms[cpu],MPI_INT,cpu,51));
			requests.push_back(MPI::COMM_WORLD.Isend(&send_cpuid_array[send_index],num_send_atoms[cpu],MPI_INT,cpu,52));
			requests.push_back(MPI::COMM_WORLD.Isend(&send_mpi_atom_num_array[send_index],num_send_atoms[cpu],MPI_INT,cpu,53));
			//std::cout << "Send complete on CPU " << vmpi::my_rank << " to CPU " << cpu << " at index " << send_index  << std::endl;
			send_index+=num_send_atoms[cpu];
		}
		if(num_recv_atoms[cpu]>0){
			//MPI::COMM_WORLD.Recv(&recv_coord_array[3*recv_index],3*num_recv_atoms[cpu],MPI_DOUBLE,cpu,50);
			//MPI::COMM_WORLD.Recv(&recv_material_array[recv_index],num_recv_atoms[cpu],MPI_INT,cpu,51);
			//MPI::COMM_WORLD.Recv(&recv_cpuid_array[recv_index],num_recv_atoms[cpu],MPI_INT,cpu,52);
			//MPI::COMM_WORLD.Recv(&recv_mpi_atom_num_array[recv_index],num_recv_atoms[cpu],MPI_INT,cpu,53);
			requests.push_back(MPI::COMM_WORLD.Irecv(&recv_coord_array[3*recv_index],3*num_recv_atoms[cpu],MPI_DOUBLE,cpu,50));
			requests.push_back(MPI::COMM_WORLD.Irecv(&recv_material_array[recv_index],num_recv_atoms[cpu],MPI_INT,cpu,51));
			requests.push_back(MPI::COMM_WORLD.Irecv(&recv_cpuid_array[recv_index],num_recv_atoms[cpu],MPI_INT,cpu,52));
			requests.push_back(MPI::COMM_WORLD.Irecv(&recv_mpi_atom_num_array[recv_index],num_recv_atoms[cpu],MPI_INT,cpu,53));
			//std::cout << "Receive complete on CPU " << vmpi::my_rank << " from CPU " << cpu << " at index " << recv_index << " at address " << &recv_mpi_atom_num_array[recv_index] << std::endl;
			recv_index+=num_recv_atoms[cpu];
		}
	}
	stati.resize(requests.size());
	MPI::Request::Waitall(requests.size(),&requests[0],&stati[0]);

	// Populate halo atoms with data
	for(int index=0;index<num_halo_atoms;index++){
		int atom = catom_array.size();
		catom_array.push_back(cs::catom_t());
		catom_array[atom].x = recv_coord_array[3*index+0];
		catom_array[atom].y = recv_coord_array[3*index+1];
		catom_array[atom].z = recv_coord_array[3*index+2];
		catom_array[atom].material = recv_material_array[index];
		catom_array[atom].mpi_cpuid = recv_cpuid_array[index];
		catom_array[atom].mpi_type = 2; // mark as halo atom
		catom_array[atom].mpi_atom_number = recv_mpi_atom_num_array[index];
	}
		
		
	// Wait for everyone to get here
	//MPI::COMM_WORLD.Barrier();
	//MPI::Finalize();
	
	//exit(0);
	
	return EXIT_SUCCESS;
}

/// @brief Set Replicated Data
///
/// @details Sets atom CPU ID for replicated data decomposition 
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    15/03/2011
///
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		15/03/2011
///	Revision:	  ---
///=====================================================================================
///
int set_replicated_data(std::vector<cs::catom_t> & catom_array){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "vmpi::set_replicated_data has been called" << std::endl;}

	// check for num_atoms > num_CPUS
	if(catom_array.size()<vmpi::num_processors){
		std::cerr << "Error! - number of atoms is less than number of CPUs - replicated data parallelisation is not possible!" << std::endl;
		err::vexit();
	}

	// arrays to store atom ranges on each CPU
	std::vector<int> rd_num_atoms(vmpi::num_processors,0);
	std::vector<int> rd_start_atom(vmpi::num_processors,0);
	std::vector<int> rd_end_atom(vmpi::num_processors,0);
	
	// Divide system according to atom numbers, replicated on all CPUs
	for(int p=0;p<vmpi::num_processors;p++){
		rd_num_atoms[p]=catom_array.size()/vmpi::num_processors;
		rd_start_atom[p] = p*rd_num_atoms[p];
		rd_end_atom[p] = (p+1)*rd_num_atoms[p]-1;
	}
	
	// add spare atoms to last CPU
	rd_end_atom[vmpi::num_processors-1]  = catom_array.size()-1;
	rd_num_atoms[vmpi::num_processors-1] = rd_end_atom[vmpi::num_processors-1]-rd_start_atom[vmpi::num_processors-1];
	
	// Populate atoms with CPU id and mpi_type
	for(int p=0;p<vmpi::num_processors;p++){
		if(p==vmpi::my_rank){
			for(int atom=rd_start_atom[p];atom<=rd_end_atom[p];atom++){
				catom_array[atom].mpi_cpuid = p;
				catom_array[atom].mpi_type = 0; // core
				catom_array[atom].mpi_atom_number=atom; // atom numbers are mirrored on all CPUs
			}
		}
		else{
			for(int atom=rd_start_atom[p];atom<=rd_end_atom[p];atom++){
				catom_array[atom].mpi_cpuid = p;
				catom_array[atom].mpi_type = 2; // halo
				catom_array[atom].mpi_atom_number=atom; // atom numbers are mirrored on all CPUs
			}
		}
	}
	
	return EXIT_SUCCESS;
}

int sort_atoms_by_mpi_type(std::vector<cs::catom_t> &,std::vector<std::vector <int> > &);

/// @brief Identify Boundary Atoms
///
/// @details Determines which atoms interact with the halo, assuming all local atoms are 
///          initially designated as core (catom_array[atom].mpi_type = 0)
///          Non-interacting halo atoms (catom_array[atom].mpi_type=3) are marked for 
///          deletion and removed after sorting atoms to core | boundary | halo 
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    15/03/2011
///
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		15/03/2011
///	Revision:	  ---
///=====================================================================================
///
int identify_boundary_atoms(std::vector<cs::catom_t> & catom_array,std::vector<std::vector <int> > & cneighbourlist){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "vmpi::identify_boundary_atoms has been called" << std::endl;}
	
	// Find and mark boundary and unneeded halo atoms 
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		int my_mpi_type=catom_array[atom].mpi_type;
		bool boundary=false;
		bool non_interacting_halo=true;
		for(unsigned int nn=0;nn<cneighbourlist[atom].size();nn++){
			int nn_mpi_type = catom_array[cneighbourlist[atom][nn]].mpi_type;
			// Test for interaction with halo
			if((my_mpi_type==0) && (nn_mpi_type == 2)){
				boundary=true;
			}
			// Test for halo interacting with non-halo
			if((my_mpi_type==2) && ((nn_mpi_type == 0) || (nn_mpi_type == 1))){
				non_interacting_halo=false;
			}
		}
		// Mark atoms appropriately
		if((my_mpi_type==2) && (non_interacting_halo==true)){
			catom_array[atom].mpi_type=3;
		}
		if(boundary==true){
			catom_array[atom].mpi_type=1;
		}
	}
	
	// Sort Arrays by MPI Type
	sort_atoms_by_mpi_type(catom_array,cneighbourlist);
	
	return EXIT_SUCCESS;
}

	// Define data type storing atom number and mpi_type
	struct data_t {
		int mpi_type;
		int atom_number;
	};
	
// comparison function
bool compare(data_t first,data_t second){
	if(first.mpi_type<second.mpi_type) return true;
	else return false;
}
	
int sort_atoms_by_mpi_type(std::vector<cs::catom_t> & catom_array,std::vector<std::vector <int> > & cneighbourlist){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "cs::sort_atoms_by_mpi_type has been called" << std::endl;}	

	// Create list object
	std::list <data_t> mpi_type_list;
	std::list <data_t>::iterator it;
	
	// copy data to list
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		data_t tmp;
		tmp.mpi_type=catom_array[atom].mpi_type;
		tmp.atom_number=atom;
		mpi_type_list.push_back(tmp);
	}

	// sort date in list
	mpi_type_list.sort(compare);
	
	// copy list to vector for ease of access
	std::vector<data_t> mpi_type_vec(mpi_type_list.size());
	copy(mpi_type_list.begin(), mpi_type_list.end(), mpi_type_vec.begin());
	
	// delete temporary list
	mpi_type_list.resize(0);
	
	//if(vmpi::my_rank==1){
	//	for (unsigned int atom=0;atom<catom_array.size();atom++){
	//		std::cout << atom << "\t" << mpi_type_vec[atom].mpi_type << "\t" << mpi_type_vec[atom].atom_number << std::endl;
	//	}
	//}
	
	//Calculate number of atoms excluding non-interacting halo atoms
	unsigned int new_num_atoms=0;
	vmpi::num_core_atoms=0;
	vmpi::num_bdry_atoms=0;
	vmpi::num_halo_atoms=0;

	// Also need inverse array for reconstructing neighbour list
	std::vector<int> inv_mpi_type_vec(catom_array.size());
	
	for (unsigned int atom=0;atom<catom_array.size();atom++){
		inv_mpi_type_vec[mpi_type_vec[atom].atom_number]=atom;

		if(mpi_type_vec[atom].mpi_type !=3) new_num_atoms++;
		if(mpi_type_vec[atom].mpi_type ==0) vmpi::num_core_atoms++;
		if(mpi_type_vec[atom].mpi_type ==1) vmpi::num_bdry_atoms++;
		if(mpi_type_vec[atom].mpi_type ==2) vmpi::num_halo_atoms++;
	}
	
		//for (unsigned int atom=0;atom<catom_array.size();atom++){
		//if(vmpi::my_rank==1){
		//	std::cout << atom << " "<< inv_mpi_type_vec[atom] << " mpi type " << catom_array[atom].mpi_type << std::endl;
		//}
		//}
		std::cout << vmpi::num_core_atoms << " " << vmpi::num_core_atoms+vmpi::num_bdry_atoms << " " << vmpi::num_core_atoms+vmpi::num_bdry_atoms + vmpi::num_halo_atoms << std::endl;
	// create temporary catom and cneighbourlist arrays for copying data
	std::vector <cs::catom_t> tmp_catom_array(new_num_atoms);
	std::vector <std::vector <int> > tmp_cneighbourlist(new_num_atoms);
	//tmp_cneighbourlist.reserve(new_num_atoms);
	
	// Populate tmp arrays (assuming all mpi_type=3 atoms are at the end of the array?)
	for (unsigned int atom=0;atom<new_num_atoms;atom++){ // new atom number
		unsigned int old_atom_num = mpi_type_vec[atom].atom_number;
		tmp_catom_array[atom]=catom_array[old_atom_num];
		tmp_catom_array[atom].mpi_old_atom_number=old_atom_num; // Store old atom numbers for translation after sorting
		//tmp_cneighbourlist.push_back(std::vector<int>());
		tmp_cneighbourlist[atom].reserve(cneighbourlist[old_atom_num].size());
		//Copy neighbourlist using new atom numbers
		//if(vmpi::my_rank==1) std::cout << vmpi::my_rank << " old " << old_atom_num << " nn: ";
		for(unsigned int nn=0;nn<cneighbourlist[old_atom_num].size();nn++){
			unsigned int old_nn_number = cneighbourlist[old_atom_num][nn];
			unsigned int new_nn_number = inv_mpi_type_vec[old_nn_number];
			tmp_cneighbourlist[atom].push_back(new_nn_number);
			//if(vmpi::my_rank==1) std::cout << cneighbourlist[old_atom_num][nn] << "\t";
		}
		//if(vmpi::my_rank==1) std::cout <<std::endl;
		//if(vmpi::my_rank==1) std::cout << vmpi::my_rank << " new " << atom << " nn: ";
		//for(unsigned int nn=0;nn<tmp_cneighbourlist[atom].size();nn++){
			//if(vmpi::my_rank==1) std::cout << tmp_cneighbourlist[atom][nn] << "\t";
		//}
		//if(vmpi::my_rank==1) std::cout <<std::endl;
	}

	// Copy tmp data over old data
	catom_array=tmp_catom_array;
	cneighbourlist=tmp_cneighbourlist; // This actually works(!) - COPIES both pointers and elements of pointers
	
	return EXIT_SUCCESS;
}

int init_mpi_comms(std::vector<cs::catom_t> & catom_array){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "vmpi::init_mpi_comms has been called" << std::endl;}
	
	// Initialise array with number of transfers from and to all CPU's
	vmpi::recv_num_array.resize(vmpi::num_processors);
	vmpi::send_num_array.resize(vmpi::num_processors);
	vmpi::recv_start_index_array.resize(vmpi::num_processors);
	vmpi::send_start_index_array.resize(vmpi::num_processors);
	
	// Calculate number of spins I need from each CPU
	for(unsigned int atom=0;atom<catom_array.size();atom++){
		if(catom_array[atom].mpi_cpuid!=vmpi::my_rank){
			vmpi::recv_num_array[catom_array[atom].mpi_cpuid]++;
		}
	}
	
	// Find total number of halo atoms I need and calculate start index
	int num_halo_swaps=0;
	for(int p=0;p<vmpi::num_processors;p++){
		vmpi::recv_start_index_array[p]=num_halo_swaps;
		num_halo_swaps+=vmpi::recv_num_array[p];
	}
	
	// Resize translation and data arrays
	vmpi::recv_atom_translation_array.resize(num_halo_swaps);
	vmpi::recv_spin_data_array.resize(3*num_halo_swaps);
	
	// Populate recv_translation_array
	std::vector<int> recv_counter_array(vmpi::num_processors);

	for(unsigned int atom=0;atom<catom_array.size();atom++){
		if(catom_array[atom].mpi_cpuid!=vmpi::my_rank){
			unsigned int p = catom_array[atom].mpi_cpuid;
			unsigned int index = vmpi::recv_start_index_array[p]+recv_counter_array[p];
			vmpi::recv_atom_translation_array[index]=atom;
			recv_counter_array[p]++;
		}
	}
	
	//for(int i=0;i<num_halo_swaps;i++){
	//	std::cout << i << "\t" << vmpi::recv_atom_translation_array[i] << std::endl;
	//}

	// Get number of spins I need to send to each CPU
	std::vector<MPI::Request> requests(0);
	std::vector<MPI::Status> stati(0);
	
	for(int cpu=0;cpu<vmpi::num_processors;cpu++){
			requests.push_back(MPI::COMM_WORLD.Isend(&vmpi::recv_num_array[cpu],1,MPI_INT,cpu,60));
			requests.push_back(MPI::COMM_WORLD.Irecv(&vmpi::send_num_array[cpu],1,MPI_INT,cpu,60));
	}

	stati.resize(requests.size());
	MPI::Request::Waitall(requests.size(),&requests[0],&stati[0]);
	
	// Find total number of boundary atoms I need to send and calculate start index

	int num_boundary_swaps=0;
	int num_send_data=0;
	for(int p=0;p<vmpi::num_processors;p++){
		vmpi::send_start_index_array[p]=num_boundary_swaps;
		num_boundary_swaps+=vmpi::send_num_array[p];
		num_send_data+=vmpi::recv_num_array[p];
	}
	
	// Resize translation and data arrays
	vmpi::send_atom_translation_array.resize(num_boundary_swaps);
	vmpi::send_spin_data_array.resize(3*num_boundary_swaps);
	std::vector<int> recv_data(num_send_data);
	// Send and receive atom numbers requested/to be sent
	requests.resize(0);

	for(int cpu=0;cpu<vmpi::num_processors;cpu++){
		// Pack remote atom number into 1D array
	  //std::vector<int> recv_data(vmpi::recv_num_array[cpu]); // This is very BAD! isend returns immediately, but local array is detroyed = memory mess!
		int si=vmpi::recv_start_index_array[cpu];
		for(int index=0;index<vmpi::recv_num_array[cpu];index++){
			int local_atom_number=vmpi::recv_atom_translation_array[si+index];
			int remote_atom_number=catom_array[local_atom_number].mpi_atom_number;
			recv_data[si+index]=remote_atom_number;
		}
		int rsi=vmpi::send_start_index_array[cpu];
		requests.push_back(MPI::COMM_WORLD.Isend(&recv_data[si],vmpi::recv_num_array[cpu],MPI_INT,cpu,61));
		requests.push_back(MPI::COMM_WORLD.Irecv(&vmpi::send_atom_translation_array[rsi],vmpi::send_num_array[cpu],MPI_INT,cpu,61));
	}

	stati.resize(requests.size());
	MPI::Request::Waitall(requests.size(),&requests[0],&stati[0]);
	
	// Translate atoms to be sent from old atom numbers
	// Find highest old atom number
	int highest=catom_array.size();
	for(unsigned int atom=0; atom<catom_array.size();atom++){
          int old_atom_num=catom_array[atom].mpi_old_atom_number;
          if(old_atom_num>highest){
	    highest=old_atom_num;
	  }
	}
	//std::cout << vmpi::my_rank << " highest " << highest << std::endl;
	// Set up atom number translation array
	std::vector <int> inv_atom_translation_array(highest+1);
	for(unsigned int atom=0; atom<catom_array.size();atom++){
	  int old_atom_num=catom_array[atom].mpi_old_atom_number;
	  //std::cout << "Rank: " << vmpi::my_rank << " Old: " << old_atom_num << " New: " << atom << " Highest: " << highest << std::endl;
	  if((old_atom_num>highest) || (old_atom_num < 0)){ // || (old_atom_num>catom_array.size())){
	    std::cout << "Old atom number out of range! on rank " << vmpi::my_rank << "; Old atom number: " << old_atom_num << " ; New atom number: " << atom << std::endl;
	    exit(1);
	  } 
	  inv_atom_translation_array[old_atom_num]=atom;
	}

	// Loop over all atoms to be sent and translate
	for(int cpu=0;cpu<vmpi::num_processors;cpu++){
	  int si=vmpi::send_start_index_array[cpu];
	  for(int index=0;index<vmpi::send_num_array[cpu];index++){
	    int old_atom_number=vmpi::send_atom_translation_array[si+index];
	    int new_atom_number=inv_atom_translation_array[old_atom_number];
	    send_atom_translation_array[si+index]=new_atom_number;
	  }
	}

	return EXIT_SUCCESS;
}
	
	
	}	//	end of namespace vmpi

#endif 

