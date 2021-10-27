//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <list>
#include <vector>
#include <fstream>

// Vampire headers
#include "atoms.hpp"
#include "create.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"


#ifdef MPICF




namespace create{

   namespace internal{

      // Structure to store key information about virtual atoms in halo region
      struct virtual_particle_t{

         int atom; // real atom number
         double x,y,z; // atomic positions
         int scx, scy, scz; // super cell coordinates
         int material; // material id
         int cpuid; // cpu where real atom resides
         int ucid; // unit cell id of atom

      };

      //-----------------------------------------------------------------
      //
      ///   Function to populate array of virtual particles including
      ///   periodic boundary conditions.
      //
      ///   (c) R F L Evans 26/04/2013
      //
      //-----------------------------------------------------------------
      void atom_needed_by_remote_cpu(int atom, // atom number
         int cpu_rank,
         double x,
         double y,
         double z,
         int scx_offset,
         int scy_offset,
         int scz_offset,
         std::vector<cs::catom_t> & catom_array,
         std::vector<std::vector<virtual_particle_t> >& virtual_particle_array,
         std::vector<double> minimax,
         bool self_interaction
      ){

         // Unpack data to constants
         const double min_x = minimax[0];
         const double min_y = minimax[1];
         const double min_z = minimax[2];
         const double max_x = minimax[3];
         const double max_y = minimax[4];
         const double max_z = minimax[5];
         /*
         std::cout << "cpu = " << cpu_rank << "\t";
         std::cout << "min_x = "  << min_x << "\t";
         std::cout << "min_y = "  << min_y << "\t";
         std::cout << "min_z = "  << min_z << "\t";
         std::cout << "max_x = "  << max_x << "\t";
         std::cout << "max_y = "  << max_y << "\t";
         std::cout << "max_z = "  << max_z << std::endl;
         */
         // temporary virtual particle
         virtual_particle_t temp_vp;

         // Check for atom within area needed by other CPU
         if(cpu_rank!=vmpi::my_rank || self_interaction==true){
            if(( (x >= min_x) && (x <= max_x) ) &&
            ( (y >= min_y) && (y <= max_y) ) &&
            ( (z >= min_z) && (z <= max_z) ) ){

               temp_vp.atom = atom;
               temp_vp.x=x;
               temp_vp.y=y;
               temp_vp.z=z;
               temp_vp.scx=catom_array[atom].scx+scx_offset;
               temp_vp.scy=catom_array[atom].scy+scy_offset;
               temp_vp.scz=catom_array[atom].scz+scz_offset;
               temp_vp.material=catom_array[atom].material;
               temp_vp.cpuid=vmpi::my_rank;
               temp_vp.ucid=catom_array[atom].uc_id;

               // pushback list of virtual particles
               virtual_particle_array[cpu_rank].push_back(temp_vp);

               return;
            }
         }

         return;
      }

      //-------------------------------------------------------------------------------------------------------
      // Generalized routine to copy atoms needed by other processors including periodic boundary conditions
      //-------------------------------------------------------------------------------------------------------
      void copy_halo_atoms(std::vector<cs::catom_t> & catom_array){

         vmpi::barrier();// wait for everyone

         zlog << zTs() << "Copying halo atoms to other processors..." << std::endl;
         std::cout << "Copying halo atoms to other processors..." << std::flush;

         // instantiate timers
         vutil::vtimer_t timer;
         vutil::vtimer_t mtimer;

         // Record initial number of atoms
         //const int num_local_atoms=catom_array.size();

         // start timers
         mtimer.start();
         timer.start();

         zlog << zTs() << "   Determining CPU ranges in x,y,z" << std::endl;

         // Populate atoms with correct cpuid
         for(unsigned int atom=0; atom < catom_array.size(); atom++){
            catom_array[atom].mpi_cpuid = vmpi::my_rank;
         }

         // Array to store all interaction ranges
         std::vector<double> cpu_range_array(6*vmpi::num_processors,0.0); // Linear Memory for MPI comms

         // Determine range+interaction range of all CPU's
         double max_interaction_range=double(cs::unit_cell.interaction_range);

         // Populate local ranges
         cpu_range_array[6*vmpi::my_rank+0]=vmpi::min_dimensions[0] - max_interaction_range*cs::unit_cell.dimensions[0]-0.01;
         cpu_range_array[6*vmpi::my_rank+1]=vmpi::min_dimensions[1] - max_interaction_range*cs::unit_cell.dimensions[1]-0.01;
         cpu_range_array[6*vmpi::my_rank+2]=vmpi::min_dimensions[2] - max_interaction_range*cs::unit_cell.dimensions[2]-0.01;
         cpu_range_array[6*vmpi::my_rank+3]=vmpi::max_dimensions[0] + max_interaction_range*cs::unit_cell.dimensions[0]+0.01;
         cpu_range_array[6*vmpi::my_rank+4]=vmpi::max_dimensions[1] + max_interaction_range*cs::unit_cell.dimensions[1]+0.01;
         cpu_range_array[6*vmpi::my_rank+5]=vmpi::max_dimensions[2] + max_interaction_range*cs::unit_cell.dimensions[2]+0.01;

         // Reduce data on all CPUs
         MPI_Allreduce(MPI_IN_PLACE, &cpu_range_array[0],6*vmpi::num_processors, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);

         // Copy ranges to 2D array
         std::vector<std::vector<double> > cpu_range_array2D(vmpi::num_processors);
         for(int cpu=0; cpu< vmpi::num_processors; cpu++){
            cpu_range_array2D[cpu].resize(6);
            cpu_range_array2D[cpu][0]=cpu_range_array[6*cpu+0];
            cpu_range_array2D[cpu][1]=cpu_range_array[6*cpu+1];
            cpu_range_array2D[cpu][2]=cpu_range_array[6*cpu+2];
            cpu_range_array2D[cpu][3]=cpu_range_array[6*cpu+3];
            cpu_range_array2D[cpu][4]=cpu_range_array[6*cpu+4];
            cpu_range_array2D[cpu][5]=cpu_range_array[6*cpu+5];
         }

         timer.stop();

         zlog << zTs() << "     Completed. Time taken: " << timer.elapsed_time() << " s" << std::endl;

         // Determine number of atoms on local CPU needed by other CPUs
         std::vector<int> num_send_atoms(vmpi::num_processors,0);
         std::vector<int> num_recv_atoms(vmpi::num_processors,0);

         // Declare array of virtual particles
         std::vector<std::vector<virtual_particle_t> >virtual_particle_array;
         virtual_particle_array.resize(vmpi::num_processors);

         // Array of minima, maxima
         std::vector<double> minimax;

         // Real offsets for periodic boundary calculations
         const double dx = cs::system_dimensions[0];
         const double dy = cs::system_dimensions[1];
         const double dz = cs::system_dimensions[2];

         // Supercell offsets for periodic boundary calculations
         const int sx = cs::total_num_unit_cells[0];
         const int sy = cs::total_num_unit_cells[1];
         const int sz = cs::total_num_unit_cells[2];

         // start timer
         timer.start();

         zlog << zTs() << "   Determining atoms needed by remote CPUs" << std::endl;

         for(unsigned int atom=0;atom<catom_array.size();atom++){
            const double x = catom_array[atom].x;
            const double y = catom_array[atom].y;
            const double z = catom_array[atom].z;
            for(int cpu=0;cpu<vmpi::num_processors;cpu++){

               // Copy minima/maxima
               minimax=cpu_range_array2D[cpu];

               // Populate virtual particles
               atom_needed_by_remote_cpu(atom, cpu, x, y, z, 0, 0, 0, catom_array, virtual_particle_array, minimax, false);

               // Move atoms by +/- unit cell dimensions in different directions to calculate periodic boundaries
               if(cs::pbc[0]==true){
                  atom_needed_by_remote_cpu(atom, cpu, x+dx, y, z, sx, 0, 0, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x-dx, y, z,-sx, 0, 0, catom_array, virtual_particle_array, minimax, true);
               }

               if(cs::pbc[1]==true){
                  atom_needed_by_remote_cpu(atom, cpu, x, y+dy, z, 0, sy, 0, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x, y-dy, z, 0,-sy, 0, catom_array, virtual_particle_array, minimax, true);
               }

               if(cs::pbc[2]==true){
                  atom_needed_by_remote_cpu(atom, cpu, x, y, z+dz, 0, 0, sz, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x, y, z-dz, 0, 0, -sz, catom_array, virtual_particle_array, minimax, true);
               }

               if( cs::pbc[0]==true && cs::pbc[1]==true ){
                  atom_needed_by_remote_cpu(atom, cpu, x+dx, y+dy, z, sx, sy, 0, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x-dx, y+dy, z,-sx, sy, 0, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x+dx, y-dy, z, sx,-sy, 0, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x-dx, y-dy, z,-sx,-sy, 0, catom_array, virtual_particle_array, minimax, true);
               }

               if( cs::pbc[0]==true && cs::pbc[2]==true ){
                  atom_needed_by_remote_cpu(atom, cpu, x+dx, y, z+dz, sx, 0, sz, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x-dx, y, z+dz,-sx, 0, sz, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x+dx, y, z-dz, sx, 0,-sz, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x-dx, y, z-dz,-sx, 0,-sz, catom_array, virtual_particle_array, minimax, true);
               }

               if( cs::pbc[1]==true && cs::pbc[2]==true ){
                  atom_needed_by_remote_cpu(atom, cpu, x, y+dy, z+dz, 0, sy, sz, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x, y-dy, z+dz, 0,-sy, sz, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x, y+dy, z-dz, 0, sy,-sz, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x, y-dy, z-dz, 0,-sy,-sz, catom_array, virtual_particle_array, minimax, true);
               }

               if( cs::pbc[0]==true && cs::pbc[1]==true && cs::pbc[2]==true){
                  atom_needed_by_remote_cpu(atom, cpu, x-dx, y-dy, z-dz, -sx, -sy, -sz, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x+dx, y-dy, z-dz, +sx, -sy, -sz, catom_array, virtual_particle_array, minimax, true);

                  atom_needed_by_remote_cpu(atom, cpu, x-dx, y+dy, z-dz, -sx, +sy, -sz, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x+dx, y+dy, z-dz, +sx, +sy, -sz, catom_array, virtual_particle_array, minimax, true);

                  atom_needed_by_remote_cpu(atom, cpu, x-dx, y-dy, z+dz, -sx, -sy, +sz, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x+dx, y-dy, z+dz, +sx, -sy, +sz, catom_array, virtual_particle_array, minimax, true);

                  atom_needed_by_remote_cpu(atom, cpu, x-dx, y+dy, z+dz, -sx, +sy, +sz, catom_array, virtual_particle_array, minimax, true);
                  atom_needed_by_remote_cpu(atom, cpu, x+dx, y+dy, z+dz, +sx, +sy, +sz, catom_array, virtual_particle_array, minimax, true);
               }
            }
         }

         timer.stop();

         zlog << zTs() << "     Completed. Time taken: " << timer.elapsed_time() << " s" << std::endl;

         // Calulate number of virtual particles for each cpu
         for(int cpu=0;cpu<vmpi::num_processors;cpu++) num_send_atoms[cpu]=virtual_particle_array[cpu].size();

         std::vector<MPI_Request> requests(0);
         std::vector<MPI_Status> stati(0);
         MPI_Request req = MPI_REQUEST_NULL;

         zlog << zTs() << "   Sharing number of atoms to be sent/received from each CPU" << std::endl;
         timer.start();

         // Send/receive number of boundary/halo atoms (manual all-to-all - better as proper all to all)
         /*for(int cpu=0;cpu<vmpi::num_processors;cpu++){
            requests.push_back(req);
            MPI_Isend(&num_send_atoms[cpu],1,MPI_INT,cpu,35, MPI_COMM_WORLD, &requests.back());
            requests.push_back(req);
            MPI_Irecv(&num_recv_atoms[cpu],1,MPI_INT,cpu,35, MPI_COMM_WORLD, &requests.back());
         }
         stati.resize(requests.size());
         MPI_Waitall(requests.size(),&requests[0],&stati[0]);*/

         MPI_Alltoall(&num_send_atoms[0], 1, MPI_INT, &num_recv_atoms[0], 1, MPI_INT, MPI_COMM_WORLD);

         timer.stop();

         zlog << zTs() << "     Completed. Time taken: " << timer.elapsed_time() << " s" << std::endl;

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
         std::vector<int> send_mpi_atom_supercell_array(3*num_bdry_atoms,0);
         std::vector<int> recv_mpi_atom_supercell_array(3*num_halo_atoms,0);
         std::vector<int> send_material_array(num_bdry_atoms,0);
         std::vector<int> recv_material_array(num_halo_atoms,0);
         std::vector<int> send_cpuid_array(num_bdry_atoms,0);
         std::vector<int> recv_cpuid_array(num_halo_atoms,0);
         std::vector<int> send_mpi_atom_num_array(num_bdry_atoms,0);
         std::vector<int> recv_mpi_atom_num_array(num_halo_atoms,0);
         std::vector<int> send_mpi_uc_id_array(num_bdry_atoms,0);
         std::vector<int> recv_mpi_uc_id_array(num_halo_atoms,0);

         // Pack up data for sending
         int counter=0;	// array index

         /* for(int cpu=0;cpu<vmpi::num_processors;cpu++){
         if(cpu!=vmpi::my_rank){
         for(unsigned int atom=0;atom<catom_array.size();atom++){
         if(   ((catom_array[atom].x >= cpu_range_array[6*cpu+0]) && (catom_array[atom].x <= cpu_range_array[6*cpu+3])) &&
         ((catom_array[atom].y >= cpu_range_array[6*cpu+1]) && (catom_array[atom].y <= cpu_range_array[6*cpu+4])) &&
         ((catom_array[atom].z >= cpu_range_array[6*cpu+2]) && (catom_array[atom].z <= cpu_range_array[6*cpu+5]))) {
         send_coord_array[3*counter+0] = catom_array[atom].x;
         send_coord_array[3*counter+1] = catom_array[atom].y;
         send_coord_array[3*counter+2] = catom_array[atom].z;
         send_mpi_atom_supercell_array[3*counter+0] = catom_array[atom].scx;
         send_mpi_atom_supercell_array[3*counter+1] = catom_array[atom].scy;
         send_mpi_atom_supercell_array[3*counter+2] = catom_array[atom].scz;
         send_material_array[counter]  = catom_array[atom].material;
         send_cpuid_array[counter]  	= vmpi::my_rank; //catom_array[atom].mpi_cpu;
         send_mpi_atom_num_array[counter] = atom;
         send_mpi_uc_id_array[counter] = catom_array[atom].uc_id;
         counter++;
      }}}*/

      for(int cpu=0;cpu<vmpi::num_processors;cpu++){
         for(int vpidx=0; vpidx<virtual_particle_array[cpu].size(); vpidx++){
            send_mpi_atom_num_array[counter]           = virtual_particle_array[cpu][vpidx].atom;
            send_coord_array[3*counter+0]              = virtual_particle_array[cpu][vpidx].x;
            send_coord_array[3*counter+1]              = virtual_particle_array[cpu][vpidx].y;
            send_coord_array[3*counter+2]              = virtual_particle_array[cpu][vpidx].z;
            send_mpi_atom_supercell_array[3*counter+0] = virtual_particle_array[cpu][vpidx].scx;
            send_mpi_atom_supercell_array[3*counter+1] = virtual_particle_array[cpu][vpidx].scy;
            send_mpi_atom_supercell_array[3*counter+2] = virtual_particle_array[cpu][vpidx].scz;
            send_material_array[counter]               = virtual_particle_array[cpu][vpidx].material;
            send_cpuid_array[counter]                  = virtual_particle_array[cpu][vpidx].cpuid;
            send_mpi_uc_id_array[counter]              = virtual_particle_array[cpu][vpidx].ucid;
            counter++;
         }
      }

      int send_index=0;
      int recv_index=0;

      zlog << zTs() << "   Sending boundary atom data to all relevant CPUs" << std::endl;
      timer.start();

      // Exchange boundary/halo data
      for(int cpu=0;cpu<vmpi::num_processors;cpu++){
         if(num_send_atoms[cpu]>0){
            requests.push_back(req);
            MPI_Isend(&send_coord_array[3*send_index],3*num_send_atoms[cpu],MPI_DOUBLE,cpu,50, MPI_COMM_WORLD, &requests.back());
            requests.push_back(req);
            MPI_Isend(&send_mpi_atom_supercell_array[3*send_index],3*num_send_atoms[cpu],MPI_INT,cpu,54, MPI_COMM_WORLD, &requests.back());
            requests.push_back(req);
            MPI_Isend(&send_material_array[send_index],num_send_atoms[cpu],MPI_INT,cpu,51, MPI_COMM_WORLD, &requests.back());
            requests.push_back(req);
            MPI_Isend(&send_cpuid_array[send_index],num_send_atoms[cpu],MPI_INT,cpu,52, MPI_COMM_WORLD, &requests.back());
            requests.push_back(req);
            MPI_Isend(&send_mpi_atom_num_array[send_index],num_send_atoms[cpu],MPI_INT,cpu,53, MPI_COMM_WORLD, &requests.back());
            requests.push_back(req);
            MPI_Isend(&send_mpi_uc_id_array[send_index],num_send_atoms[cpu],MPI_INT,cpu,55, MPI_COMM_WORLD, &requests.back());
            //std::cout << "Send complete on CPU " << vmpi::my_rank << " to CPU " << cpu << " at index " << send_index  << std::endl;
            send_index+=num_send_atoms[cpu];
         }
         if(num_recv_atoms[cpu]>0){
            requests.push_back(req);
            MPI_Irecv(&recv_coord_array[3*recv_index],3*num_recv_atoms[cpu],MPI_DOUBLE,cpu,50, MPI_COMM_WORLD, &requests.back());
            requests.push_back(req);
            MPI_Irecv(&recv_mpi_atom_supercell_array[3*recv_index],3*num_recv_atoms[cpu],MPI_INT,cpu,54, MPI_COMM_WORLD, &requests.back());
            requests.push_back(req);
            MPI_Irecv(&recv_material_array[recv_index],num_recv_atoms[cpu],MPI_INT,cpu,51, MPI_COMM_WORLD, &requests.back());
            requests.push_back(req);
            MPI_Irecv(&recv_cpuid_array[recv_index],num_recv_atoms[cpu],MPI_INT,cpu,52, MPI_COMM_WORLD, &requests.back());
            requests.push_back(req);
            MPI_Irecv(&recv_mpi_atom_num_array[recv_index],num_recv_atoms[cpu],MPI_INT,cpu,53, MPI_COMM_WORLD, &requests.back());
            requests.push_back(req);
            MPI_Irecv(&recv_mpi_uc_id_array[recv_index],num_recv_atoms[cpu],MPI_INT,cpu,55, MPI_COMM_WORLD, &requests.back());
            //std::cout << "Receive complete on CPU " << vmpi::my_rank << " from CPU " << cpu << " at index " << recv_index << " at address " << &recv_mpi_atom_num_array[recv_index] << std::endl;
            recv_index+=num_recv_atoms[cpu];
         }
      }
      stati.resize(requests.size());
      MPI_Waitall(requests.size(),&requests[0],&stati[0]);

      timer.stop();

      zlog << zTs() << "     Completed. Time taken: " << timer.elapsed_time() << " s" << std::endl;

      // Populate halo atoms with data
      for(int index=0;index<num_halo_atoms;index++){
         int atom = catom_array.size();
         catom_array.push_back(cs::catom_t());
         catom_array[atom].x = recv_coord_array[3*index+0];
         catom_array[atom].y = recv_coord_array[3*index+1];
         catom_array[atom].z = recv_coord_array[3*index+2];
         catom_array[atom].scx = recv_mpi_atom_supercell_array[3*index+0];
         catom_array[atom].scy = recv_mpi_atom_supercell_array[3*index+1];
         catom_array[atom].scz = recv_mpi_atom_supercell_array[3*index+2];
         catom_array[atom].material = recv_material_array[index];
         catom_array[atom].mpi_cpuid = recv_cpuid_array[index];
         catom_array[atom].mpi_type = 2; // mark as halo atom
         catom_array[atom].mpi_atom_number = recv_mpi_atom_num_array[index];
         catom_array[atom].uc_id = recv_mpi_uc_id_array[index];
      }

      // wait for everyone
      vmpi::barrier();

      mtimer.stop();

      zlog << zTs() << "\tdone! Total time taken: " << mtimer.elapsed_time() << std::endl;

      if(vmpi::my_rank == 0){
         terminaltextcolor(GREEN);
         std::cout << "done!" << std::endl;
         terminaltextcolor(WHITE);
      }

      return;

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
         terminaltextcolor(RED);
         std::cerr << "Error! - number of atoms is less than number of CPUs - replicated data parallelisation is not possible!" << std::endl;
         terminaltextcolor(WHITE);
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
   void identify_mpi_boundary_atoms(std::vector<cs::catom_t>& catom_array, neighbours::list_t& cneighbourlist){

         // Find and mark boundary and unneeded halo atoms
         for( unsigned int atom = 0; atom < catom_array.size(); atom++ ){

            // define mpi type of local atom
            const int my_mpi_type = catom_array[atom].mpi_type;

            // loop over all neighbours for atom
            for( unsigned int nn = 0; nn < cneighbourlist.list[atom].size(); nn++ ){

               // identify neighbour atom
               const uint64_t natom = cneighbourlist.list[atom][nn].nn;

               // define nearest neighbour MPI type
               int nn_mpi_type = catom_array[natom].mpi_type;

               // Test for interaction with halo
               if( (my_mpi_type == 0 || my_mpi_type == 1) && (nn_mpi_type == 2) ){
                  // a core atom interacting with the halo -> boundary
                  catom_array[atom].boundary = true;
                  // identify halo atom as interacting
                  catom_array[natom].non_interacting_halo = false;
               }

            }

            // Mark atoms appropriately
            if( catom_array[atom].boundary == true ){
               catom_array[atom].mpi_type=1;
            }

         }

         return;
      }


      //----------------------------------------------------------------------------------
      // Simple function to identify non-interacting halo atoms for deletion
      //----------------------------------------------------------------------------------
      void mark_non_interacting_halo(std::vector<cs::catom_t>& catom_array){

         // Find and mark boundary and unneeded halo atoms
         for( unsigned int atom = 0; atom < catom_array.size(); atom++ ){

            // define mpi type of local atom
            const int my_mpi_type = catom_array[atom].mpi_type;

            if( ( my_mpi_type == 2 ) && ( catom_array[atom].non_interacting_halo == true ) ){
               catom_array[atom].mpi_type = 3;
            }
         }

         return;
      }

      /// Define data type storing atom number and mpi_type
      struct data_t {
         int mpi_type;
         int atom_number;
      };

      /// comparison function
      bool compare(data_t first,data_t second){
         if(first.mpi_type<second.mpi_type) return true;
         else return false;
      }

      //------------------------------------------------------------------------
      // Sort atoms accoriding to order core | boundary | halo
      //------------------------------------------------------------------------
      void sort_atoms_by_mpi_type(
         std::vector<cs::catom_t> & catom_array,
         neighbours::list_t& bilinear,
         neighbours::list_t& biquadratic
      ){

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

         // Also need inverse array of atoms for reconstructing neighbour list
         std::vector<int> inv_mpi_type_vec(catom_array.size());

         // loop over new atom list
         for (unsigned int atom=0;atom<catom_array.size();atom++){
            // store new atom number in array of old atom numbers
            inv_mpi_type_vec[mpi_type_vec[atom].atom_number]=atom;

            if(mpi_type_vec[atom].mpi_type !=3) new_num_atoms++;
            if(mpi_type_vec[atom].mpi_type ==0) vmpi::num_core_atoms++;
            if(mpi_type_vec[atom].mpi_type ==1) vmpi::num_bdry_atoms++;
            if(mpi_type_vec[atom].mpi_type ==2) vmpi::num_halo_atoms++;

         }

         // Print out neighbourlist before sorting
         //for (unsigned int atom=0;atom<catom_array.size();atom++){
         //	std::cout << atom << " MPI_t: " <<  catom_array[atom].mpi_type << "NL: ";
         //	for(int i=0;i<cneighbourlist[atom].size();i++) std::cout << cneighbourlist[atom][i].nn << "\t";
         //	std::cout << " on rank " << vmpi::my_rank << std::endl;
         //}

         //for (unsigned int atom=0;atom<catom_array.size();atom++){
         //if(vmpi::my_rank==1){
         //	std::cout << atom << " "<< inv_mpi_type_vec[atom] << " mpi type " << catom_array[atom].mpi_type << std::endl;
         //}
         //}
         //std::cout << vmpi::num_core_atoms << " " << vmpi::num_core_atoms+vmpi::num_bdry_atoms << " " << vmpi::num_core_atoms+vmpi::num_bdry_atoms + vmpi::num_halo_atoms << std::endl;
         zlog << zTs() << "Number of core  atoms: " << vmpi::num_core_atoms << std::endl;
         zlog << zTs() << "Number of local atoms: " << vmpi::num_core_atoms +vmpi::num_bdry_atoms << std::endl;
         zlog << zTs() << "Number of total atoms: " << vmpi::num_core_atoms +vmpi::num_bdry_atoms + vmpi::num_halo_atoms << std::endl;

         // create temporary catom and cneighbourlist arrays for copying data
         std::vector <cs::catom_t> tmp_catom_array(new_num_atoms);
         std::vector <std::vector <neighbours::neighbour_t> > tmp_bilinear(new_num_atoms);
         std::vector <std::vector <neighbours::neighbour_t> > tmp_biquadratic(new_num_atoms);

         // Populate tmp arrays (assuming all mpi_type=3 atoms are at the end of the array?)
         for (unsigned int atom=0;atom<new_num_atoms;atom++){ // new atom number
            unsigned int old_atom_num = mpi_type_vec[atom].atom_number;
            tmp_catom_array[atom]=catom_array[old_atom_num];
            tmp_catom_array[atom].mpi_old_atom_number=old_atom_num; // Store old atom numbers for translation after sorting
            //---//tmp_cneighbourlist[atom].reserve(cneighbourlist[old_atom_num].size());

            //Copy neighbourlist using new atom numbers
            //if(vmpi::my_rank==0) std::cout << vmpi::my_rank << " old " << old_atom_num << " nn: ";
            /*for(unsigned int nn=0;nn<cneighbourlist[old_atom_num].size();nn++){
               unsigned int old_nn_number = cneighbourlist[old_atom_num][nn].nn;
               unsigned int interaction_id= cneighbourlist[old_atom_num][nn].i;
               unsigned int new_nn_number = inv_mpi_type_vec[old_nn_number];
               cs::neighbour_t temp_nt;
               temp_nt.nn=new_nn_number;
               temp_nt.i=interaction_id;
               // Actual neighbours stay the same so simply copy separation vectors
               temp_nt.vx=cneighbourlist[old_atom_num][nn].vx;
               temp_nt.vy=cneighbourlist[old_atom_num][nn].vy;
               temp_nt.vz=cneighbourlist[old_atom_num][nn].vz;
               // ignore all halo-x interactions but not x-halo
               //if(!((mpi_type_vec[atom].mpi_type==2) && (mpi_type_vec[new_nn_number].mpi_type==2)))
               if(!(mpi_type_vec[atom].mpi_type==2))
               tmp_cneighbourlist[atom].push_back(temp_nt);
               //else std::cerr << "ignoring interaction " << atom << "\t" << new_nn_number << " on " << vmpi::my_rank << std::endl;
               //if(vmpi::my_rank==0) std::cout << cneighbourlist[old_atom_num][nn].nn << "\t";
            }*/

            tmp_bilinear[atom].reserve(bilinear.list[old_atom_num].size());

            // Copy bilinear neighbourlist using new atom numbers
            for(unsigned int nn = 0; nn < bilinear.list[old_atom_num].size(); nn++){

               unsigned int old_nn_number = bilinear.list[old_atom_num][nn].nn;
               unsigned int interaction_id= bilinear.list[old_atom_num][nn].i;
               unsigned int new_nn_number = inv_mpi_type_vec[old_nn_number];

               neighbours::neighbour_t temp_nt;

               temp_nt.nn=new_nn_number;
               temp_nt.i=interaction_id;

               // Actual neighbours stay the same so simply copy separation vectors
               temp_nt.vx=bilinear.list[old_atom_num][nn].vx;
               temp_nt.vy=bilinear.list[old_atom_num][nn].vy;
               temp_nt.vz=bilinear.list[old_atom_num][nn].vz;

               // ignore all halo-x interactions but not x-halo
               if(!(mpi_type_vec[atom].mpi_type==2)) tmp_bilinear[atom].push_back(temp_nt);

            }

            // optionally replace biquadratic list
            if(exchange::biquadratic){
               tmp_biquadratic[atom].reserve(biquadratic.list[old_atom_num].size());

               // Copy bilinear neighbourlist using new atom numbers
               for(unsigned int nn = 0; nn < biquadratic.list[old_atom_num].size();nn++){

                  unsigned int old_nn_number = biquadratic.list[old_atom_num][nn].nn;
                  unsigned int interaction_id= biquadratic.list[old_atom_num][nn].i;
                  unsigned int new_nn_number = inv_mpi_type_vec[old_nn_number];

                  neighbours::neighbour_t temp_nt;

                  temp_nt.nn = new_nn_number;
                  temp_nt.i = interaction_id;

                  // Actual neighbours stay the same so simply copy separation vectors
                  temp_nt.vx = biquadratic.list[old_atom_num][nn].vx;
                  temp_nt.vy = biquadratic.list[old_atom_num][nn].vy;
                  temp_nt.vz = biquadratic.list[old_atom_num][nn].vz;

                  // ignore all halo-x interactions but not x-halo
                  if(!(mpi_type_vec[atom].mpi_type==2)) tmp_biquadratic[atom].push_back(temp_nt);

               }
            } // end of biquadratic

         }



         // Swap tmp data over old data more efficient and saves memory
         catom_array.swap(tmp_catom_array);
         //----// cneighbourlist.swap(tmp_cneighbourlist); // This actually works(!) - COPIES both pointers and elements of pointers

         bilinear.list.swap(tmp_bilinear); // This actually works(!) - COPIES both pointers and elements of pointers
         biquadratic.list.swap(tmp_biquadratic);

         // Print out final neighbourlist
         //for (unsigned int atom=0;atom<new_num_atoms;atom++){
         //	std::cout << "Atom: " << atom << " MPI_type: " << catom_array[atom].mpi_type;
         //	for(unsigned int nn=0;nn<cneighbourlist[atom].size();nn++){
         //		std::cout << "\t" << cneighbourlist[atom][nn].nn;
         //	}
         //	std::cout << std::endl;
         //}

         return;
      }

      //------------------------------------------------------------------------
      // initialise MPI comms
      //------------------------------------------------------------------------
      void init_mpi_comms(std::vector<cs::catom_t> & catom_array){

         // check calling of routine if error checking is activated
         if(err::check==true){std::cout << "vmpi::init_mpi_comms has been called" << std::endl;}

         zlog << zTs() << "Initialising MPI communications..." << std::endl;
         std::cout << "Initialising MPI communications..." << std::flush;

         // instantiate timers
         vutil::vtimer_t mtimer;

         // Record initial number of atoms
         //const int num_local_atoms=catom_array.size();

         // start timers
         mtimer.start();

         // Initialise array with number of transfers from and to all CPU's
         vmpi::recv_num_array.resize(vmpi::num_processors);
         vmpi::send_num_array.resize(vmpi::num_processors);
         vmpi::recv_start_index_array.resize(vmpi::num_processors);
         vmpi::send_start_index_array.resize(vmpi::num_processors);

         // Calculate number of spins I need from each CPU
         for(unsigned int atom=0;atom<catom_array.size();atom++){
            // Only receive for halo atoms
            if(catom_array[atom].mpi_type==2){
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
            // Only receive for halo atoms
            if(catom_array[atom].mpi_type==2){
               unsigned int p = catom_array[atom].mpi_cpuid;
               unsigned int index = vmpi::recv_start_index_array[p]+recv_counter_array[p];
               vmpi::recv_atom_translation_array[index]=atom;
               recv_counter_array[p]++;
            }
         }

         // Get number of spins I need to send to each CPU
         std::vector<MPI_Request> requests(0);
         std::vector<MPI_Status> stati(0);
         MPI_Request req = MPI_REQUEST_NULL;

         /*for(int cpu=0;cpu<vmpi::num_processors;cpu++){
            requests.push_back(req);
            MPI_Isend(&vmpi::recv_num_array[cpu],1,MPI_INT,cpu,60, MPI_COMM_WORLD, &requests.back());
            requests.push_back(req);
            MPI_Irecv(&vmpi::send_num_array[cpu],1,MPI_INT,cpu,60, MPI_COMM_WORLD, &requests.back());
         }

         stati.resize(requests.size());
         MPI_Waitall(requests.size(),&requests[0],&stati[0]);*/

         // Get number of spins I need to send to each CPU
         MPI_Alltoall(&vmpi::recv_num_array[0], 1, MPI_INT, &vmpi::send_num_array[0], 1, MPI_INT, MPI_COMM_WORLD);

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

         // post receieves from all processors in advance of data sending
         for(int cpu=0;cpu<vmpi::num_processors;cpu++){
            // check that i have at least one data point to receive
            if(vmpi::send_num_array[cpu] > 0 ){
               int rsi=vmpi::send_start_index_array[cpu];
               requests.push_back(req);
               MPI_Irecv(&vmpi::send_atom_translation_array[rsi],vmpi::send_num_array[cpu],MPI_INT,cpu,61, MPI_COMM_WORLD, &requests.back());
            }
         }

         for(int cpu=0;cpu<vmpi::num_processors;cpu++){
            // Pack remote atom number into 1D array
            //std::vector<int> recv_data(vmpi::recv_num_array[cpu]); // This is very BAD! isend returns immediately, but local array is detroyed = memory mess!
            int si=vmpi::recv_start_index_array[cpu];
            for(int index=0;index<vmpi::recv_num_array[cpu];index++){
               int local_atom_number=vmpi::recv_atom_translation_array[si+index];
               int remote_atom_number=catom_array[local_atom_number].mpi_atom_number;
               recv_data[si+index]=remote_atom_number;
            }
            // check that i have at least one data point to send
            if(vmpi::recv_num_array[cpu] > 0 ){
               requests.push_back(req);
               MPI_Isend(&recv_data[si],vmpi::recv_num_array[cpu],MPI_INT,cpu,61, MPI_COMM_WORLD, &requests.back());
            }
            // check that i have at least one data point to receive
            /*if(vmpi::send_num_array[cpu] > 0 ){
               requests.push_back(req);
               MPI_Irecv(&vmpi::send_atom_translation_array[rsi],vmpi::send_num_array[cpu],MPI_INT,cpu,61, MPI_COMM_WORLD, &requests.back());
            }*/
         }

         stati.resize(requests.size());
         MPI_Waitall(requests.size(),&requests[0],&stati[0]);

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
               terminaltextcolor(RED);
               std::cerr << "Old atom number out of range! on rank " << vmpi::my_rank << "; Old atom number: " << old_atom_num << " ; New atom number: " << atom << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
            inv_atom_translation_array[old_atom_num]=atom;
         }

         // Loop over all atoms to be sent and translate
         for(int cpu=0;cpu<vmpi::num_processors;cpu++){
            int si=vmpi::send_start_index_array[cpu];
            for(int index=0;index<vmpi::send_num_array[cpu];index++){
               int old_atom_number=vmpi::send_atom_translation_array[si+index];
               int new_atom_number=inv_atom_translation_array[old_atom_number];
               vmpi::send_atom_translation_array[si+index]=new_atom_number;
            }
         }

         // wait for everyone
         vmpi::barrier();

         // stop the timer
         mtimer.stop();

         zlog << zTs() << "\tdone! Total time taken: " << mtimer.elapsed_time() << std::endl;
         if(vmpi::my_rank == 0){
            std::cout << "done!" << std::endl;
         }

         return;

      }

   } // end of internal namespace

} // end of create namespace

#endif
