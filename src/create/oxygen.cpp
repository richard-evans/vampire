//-----------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <list>

// Vampire headers
#include "create.hpp"
#include "vio.hpp"

// Internal sim header
#include "internal.hpp"

namespace create{
   namespace internal{

      //-----------------------------------------------------------------------------
      //
      // Function to generate a random or partially random alloy
      //
      // (c) R F L Evans 2016. All rights reserved.
      //
      //-----------------------------------------------------------------------------
      void oxygen(std::vector<cs::catom_t> & catom_array){

      	// Material parameters not guaranteed to exist
      	// return here if unused to avoid segmentation fault
      	if(create::internal::mp.size() == 0 || create::internal::calculate_oxygen_termination == false) return;

      	// Print informative message to screen
      	zlog << zTs() << "Calculating oxygen termination for system" << std::endl;

         // Loop over unit cell interactions and determine bulk oxygen coordination
         std::vector<int> oxygen_coordination(cs::unit_cell.atom.size(), 0);
         for(unsigned int i = 0; i < cs::unit_cell.interaction.size(); i++){

   			const int atom_i = cs::unit_cell.interaction[i].i;
   			const int atom_j = cs::unit_cell.interaction[i].j;
            const int i_mat = cs::unit_cell.atom[atom_i].mat;
            const int j_mat = cs::unit_cell.atom[atom_j].mat;

            double rx = cs::unit_cell.atom[atom_j].x - cs::unit_cell.atom[atom_i].x + cs::unit_cell.interaction[i].dx;
            double ry = cs::unit_cell.atom[atom_j].y - cs::unit_cell.atom[atom_i].y + cs::unit_cell.interaction[i].dy;
            double rz = cs::unit_cell.atom[atom_j].z - cs::unit_cell.atom[atom_i].z + cs::unit_cell.interaction[i].dz;
            double rij2 = rx*rx + ry*ry + rz*rz;
            double range = create::internal::mp[i_mat].oxygen_coordination_range;

            // check oxygen is in range and the neighbouring atom is oxygen
            bool j_is_oxygen = create::internal::mp[j_mat].oxygen;
            bool in_range = rij2 <= range*range;

            if(j_is_oxygen && in_range) oxygen_coordination[atom_i] += 1;

         }

         /*for(int idx=0; idx < oxygen_coordination.size(); idx++){
            bool idx_is_not_oxygen = !create::internal::mp[cs::unit_cell.atom[idx].mat].oxygen;
            bool is_oxygen_terminated = create::internal::mp[cs::unit_cell.atom[idx].mat].oxygen_terminated;
            if(idx_is_not_oxygen) std::cout << idx << "\t" << is_oxygen_terminated << "\t" << oxygen_coordination[idx] << std::endl;
         }*/

         // copy atom array and make a neighbour list
         std::vector<cs::catom_t> oatom_array = catom_array;
      	std::vector<std::vector<cs::neighbour_t> > oneighbourlist;

         // Copy atoms for interprocessor communications
      	#ifdef MPICF
      	if(vmpi::mpi_mode==0){
      		MPI::COMM_WORLD.Barrier(); // wait for everyone
      		vmpi::copy_halo_atoms(oatom_array);
      		MPI::COMM_WORLD.Barrier(); // sync after halo atoms copied
      	}
      	else if(vmpi::mpi_mode==1){
      		vmpi::set_replicated_data(oatom_array);
      	}
      	#endif

         std::cerr << "Oxygen termination ";
      	// Create Neighbour list for system
      	cs::create_neighbourlist(oatom_array,oneighbourlist);

         int num_trimmed_atoms = 0;

         // loop over all atoms
         for(unsigned int atom = 0; atom < oneighbourlist.size(); atom++){

            int atom_oxygen_coordination = 0;
            const int atom_uc_id = oatom_array[atom].uc_id;
            const int i_mat = oatom_array[atom].material;
            bool is_oxygen_terminated = create::internal::mp[i_mat].oxygen_terminated;
            bool i_included = oatom_array[atom].include;
            // loop over all interactions for oxygen terminated atoms
            if(i_included && is_oxygen_terminated){

               for(unsigned int nn=0;nn<oneighbourlist[atom].size();nn++){

                  // determine if interaction is M-oxygen and in range
                  const int atom_j = oneighbourlist[atom][nn].nn;
                  bool j_included = oatom_array[atom_j].include;
                  const int j_mat = oatom_array[atom_j].material;

                  double rx = oneighbourlist[atom][nn].vx;
                  double ry = oneighbourlist[atom][nn].vy;
                  double rz = oneighbourlist[atom][nn].vz;
                  double rij2 = rx*rx + ry*ry + rz*rz;
                  double range = cs::unit_cell.dimensions[0]*create::internal::mp[i_mat].oxygen_coordination_range;

                  // check oxygen is in range and the neighbouring atom is oxygen
                  bool j_is_oxygen = create::internal::mp[j_mat].oxygen;
                  bool in_range = rij2 <= range*range;

                  if(j_is_oxygen && in_range && j_included) atom_oxygen_coordination += 1;
               }

               //std::cout << "Oxy " << atom << "\t" << atom_uc_id << "\t" << i_mat << "\t" << is_oxygen_terminated << "\t"
               //          << atom_oxygen_coordination << "\t" << oxygen_coordination[atom_uc_id] << std::endl;

               // now trim non-fully coordinated atoms
               if(atom_oxygen_coordination < oxygen_coordination[atom_uc_id]){
                  catom_array[atom].include = false;
                  num_trimmed_atoms++;
               }

            }

         }

         zlog << zTs() << "Trimmed " << num_trimmed_atoms << " non-fully oxygen coordinated atoms from structure" << std::endl;

      	return;

      }

   } // end of namespace internal
} // end of namespace create
