//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

namespace stats{

   void initialize(const int num_atoms,
                   const int num_materials,
                   const int num_grains,
                   const std::vector<double>& magnetic_moment_array,
                   const std::vector<int>& material_type_array,
                   const std::vector<int>& grain_array,
                   const std::vector<int>& height_category_array,
                   const std::vector<bool>& non_magnetic_materials_array){

      zlog << zTs() << "Initialising statistics module" << std::endl;

      //--------------------------------------------------------------
      // Set up statistics masks for different data sets
      //--------------------------------------------------------------
      stats::num_atoms = num_atoms;

      // define vector mask
      std::vector<int> mask(stats::num_atoms,0);

      //------------------------------------------------------------------------
      // system energy
      //------------------------------------------------------------------------
      if(stats::calculate_system_energy){
         for(int atom=0; atom < stats::num_atoms; ++atom){
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] = 1;
            // all other atoms are included
            else mask[atom] = 0;
         }
         stats::system_energy.set_mask(1+1,mask);
      }

      //------------------------------------------------------------------------
      // material energy
      //------------------------------------------------------------------------
      if(stats::calculate_grain_energy){
         for(int atom=0; atom < stats::num_atoms; ++atom){
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] = num_grains;
            // other atoms assigned to material level masks
            else mask[atom] = grain_array[atom];
         }
         stats::grain_energy.set_mask(num_grains+1,mask);
      }

      //------------------------------------------------------------------------
      // material energy
      //------------------------------------------------------------------------
      if(stats::calculate_material_energy){
         for(int atom=0; atom < stats::num_atoms; ++atom){
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] = num_materials;
            // other atoms assigned to material level masks
            else mask[atom] = material_type_array[atom];
         }
         stats::material_energy.set_mask(num_materials+1,mask);
      }

      //------------------------------------------------------------------------
      // system magnetization
      //------------------------------------------------------------------------
      if(stats::calculate_system_magnetization){
         for(int atom=0; atom < stats::num_atoms; ++atom){
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] = 1;
            // all other atoms are included
            else mask[atom] = 0;
         }
         stats::system_magnetization.set_mask(1+1,mask,magnetic_moment_array);
      }

      //------------------------------------------------------------------------
      // grain magnetization
      //------------------------------------------------------------------------
      if(stats::calculate_grain_magnetization){
         for(int atom=0; atom < stats::num_atoms; ++atom){
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] = num_grains;
            // all other atoms are included
            else mask[atom] = grain_array[atom];
         }
         stats::grain_magnetization.set_mask(num_grains+1, mask, magnetic_moment_array);
      }

      //------------------------------------------------------------------------
      // material magnetization
      //------------------------------------------------------------------------
      if(stats::calculate_material_magnetization){
         for(int atom=0; atom < stats::num_atoms; ++atom){
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] =  num_materials;
            // other atoms assigned to material level masks
            else mask[atom] = material_type_array[atom];
         }
         stats::material_magnetization.set_mask(num_materials+1,mask,magnetic_moment_array);
      }

      //------------------------------------------------------------------------
      // material grain magnetization
      //------------------------------------------------------------------------
      if(stats::calculate_material_grain_magnetization){
         // store as blocks of material magnetisation for each grain [ m1x m1y m1z m1m m2x m2y m2z 2m2 ] [ m1x m1y m1z m1m m2x m2y m2z 2m2 ] ...
         // num masks = num_materials*num_grains

         for(int atom=0; atom < stats::num_atoms; ++atom){
            int grain = grain_array[atom];
            int mat = material_type_array[atom];
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] = (num_grains)*num_materials;
            else mask[atom] = num_materials * grain + mat;
         }

         stats::material_grain_magnetization.set_mask(num_materials*(num_grains)+1,mask,magnetic_moment_array);

      }

      //------------------------------------------------------------------------
      // material grain height magnetization
      //------------------------------------------------------------------------
      if(stats::calculate_material_grain_height_magnetization){
         // store as blocks of material magnetisation for each grain and height
         // g1h1 [ m1x m1y m1z m1m m2x m2y m2z 2m2 ]
         // g1h2 [ m1x m1y m1z m1m m2x m2y m2z 2m2 ] ...
         // num masks = num_materials*num_grains*num_heights

         //-----------------------------------------
         // work out maximum height across all CPUs
         //-----------------------------------------
         int max_height=0;
         for(int atom=0; atom < stats::num_atoms; ++atom){
            int hc = height_category_array[atom];
            if(hc>max_height) max_height=hc;
         }
         // Reduce maximum height on all CPUS
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &max_height, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
         #endif

         // calculate num masks
         const int num_heights = max_height+1;
         int num_masks = num_materials*num_heights*num_grains;

         for(int atom=0; atom < stats::num_atoms; ++atom){
            int grain = grain_array[atom];
            int mat = material_type_array[atom];
            int height = height_category_array[atom];
            int mask_id = height * num_grains * num_materials + grain * num_materials + mat;
            int nm_mask_id = num_masks+mat; // group all non-magnetic materials together
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] = nm_mask_id;
            else mask[atom] = mask_id;
         }

         stats::material_grain_height_magnetization.set_mask(1+num_masks+num_materials,mask,magnetic_moment_array);

      }

      //------------------------------------------------------------------------
      // height magnetization
      //------------------------------------------------------------------------
      if(stats::calculate_height_magnetization){
         int max_height=0;
         for(int atom=0; atom < stats::num_atoms; ++atom){
            mask[atom] = height_category_array[atom];
            if(mask[atom]>max_height) max_height=mask[atom];
         }
         // Reduce maximum height on all CPUS
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &max_height, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
         #endif

         // reassign all non-magnetic atoms to last mask
         for(int atom=0; atom < stats::num_atoms; ++atom){
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] =  max_height+1;
         }

         stats::height_magnetization.set_mask(max_height+2,mask,magnetic_moment_array);
      }

      //------------------------------------------------------------------------
      // material height magnetization
      //------------------------------------------------------------------------
      if(stats::calculate_material_height_magnetization){
         // store as blocks of material magnetisation for each height [ m1x m1y m1z m1m m2x m2y m2z 2m2 ] [ m1x m1y m1z m1m m2x m2y m2z 2m2 ] ...
         // num masks = num_materials*num_heights
         int max_height=0;
         for(int atom=0; atom < stats::num_atoms; ++atom){
            int height = height_category_array[atom];
            int mat = material_type_array[atom];
            mask[atom] = num_materials*height+mat;
            if(height>max_height) max_height=height;
         }
         // Reduce maximum height on all CPUS
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &max_height, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
         #endif

         // reassign all non-magnetic atoms to last mask
         for(int atom=0; atom < stats::num_atoms; ++atom){
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] = num_materials*(max_height+1);
         }

         stats::material_height_magnetization.set_mask(num_materials*(max_height+1)+1,mask,magnetic_moment_array);
      }

      //------------------------------------------------------------------------
      // system torque
      //------------------------------------------------------------------------
      if(stats::calculate_system_torque){
         for(int atom=0; atom < stats::num_atoms; ++atom){
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] = 1;
            // all other atoms are included
            else mask[atom] = 0;
         }
         stats::system_torque.set_mask(1+1,mask,magnetic_moment_array);
      }

      //------------------------------------------------------------------------
      // grain torque
      //------------------------------------------------------------------------
      if(stats::calculate_grain_torque){
         for(int atom=0; atom < stats::num_atoms; ++atom){
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] = num_grains;
            // all other atoms are included
            else mask[atom] = grain_array[atom];
         }
         stats::grain_torque.set_mask(num_grains+1, mask, magnetic_moment_array);
      }

      //------------------------------------------------------------------------
      // material torque
      //------------------------------------------------------------------------
      if(stats::calculate_material_torque){
         for(int atom=0; atom < stats::num_atoms; ++atom){
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] =  num_materials;
            // other atoms assigned to material level masks
            else mask[atom] = material_type_array[atom];
         }
         stats::material_torque.set_mask(num_materials+1,mask,magnetic_moment_array);
      }

      //------------------------------------------------------------------------
      // system specific heat
      //------------------------------------------------------------------------
      if(stats::calculate_system_specific_heat)   stats::system_specific_heat.initialize(stats::system_energy);
      if(stats::calculate_grain_specific_heat)    stats::grain_specific_heat.initialize(stats::grain_energy);
      if(stats::calculate_material_specific_heat) stats::material_specific_heat.initialize(stats::material_energy);

      //------------------------------------------------------------------------
      // standard deviation in time
      //------------------------------------------------------------------------
      if(stats::calculate_material_standard_deviation) stats::material_standard_deviation.initialize(stats::system_magnetization);

      //------------------------------------------------------------------------
      // system susceptibility
      //------------------------------------------------------------------------
      if(stats::calculate_system_susceptibility)   stats::system_susceptibility.initialize(stats::system_magnetization);
      if(stats::calculate_grain_susceptibility)    stats::grain_susceptibility.initialize(stats::grain_magnetization);
      if(stats::calculate_material_susceptibility) stats::material_susceptibility.initialize(stats::material_magnetization);

      //------------------------------------------------------------------------
      // binder cumulant
      //------------------------------------------------------------------------
      if(stats::calculate_system_binder_cumulant) stats::system_binder_cumulant.initialize(stats::system_magnetization);
      if(stats::calculate_material_binder_cumulant) stats::material_binder_cumulant.initialize(stats::material_magnetization);

      //------------------------------------------------------------------------
      // system spin length
      //------------------------------------------------------------------------
      if(stats::calculate_system_spin_length){
         for(int atom=0; atom < stats::num_atoms; ++atom){
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] = 1;
            // all other atoms are included
            else mask[atom] = 0;
         }
         stats::system_spin_length.set_mask(1+1,mask);
      }

      //------------------------------------------------------------------------
      // material spin length
      //------------------------------------------------------------------------
      if(stats::calculate_material_spin_length){
         for(int atom=0; atom < stats::num_atoms; ++atom){
            // ignore non-magnetic atoms in stats calculation by assigning them to last mask
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] =  num_materials;
            // other atoms assigned to material level masks
            else mask[atom] = material_type_array[atom];
         }
         stats::material_spin_length.set_mask(num_materials+1,mask);
      }

      //------------------------------------------------------------------------
      // height spin length
      //------------------------------------------------------------------------
      if(stats::calculate_height_spin_length){
         int max_height=0;
         for(int atom=0; atom < stats::num_atoms; ++atom){
            mask[atom] = height_category_array[atom];
            if(mask[atom]>max_height) max_height=mask[atom];
         }
         // Reduce maximum height on all CPUS
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &max_height, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
         #endif

         // reassign all non-magnetic atoms to last mask
         for(int atom=0; atom < stats::num_atoms; ++atom){
            if(non_magnetic_materials_array[material_type_array[atom]]) mask[atom] =  max_height+1;
         }

         stats::height_spin_length.set_mask(max_height+2,mask);
      }

      return;

   }
} // end of namespace stats
