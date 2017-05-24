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
                   const std::vector<double>& magnetic_moment_array,
                   const std::vector<int>& material_type_array,
                   const std::vector<int>& height_category_array){

      zlog << zTs() << "Initialising statistics module" << std::endl;

      //--------------------------------------------------------------
      // Set up statistics masks for different data sets
      //--------------------------------------------------------------
      stats::num_atoms = num_atoms;

      // define vector mask
      std::vector<int> mask(stats::num_atoms,0);

      // system magnetization
      if(stats::calculate_system_magnetization){
         stats::system_magnetization.set_mask(1,mask,magnetic_moment_array);
      }

      // material magnetization
      if(stats::calculate_material_magnetization){
         for(int atom=0; atom < stats::num_atoms; ++atom) mask[atom] = material_type_array[atom];
         stats::material_magnetization.set_mask(num_materials,mask,magnetic_moment_array);
      }

      // height magnetization
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
         stats::height_magnetization.set_mask(max_height+1,mask,magnetic_moment_array);
      }

      // material height magnetization
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
         stats::material_height_magnetization.set_mask(num_materials*(max_height+1),mask,magnetic_moment_array);
      }

      // system susceptibility
      if(stats::calculate_system_susceptibility) stats::system_susceptibility.initialize(stats::system_magnetization);
      if(stats::calculate_material_susceptibility) stats::material_susceptibility.initialize(stats::material_magnetization);

      return;
   }
} // end of namespace stats
