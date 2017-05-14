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
#include "create.hpp"
#include "material.hpp"

// Internal create header
#include "internal.hpp"

namespace cs{

   //--------------------------------------------------------------------------
   // Function to create a multilayer system replicating a standard stack
   // structure cs::num_multilayers times
   //
   //   |-----------------------|  1.0          |-----------------------| L2 M2
   //   |      Material 2       |               |-----------------------|
   //   |-----------------------|  0.6   x2     |                       | L2 M2
   //   |                       |        -->    |-----------------------| L1 M1
   //   |      Material 1       |               |-----------------------|
   //   |                       |               |                       | L1 M2
   //   |-----------------------|  0.0          |-----------------------|
   //
   // The multilayers code divides the total system height into n multiples
   // of the defined material heights.
   //--------------------------------------------------------------------------
   void generate_multilayers(std::vector<cs::catom_t>& catom_array){

      // Load number of multilayer repeats to temporary constant
      const int num_layers = cs::num_multilayers;

      // Determine fractional system height for each layer
      const double fractional_system_height = cs::system_dimensions[2]/double(num_layers);

      // Unroll min, max and fill for performance
      std::vector<double> mat_min(mp::num_materials);
      std::vector<double> mat_max(mp::num_materials);
      std::vector<bool> mat_fill(mp::num_materials);

      for(int mat=0;mat<mp::num_materials;mat++){
         mat_min[mat]=create::internal::mp[mat].min;
         mat_max[mat]=create::internal::mp[mat].max;
         // alloys generally are not defined by height, and so have max = 0.0
         if(mat_max[mat]<0.0000001) mat_max[mat]=-0.1;
         mat_fill[mat]=mp::material[mat].fill;
      }

      // Assign materials to generated atoms
      for(unsigned int atom=0;atom<catom_array.size();atom++){
         // loop over multilayers
         for(int multi=0; multi < num_layers; ++multi){
            const double multilayer_num = double(multi);
            for(int mat=0;mat<mp::num_materials;mat++){
               // determine mimimum and maximum heigh for this layer
               const double mat_min_z = (mat_min[mat] + multilayer_num)*fractional_system_height;
               const double mat_max_z = (mat_max[mat] + multilayer_num)*fractional_system_height;
               const double cz=catom_array[atom].z;
               const int atom_uc_cat = catom_array[atom].uc_category;
               const int mat_uc_cat = create::internal::mp[mat].unit_cell_category;
               // if in range then allocate to material
               if((cz>=mat_min_z) && (cz<mat_max_z) && (mat_fill[mat]==false) && (atom_uc_cat == mat_uc_cat) ){
                  catom_array[atom].material=mat;
                  catom_array[atom].include=true;
                  // Optionally recategorize heigh magnetization by layer in multilayer
                  if(cs::multilayer_height_category) catom_array[atom].lh_category = multi;
               }
            }
         }
      }

      return;

   }

} // end of namespace cs
