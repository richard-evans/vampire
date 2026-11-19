//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2017. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "create.hpp"

// Internal create header
#include "internal.hpp"

namespace create{

   namespace internal{

      // comparison function for reverse order sorting: decreasing radius, and materials with equal radii (the default
      // core-shell-size of 1 for all materials) in decreasing material index. This must be a strict weak ordering for
      // std::list::sort (false for equal elements); the index tie-break keeps the order the earlier non-strict
      // comparison gave with both libstdc++ and libc++, so material assignment is unchanged.
      bool compare_radius(core_radius_t first,core_radius_t second){
         if(first.radius != second.radius) return first.radius > second.radius;
         return first.mat > second.mat;
      }

   } // end of internal namespace

   //--------------------------------------------------------------------
   // Public equivalents of the two functions above, for use by other
   // modules (e.g. grains) that need core-shell material ordering but
   // must not reach into create::internal directly.
   //--------------------------------------------------------------------
   bool compare_radius(core_radius_t first, core_radius_t second){
      if(first.radius != second.radius) return first.radius > second.radius;
      return first.mat > second.mat;
   }

   std::list<core_radius_t> sorted_core_shell_materials(){

      std::list<core_radius_t> material_order;
      for(int mat=0; mat<mp::num_materials; mat++){
         core_radius_t tmp;
         tmp.mat = mat;
         tmp.radius = mp::material[mat].core_shell_size;
         material_order.push_back(tmp);
      }
      material_order.sort(compare_radius);

      return material_order;

   }

} // end of create namespace
