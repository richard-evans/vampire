//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "errors.hpp"
#include "create.hpp"

// Internal create header
#include "internal.hpp"

namespace create{
namespace internal{

//------------------------------------------------------------------------------
//					Function to cut single particle from lattice
//------------------------------------------------------------------------------
void particle(std::vector<cs::catom_t> & catom_array){

   //---------------------------------------------------
   // Set particle origin to atom at centre of lattice
   //---------------------------------------------------
   std::vector<double> particle_origin(3,0.0);

   particle_origin[0] = cs::system_dimensions[0]*0.5;
   particle_origin[1] = cs::system_dimensions[1]*0.5;
   particle_origin[2] = cs::system_dimensions[2]*0.5;

   internal::centre_particle_on_atom(particle_origin, catom_array);

   // check for move in particle origin and that unit cell size < 0.5 system size
   if(cs::particle_creation_parity==1 &&
   	(2.0*cs::unit_cell.dimensions[0]<cs::system_dimensions[0]) &&
   	(2.0*cs::unit_cell.dimensions[1]<cs::system_dimensions[1]) &&
   	(2.0*cs::unit_cell.dimensions[2]<cs::system_dimensions[2])){
   	particle_origin[0]+=cs::unit_cell.dimensions[0]*0.5;
   	particle_origin[1]+=cs::unit_cell.dimensions[1]*0.5;
   	particle_origin[2]+=cs::unit_cell.dimensions[2]*0.5;
   }

   // Use particle type flags to determine which particle shape to cut
   switch(cs::system_creation_flags[1]){
      case 0: // Bulk
         create::internal::bulk(catom_array);
         break;
      case 1: // Cube
         create::internal::cube(particle_origin,catom_array,0);
         break;
      case 2: // Cylinder
         create::internal::cylinder(particle_origin,catom_array,0);
         break;
      case 3: // Ellipsoid
         create::internal::ellipsoid(particle_origin,catom_array,0);
         break;
      case 4: // Sphere
         create::internal::sphere(particle_origin,catom_array,0);
         break;
      case 5: // Truncated Octahedron
         create::internal::truncated_octahedron(particle_origin,catom_array,0);
         break;
      case 6: // Teardrop
         create::internal::teardrop(particle_origin,catom_array,0);
         break;
      case 7: // Faceted particle
         create::internal::faceted(particle_origin,catom_array,0);
         break;
      case 8: // Cone
         create::internal::cone(particle_origin,catom_array,0);
         break;
      case 9: // Bubble
         create::internal::bubble(particle_origin,catom_array,0);
         break;
      case 10: // Ellipse
         create::internal::ellipse(particle_origin,catom_array,0);
         break;
      default:
         std::cout << "Unknown particle type requested for single particle system" << std::endl;
         err::vexit();
   }

	return;
}

} // end of namespace internal
} // end of namespace create
