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
#include "grains.hpp"
#include "vio.hpp"
#include "vmath.hpp"

// Internal create header
#include "internal.hpp"

namespace create{
namespace internal{

//------------------------------------------------------------------------------
// Function to cut many particles from lattice in a cubic configuration
//------------------------------------------------------------------------------
void particle_array(std::vector<cs::catom_t> & catom_array){

   // Set number of particles in x and y directions
   const double repeat_size = cs::particle_scale+cs::particle_spacing;
   int num_x_particle = vmath::iceil(cs::system_dimensions[0]/repeat_size);
   int num_y_particle = vmath::iceil(cs::system_dimensions[1]/repeat_size);

   // Loop to generate cubic lattice points
   int particle_number=0;

   std::vector<double> particle_origin(3,0.0);

   for (int x_particle=0;x_particle < num_x_particle;x_particle++){
   	for (int y_particle=0;y_particle < num_y_particle;y_particle++){

         // Determine particle origin
         particle_origin[0] = double(x_particle)*repeat_size + cs::particle_scale*0.5 + cs::particle_array_offset_x;
         particle_origin[1] = double(y_particle)*repeat_size + cs::particle_scale*0.5 + cs::particle_array_offset_y;
         particle_origin[2] = double(vmath::iround(cs::system_dimensions[2]/(2.0*cs::unit_cell.dimensions[2])))*cs::unit_cell.dimensions[2];

         create::internal::centre_particle_on_atom(particle_origin, catom_array);

         if(cs::particle_creation_parity==1){
         	particle_origin[0]+=cs::unit_cell.dimensions[0]*0.5;
         	particle_origin[1]+=cs::unit_cell.dimensions[1]*0.5;
         	particle_origin[2]+=cs::unit_cell.dimensions[2]*0.5;
         }

         // Check to see if a complete particle fits within the system bounds
         if((particle_origin[0]<=(cs::system_dimensions[0]-cs::particle_scale*0.5)) &&
         (particle_origin[1]<=(cs::system_dimensions[1]-cs::particle_scale*0.5))){

            // Use particle type flags to determine which particle shape to cut
            switch(cs::system_creation_flags[1]){
               case 0: // Bulk
                  create::internal::bulk(catom_array);
                  break;
               case 1: // Cube
                  create::internal::cube(particle_origin,catom_array,particle_number);
                  break;
               case 2: // Cylinder
                  create::internal::cylinder(particle_origin,catom_array,particle_number);
                  break;
               case 3: // Ellipsoid
                  create::internal::ellipsoid(particle_origin,catom_array,particle_number);
                  break;
               case 4: // Sphere
                  create::internal::sphere(particle_origin,catom_array,particle_number);
                  break;
               case 5: // Truncated Octahedron
                  create::internal::truncated_octahedron(particle_origin,catom_array,particle_number);
                  break;
               case 6: // Teardrop
                  create::internal::teardrop(particle_origin,catom_array,particle_number);
                  break;
               case 7: // Faceted particle
                  create::internal::faceted(particle_origin,catom_array,particle_number);
                  break;
               case 8: // Cone
                  create::internal::cone(particle_origin,catom_array,particle_number);
                  break;
               case 9: // Bubble
                  create::internal::bubble(particle_origin,catom_array,particle_number);
                  break;
               case 10: // Ellipse
                  create::internal::ellipse(particle_origin,catom_array,0);
                  break;
               default:
                  std::cout << "Unknown particle type requested for single particle system" << std::endl;
                  err::vexit();
            }
            // Increment Particle Number Counter
            particle_number++;
         }
      }
   }

   grains::num_grains = particle_number;

   // Check for no generated particles and print error message
   if(particle_number==0){
     zlog << zTs() << "Error: no particles generated in particle array." << std::endl;
     zlog << zTs() << "Info: Particle arrays require that at least 1 complete particle fits within the system dimensions." << std::endl;
     zlog << zTs() << "Info: Increase x and y system dimensions to at least one particle-scale." << std::endl;
   }

   // Re-order atoms by particle number
   create::internal::sort_atoms_by_grain(catom_array);

   return;

}

} // end of namespace internal
} // end of namespace create
