//-----------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "create.hpp"
#include "errors.hpp"
#include "grains.hpp"
#include "vio.hpp"

// Internal create header
#include "internal.hpp"

namespace create{

namespace internal{

   //------------------------------------------------------------
   // Function to create an agglomorate of nanoparticles
   //------------------------------------------------------------
   void agglomorate(std::vector<cs::catom_t> & catom_array){

      zlog << zTs() << "Generating agglomorate system of grains..." << std::endl;

      //------------------------------------------------
      // load text file with grain positions and radii
      //------------------------------------------------

      // open input file
      std::string agg_file_name = "agglomorate.txt";
      std::ifstream ifile;
      ifile.open(agg_file_name);

      // check file is open
      if(!ifile.is_open()){
         std::cerr << "Error opening file with agglomorate data in file " << agg_file_name << std::endl;
         zlog << zTs() << "Error opening file with agglomorate data in file " << agg_file_name << std::endl;
         err::vexit();
      }

      // Read number of grains from agglomorate file
      int num_grains = 0;
      ifile >> num_grains;

      // check number of grains is sensible
      if( !( num_grains > 0 ) ){
         std::cerr << "Error with agglomorate file " << agg_file_name << " : number of grain must be at least 1" << std::endl;
         zlog << zTs() << "Error with agglomorate file " << agg_file_name << " : number of grain must be at least 1" << std::endl;
         err::vexit();
      }

      // Read data from agglomorate file
      int ln = 0;
      double gx,gy,gz,gr;
      std::vector<double> x,y,z,r; // vectors to store data
      while(!ifile.eof()){

         // read in whole line
         std::string line;
         getline(ifile,line);

         std::stringstream ss(line);
         // check for empty lines
         if(line != ""){

            // read variables in each column
            ss >> gx >> gy >> gz >> gr;

            // save to vectors
            x.push_back(gx);
            y.push_back(gy);
            z.push_back(gz);
            r.push_back(gr);

            std::cout << ln << "\t" << gx << "\t" << gy << "\t" << gz << std::endl;
            ln++;

         }
      }

      // check actual and expected number of grains is the same
      if( num_grains != x.size() ){
         std::cerr << "Error - number of grains in agglomorate file " << x.size() << " is different from expected number of " << num_grains << ". Exiting" << std::endl;
         zlog << zTs() << "Error - number of grains in agglomorate file " << x.size() << " is different from expected number of " << num_grains << ". Exiting" << std::endl;
         err::vexit();
      }

      // set vector to store particle origin
      std::vector<double> particle_origin( 3, 0.0 );

      // initialise counter to count number of generated particles
      int particle_number = 0;

      // loop over all particles to generate local structures
      for(int p = 0; p < num_grains; p++){

         // Determine particle origin, multiplying by system size
         particle_origin[0] = x[p] * cs::system_dimensions[0];
         particle_origin[1] = y[p] * cs::system_dimensions[1];
         particle_origin[2] = z[p] * cs::system_dimensions[2];

         // Determine particle scale
         double agg_particle_scale = r[p] * cs::particle_scale;

         // Check to see if a complete particle fits within the system bounds to avoid partial particles
         const bool in_px = particle_origin[0] <= (cs::system_dimensions[0] - agg_particle_scale * 0.5);
         const bool in_mx = particle_origin[0] >= agg_particle_scale * 0.5;
         const bool in_py = particle_origin[1] <= (cs::system_dimensions[1] - agg_particle_scale * 0.5);
         const bool in_my = particle_origin[1] >= agg_particle_scale * 0.5;
         const bool in_pz = particle_origin[2] <= (cs::system_dimensions[2] - agg_particle_scale * 0.5);
         const bool in_mz = particle_origin[2] >= agg_particle_scale * 0.5;
         const bool inside = in_px && in_mx && in_py && in_my && in_pz && in_mz;

         // now generate particle of the desired shape if its fully insid the system
         // note: could also potentially define agglomorate particles as
         //       fraction of normal particale size
         if(inside){

            // Use particle type flags to determine which particle shape to cut
            // when particla radius is a function parameter
            /*switch(cs::system_creation_flags[1]){
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
                  std::cout << "Programmer error: unknown particle type requested for agglomorate particle system" << std::endl;
                  err::vexit();
            } // end of case statement*/

            create::internal::sphere_size(particle_origin, catom_array, particle_number, agg_particle_scale);


            // Increment Particle Number Counter
            particle_number++;

         } // end of inside check

      } // end of particle generation loop

      // set number ofgrains to particle number
      grains::num_grains = particle_number;

      // Check for no generated particles and print error message
      if( particle_number == 0 ){
         zlog << zTs() << "Error: no particles generated in agglomorate." << std::endl;
         zlog << zTs() << "Info: Agglomorates require that at least 1 complete particle fits within the system dimensions." << std::endl;
         zlog << zTs() << "Info: Increase x and y system dimensions to at least one particle-scale." << std::endl;
         err::vexit();
      }

      // Re-order atoms by particle number
      create::internal::sort_atoms_by_grain(catom_array);

      zlog << zTs() << "   done!" << std::endl;

      return;
   }

} // end of internal namespace

} // end of create namespace
