//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
// Per-grain magnetic-property bookkeeping: computes per-grain atom count,
// mean coordinates and summed saturation moment from atoms::grain_array.
//-----------------------------------------------------------------------------

#include "atoms.hpp"
#include "constants.hpp"
#include "grains.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

// grains module headers
#include "internal.hpp"

#include <cmath>
#include <iostream>

namespace grains{

   int num_grains=1; // always assume 1 grain
   bool random_anisotropy = false; // flag to control randomly oriented uniaxial anisotropy

   namespace internal{

      std::vector<int> grain_size_array(0);

      std::vector<double> x_coord_array(0);
      std::vector<double> y_coord_array(0);
      std::vector<double> z_coord_array(0);

      std::vector<double> sat_mag_array(0);

      std::vector<double> orientation_array(0); // per-grain Wulff facet orientation, radians

      // RAW (pre-remix) per-grain orientation, as last registered by
      // set_raw_orientation(); consumed once, by the next
      // remix_grain_numbers() call.
      namespace{
         std::vector<double> raw_orientation_buffer(0);
      }

      void set_raw_orientation(const std::vector<double>& raw_orientation){
         raw_orientation_buffer = raw_orientation;
         return;
      }

      //---------------------------------------------------------------------------
      // Remix grain numbers so no empty grains exist
      //---------------------------------------------------------------------------
      void remix_grain_numbers(){

         #ifdef MPICF
            const unsigned int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
         #else
            const unsigned int num_local_atoms = atoms::num_atoms;
         #endif

         std::vector<int> new_grain_numbers(grains::num_grains,0); // array to store new grain number for each grain
         std::vector<int> old_grain_numbers(0); // array to store old grain number for each grain

         //---------------------------------------------------------------------------
         // determine which grains contain atoms
         //---------------------------------------------------------------------------
         std::vector<int> atoms_per_grain(grains::num_grains, 0); // store total atoms per grain on all processors

         for(unsigned int atom = 0; atom < num_local_atoms; atom++){

            // get grain ID
            const int grain = atoms::grain_array[atom];

            // increment grain atom counter
            atoms_per_grain[grain]++;

         }

         #ifdef MPICF
            // add up atoms per grain on all processors
            MPI_Allreduce(MPI_IN_PLACE, atoms_per_grain.data(), grains::num_grains, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         #endif

         // loop over all grains to find unique grain numbers
         for(int grain = 0; grain < grains::num_grains; grain++){

            // if grain has more than one atom, then add to list of grains to keep
            if(atoms_per_grain[grain] > 0){
               old_grain_numbers.push_back(grain);
            }
         }

         // now generate new numbers based on index of old_grain_numbers array
         for(unsigned int grain = 0; grain < old_grain_numbers.size(); grain++){

            const int old_grain_number = old_grain_numbers[grain];
            const int new_grain_number = grain;

            // do a bounds check for sanity
            if( old_grain_number < grains::num_grains){
               new_grain_numbers[old_grain_number] = new_grain_number;
            }
            else{
               std::cerr << "Programmer error! grain index of " << old_grain_number << " is greater than total number of grains " << grains::num_grains << std::endl;
               err::vexit();
            }

         }

         // Now lets renumber the grains
         for(unsigned int atom=0;atom< num_local_atoms;atom++){
            const int old_grain_number = atoms::grain_array[atom];
            const int new_grain_number = new_grain_numbers[old_grain_number];

            // set atom grain ID to new grain number
            atoms::grain_array[atom] = new_grain_number;
         }

         // Remap per-grain orientation through the same old->new compaction
         // just built above. A size mismatch means no granular-film path
         // called set_raw_orientation() for this system, so
         // orientation_array simply defaults to all-zero.
         if(raw_orientation_buffer.size() == static_cast<size_t>(grains::num_grains)){
            orientation_array.assign(old_grain_numbers.size(), 0.0);
            for(size_t grain=0; grain<old_grain_numbers.size(); grain++){
               orientation_array[grain] = raw_orientation_buffer[old_grain_numbers[grain]];
            }
         }
         else{
            orientation_array.assign(old_grain_numbers.size(), 0.0);
         }

         // Finally reset number of grains
         grains::num_grains = old_grain_numbers.size();

         // print helpful message to user about the number of grains generated (if more than 1)
         if(grains::num_grains > 1){
            std::cout << "Generated " << grains::num_grains << " grains" << std::endl;
            zlog << zTs() << "Generated " << grains::num_grains << " grains" << std::endl;
         }

         return;

      }

   } // end of internal namespace

   int set_properties(){
      //========================================================================================================
      ///		 				Function to calculate number of atoms in each grain
      //
      ///														Version 1.0
      //
      ///												R F Evans 15/07/2009
      //
      //========================================================================================================

      //----------------------------------------------------------
      // check calling of routine if error checking is activated
      //----------------------------------------------------------
      if(err::check==true){std::cout << "grains::set_properties has been called" << std::endl;}
      #ifdef MPICF
         const unsigned int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
      #else
         const unsigned int num_local_atoms = atoms::num_atoms;
      #endif

      internal::remix_grain_numbers();

      //---------------------------------------------------------------------------
      // resize and initialise grain arrays
      //---------------------------------------------------------------------------
      internal::grain_size_array.resize(grains::num_grains,0);
      internal::x_coord_array.resize(grains::num_grains,0.0);
      internal::y_coord_array.resize(grains::num_grains,0.0);
      internal::z_coord_array.resize(grains::num_grains,0.0);
      internal::sat_mag_array.resize(grains::num_grains,0.0);

      // loop over atoms to determine grain properties
      for(unsigned int atom=0;atom < num_local_atoms;atom++){

         const int grain = atoms::grain_array[atom]; // grain number
         const int mat   = atoms::type_array[atom];  // material

         // check grain is within allowable bounds
         if( (grain >= 0) && ( grain < grains::num_grains) ){
            internal::grain_size_array[grain]+=1;
            internal::x_coord_array[grain] += atoms::x_coord_array[atom];
            internal::y_coord_array[grain] += atoms::y_coord_array[atom];
            internal::z_coord_array[grain] += atoms::z_coord_array[atom];
            internal::sat_mag_array[grain] += mp::material[mat].mu_s_SI;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - atom " << atom << " belongs to grain " << grain << " which is greater than maximum number of grains ";
            std::cerr << grains::num_grains << ". Exiting" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }

      // Reduce grain properties on all CPUs
      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, &internal::grain_size_array[0],grains::num_grains, MPI_INT,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &internal::x_coord_array[0],grains::num_grains, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &internal::y_coord_array[0],grains::num_grains, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &internal::z_coord_array[0],grains::num_grains, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &internal::sat_mag_array[0],grains::num_grains, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      #endif

      // Calculate mean grains coordinates
      for(int grain=0;grain<grains::num_grains;grain++){
         // check for grains with zero atoms
         if(internal::grain_size_array[grain]==0){
            // grain has no constituent atoms
         }
         else{
            internal::x_coord_array[grain]/=double(internal::grain_size_array[grain]);
            internal::y_coord_array[grain]/=double(internal::grain_size_array[grain]);
            internal::z_coord_array[grain]/=double(internal::grain_size_array[grain]);
         }
      }

      //--------------------------------------------------
      // output grain coordinates to disk on root process
      //--------------------------------------------------
      if( vmpi::my_rank == 0 && grains::num_grains > 1){

         std::ofstream file4;
         file4.open("grain-coordinates.txt");
         file4 << "#-----------------------------------------------------------------------" << std::endl;
         file4 << "# Grain coordinate file for vampire" << std::endl;
         file4 << "# Number of grains: " << grains::num_grains << std::endl;
         file4 << "# Format: grainID x, y, num_atoms, Ms (muB)" << std::endl;
         file4 << "#-----------------------------------------------------------------------" << std::endl;
         for(unsigned int grain=0;grain < internal::x_coord_array.size(); grain++){
            file4 << grain << "\t" << internal::x_coord_array[grain] << "\t" <<
                                       internal::y_coord_array[grain] << "\t" <<
                                       internal::grain_size_array[grain] << "\t" <<
                                       internal::sat_mag_array[grain]/constants::muB << std::endl;
         }
         file4.close();

      }

      return EXIT_SUCCESS;

   }

} // End of namespace grains
