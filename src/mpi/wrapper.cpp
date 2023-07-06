//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "errors.hpp"
#include "vmpi.hpp"

// Internal vmpi header

namespace vmpi{

//------------------------------------------------------------------------------
// Wrapper function for MPI barrier
//------------------------------------------------------------------------------
void barrier(){

   // Wait for all processors just in case anyone else times out
   #ifdef MPICF
      MPI_Barrier(MPI_COMM_WORLD);
   #endif

   return;

}

//------------------------------------------------------------------------------
// Wrapper function(s) for MPI reduce (to master) operation
//------------------------------------------------------------------------------
// **note: return value is only sensible on MASTER process**
//------------------------------------------------------------------------------
uint64_t reduce_sum(uint64_t local){

   uint64_t global = 0;

   #ifdef MPICF
      // Perform MPI reduce for MPI code
      MPI_Reduce(&local, &global, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
   #else
      // set global variable equal to local for serial calls
      global = local;
   #endif

   return global;

}

double reduce_sum(double local){

   double global = 0.0;

   #ifdef MPICF
      // Perform MPI reduce for MPI code
      MPI_Reduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   #else
      // set global variable equal to local for serial calls
      global = local;
   #endif

   return global;

}

//------------------------------------------------------------------------------
// Wrapper function(s) for MPI all reduce operation
//------------------------------------------------------------------------------
uint64_t all_reduce_sum(uint64_t local){

   uint64_t global = 0;

   #ifdef MPICF
      // Perform MPI reduce for MPI code
      MPI_Allreduce(&local, &global, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
   #else
      // set global variable equal to local for serial calls
      global = local;
   #endif

   return global;

}

double all_reduce_sum(double local){

   double global = 0.0;

   #ifdef MPICF
      // Perform MPI reduce for MPI code
      MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #else
      // set global variable equal to local for serial calls
      global = local;
   #endif

   return global;

}

void all_reduce_sum(std::vector<double>& array){

   #ifdef MPICF
      // Perform MPI reduce for MPI code
      MPI_Allreduce(MPI_IN_PLACE, &array[0], array.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

}

void all_reduce_sum(std::vector<int>& array){

   #ifdef MPICF
      // Perform MPI reduce for MPI code
      MPI_Allreduce(MPI_IN_PLACE, &array[0], array.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif

}

//--------------------------------------------------------------------------------------
// Function to collate an array on master from distributed components on all processors
//--------------------------------------------------------------------------------------
void collate(std::vector<double>& input, std::vector<double>& output){

#ifdef MPICF

   //--------------------------------------------------------
   // Parallel version
   //--------------------------------------------------------

   std::vector<int> counts(vmpi::num_processors); // number of components received from each processor
   std::vector<int> displacements(vmpi::num_processors); // offset of each array received from each processor

   const int num_local_data = input.size();

   // gather number of data to be received from each processor
   MPI_Gather(&num_local_data, 1, MPI_INT, &counts[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

   // calculate displacements for gatherv and total number of points
   if(vmpi::master){
      int disp = 0; // initial displacement
      for(int p = 0; p < vmpi::num_processors; p++ ){
         displacements[p] = disp;
         disp += counts[p]; // add number of counts to be recieved
      }
      // Check output array is correct size on master to hold collated data, if not then exit disgracefully
      if( output.size() != disp ){
         std::cerr << "Programmer error in vmpi::collate() due to output array being the wrong size." << std::endl;
         std::cerr << "\t Function call requires " << disp << " elements in output buffer which only has " << output.size() << " elements on master process" << std::endl;
         err::vexit();
      }
   }

   // wait here for everyone to allow master to check memory size
   MPI_Barrier(MPI_COMM_WORLD);

   // Now collate data on master process
   MPI_Gatherv(&input[0], input.size(), MPI_DOUBLE, &output[0], &counts[0], &displacements[0], MPI_DOUBLE, vmpi::master_id, MPI_COMM_WORLD);

#else

   //--------------------------------------------------------
   // Serial version
   //--------------------------------------------------------

   // Check input and output are the same size
   if( output.size() != input.size() ){
      std::cerr << "Programmer error in serial vmpi::collate() due to output array being the wrong size." << std::endl;
      std::cerr << "\t Function call requires " << input.size() << " elements in output buffer which only has " << output.size() << " elements on master process" << std::endl;
      err::vexit();
   }

   // Now hard copy array elements
   for(size_t i = 0 ; i < input.size(); i++){
      output[i] = input[i];
   }

#endif

   // all done!
   return;

}

//---------------------------------------------------------------------------------------
// Function to calculate counts and displacements for an MPI gatherv call on input array
//---------------------------------------------------------------------------------------
void counts_and_displacements(std::vector<double>& input, std::vector<double>& output, std::vector<int>& counts, std::vector<int>& displacements){

#ifdef MPICF

   //--------------------------------------------------------
   // Parallel version only
   //--------------------------------------------------------

   // resize counts and displacements array
   if(vmpi::master) counts.resize(vmpi::num_processors); // number of components received from each processor
   if(vmpi::master) displacements.resize(vmpi::num_processors); // offset of each array received from each processor

   const int num_local_data = input.size();

   // gather number of data to be received from each processor
   MPI_Gather(&num_local_data, 1, MPI_INT, &counts[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

   // calculate displacements for gatherv and total number of points
   if(vmpi::master){
      int disp = 0; // initial displacement
      for(int p = 0; p < vmpi::num_processors; p++ ){
         displacements[p] = disp;
         disp += counts[p]; // add number of counts to be recieved
      }
      // Check output array is correct size on master to hold collated data, if not then exit disgracefully
      if( output.size() != disp ){
         std::cerr << "Programmer error in vmpi::collate() due to output array being the wrong size." << std::endl;
         std::cerr << "\t Function call requires " << disp << " elements in output buffer which only has " << output.size() << " elements on master process" << std::endl;
         err::vexit();
      }
   }

   // wait here for everyone to allow master to check memory size
   MPI_Barrier(MPI_COMM_WORLD);

#endif

   return;

}

//-------------------------------------------------------------------------------------------------------------------------------
// Function to collate an array on master from distributed components on all processors with predefined counts and displacements
//-------------------------------------------------------------------------------------------------------------------------------
void fast_collate(std::vector<double>& input, std::vector<double>& output, std::vector<int>& counts, std::vector<int>& displacements){

#ifdef MPICF

   //--------------------------------------------------------
   // Parallel version
   //--------------------------------------------------------

   // Now collate data on master process
   MPI_Gatherv(&input[0], input.size(), MPI_DOUBLE, &output[0], &counts[0], &displacements[0], MPI_DOUBLE, vmpi::master_id, MPI_COMM_WORLD);

#else

   //--------------------------------------------------------
   // Serial version
   //--------------------------------------------------------

   // Check input and output are the same size
   if( output.size() != input.size() ){
      std::cerr << "Programmer error in serial vmpi::collate() due to output array being the wrong size." << std::endl;
      std::cerr << "\t Function call requires " << input.size() << " elements in output buffer which only has " << output.size() << " elements on master process" << std::endl;
      err::vexit();
   }

   // Now hard copy array elements
   for(size_t i = 0 ; i < input.size(); i++){
      output[i] = input[i];
   }

#endif

   return;

}


//MPI_Allreduce(MPI_IN_PLACE,&stats::sublattice_mean_torque_x_array[0],mp::num_materials,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);

}
