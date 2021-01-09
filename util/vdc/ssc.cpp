//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2020. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

// program header
#include "vdc.hpp"

// openmp header
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

namespace vdc{

//------------------------------------------------------------------------------
// Function to initialise ssc data for averaging
//------------------------------------------------------------------------------
void initialise_ssc(){

   const double bin_width = 1.0; // 1 Angstrom
   const double inv_bin_width = 1.0 / bin_width;
   const double max_range = sqrt( vdc::system_size[0] * vdc::system_size[0] + vdc::system_size[1] * vdc::system_size[1] + vdc::system_size[2] * vdc::system_size[2]);
   const int num_bins = 1+int(max_range * inv_bin_width); // bin width of 1 Angstrom

   vdc::ssc_num_bins = num_bins;
   vdc::ssc_bin_width = bin_width;
   vdc::ssc_inv_bin_width = inv_bin_width;

   // set data for storing averaged spin-spin correlations
   ssc_counts.resize(ssc_num_bins);
   ssc_correl.resize(ssc_num_bins);

}

//------------------------------------------------------------------------------
// Function to output ssc-xxxxxxxx.txt file in plaintext format
//------------------------------------------------------------------------------
//
// Writes a single file for each snapshot ID formatted as:
//  rij  < sum (si . sj) >
//
//------------------------------------------------------------------------------
void output_ssc_file(unsigned int spin_file_id){

   // Open Povray Include File
	std::stringstream txt_file_sstr;
	txt_file_sstr << "ssc-";
	txt_file_sstr << std::setfill('0') << std::setw(8) << spin_file_id;
	txt_file_sstr << ".txt";
	std::string txt_file = txt_file_sstr.str();

   //--------------------------------------------------------------------------------------------------------------
   // Spin-spin correlation calculation (radially symmetric)
   //
   // https://en.wikipedia.org/wiki/Correlation_function_(statistical_mechanics)
   //
   // C(r) = < S(r) . S(r+dr) > - < M >
   //
   // Data is binned according intratomic distance dr
   //--------------------------------------------------------------------------------------------------------------

   const double bin_width = vdc::ssc_bin_width;
   const double inv_bin_width = vdc::ssc_inv_bin_width;
   const int num_bins = vdc::ssc_num_bins;

   //--------------------------------------------------
   // calculate total magnetization <M> = <Sxyz> . <Sxyz>
   //--------------------------------------------------
   double total_sx = 0.0;
   double total_sy = 0.0;
   double total_sz = 0.0;
   for(unsigned int atomi = 0; atomi < vdc::num_atoms; atomi++){
      total_sx += spins[3*atomi+0];
      total_sy += spins[3*atomi+1];
      total_sz += spins[3*atomi+2];
   }
   const double total_s = ( total_sx * total_sx +
                            total_sy * total_sy +
                            total_sz * total_sz )  /  ( double(vdc::num_atoms) * double(vdc::num_atoms) );

   // data to store counts and summations of s.s
   std::vector<double> counts(num_bins); // number of counts
   std::vector<double> correl(num_bins); // sum of correlations

   // output informative message to user
   if(vdc::verbose) std::cout << "   Generating spin-spin correlation data " << txt_file << "..." << std::flush;

   //---------------------------------------------------------------------------
   // parallelise calculation for better performance
   //---------------------------------------------------------------------------
   #pragma omp parallel
   {

      // thread private data to store counts and summations of s.s
      std::vector<double> thread_counts(num_bins); // number of counts
      std::vector<double> thread_correl(num_bins); // sum of correlations

      #pragma omp for
      for(unsigned int atomi = 0; atomi < vdc::num_atoms; atomi++){
         for(unsigned int atomj = 0; atomj < vdc::num_atoms; atomj++){

            const double dx = coordinates[3*atomj+0] - coordinates[3*atomi+0];
            const double dy = coordinates[3*atomj+1] - coordinates[3*atomi+1];
            const double dz = coordinates[3*atomj+2] - coordinates[3*atomi+2];

            const double sx = spins[3*atomi+0] * spins[3*atomj+0];
            const double sy = spins[3*atomi+1] * spins[3*atomj+1];
            const double sz = spins[3*atomi+2] * spins[3*atomj+2];

            const double rij = sqrt(dx*dx + dy*dy + dz*dz);
            const int index = int(rij * inv_bin_width);

            thread_correl[index] += sx+sy+sz;
            thread_counts[index] += 1.0;

         }
      }

      // now reduce thread private arrays
      #pragma omp critical
      for(int i=0; i<num_bins; i++) counts[i] += thread_counts[i];

      #pragma omp critical
      for(int i=0; i<num_bins; i++) correl[i] += thread_correl[i];

   } // end of parallel region

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   //--------------------------------------------------------------------------------------------------------------

   // output informative message to user
   if(vdc::verbose) std::cout << "   Writing spin-spin correlation file " << txt_file << "..." << std::flush;

   // open incfile
   std::ofstream txtfile;
   txtfile.open(txt_file.c_str());

   for(int i=0; i<num_bins; i++){
      int count = counts[i];
      double sum = correl[i];
      double corr = count > 0 ? sum/count : 0.0;
      if(count > 0) txtfile << (0.5 + double(i)) * bin_width << "\t" << count << "\t" << sum << "\t" << corr << "\t" << total_s << "\t" << corr - total_s << std::endl;
   }

   // flush data to include file and close
   txtfile << std::flush;
   txtfile.close();

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   //-----------------------------------------
   // save data to average arrays
   //-----------------------------------------
   for(int i=0; i<num_bins; i++) vdc::ssc_counts[i] += counts[i];
   for(int i=0; i<num_bins; i++) vdc::ssc_correl[i] += correl[i];
   vdc::ssc_magnetization += total_s;
   vdc::ssc_snapshots += 1.0;

   return;

}

//------------------------------------------------------------------------------
// Function to output time-averaged ssc.txt file in plaintext format
//------------------------------------------------------------------------------
void output_average_ssc_file(){

   // Open Povray Include File
	std::string txt_file = "ssc-average.txt";

   //--------------------------------------------------------------------------------------------------------------
   // Spin-spin correlation calculation (radially symmetric)
   //--------------------------------------------------------------------------------------------------------------
   const int num_bins = ssc_num_bins;
   const double bin_width = ssc_bin_width;

   // output informative message to user
   if(vdc::verbose) std::cout << "   Writing average spin-spin correlation file " << txt_file << "..." << std::flush;

   // open incfile
   std::ofstream txtfile;
   txtfile.open(txt_file.c_str());

   const double inv_snapshots = 1.0 / vdc::ssc_snapshots;
   const double total_s = vdc::ssc_magnetization * inv_snapshots;

   for(int i=0; i<num_bins; i++){
      int count = vdc::ssc_counts[i] * inv_snapshots;
      double sum = vdc::ssc_correl[i] * inv_snapshots;
      double corr = count > 0 ? sum/count : 0.0;
      // only output values with non-zero counts
      if(count > 0) txtfile << (0.5 + double(i)) * bin_width << "\t" << count << "\t" << sum << "\t" << corr << "\t" << total_s << "\t" << corr - total_s << std::endl;
   }

   // flush data to include file and close
   txtfile << std::flush;
   txtfile.close();

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;

}

}
