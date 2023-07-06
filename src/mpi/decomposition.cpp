//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iomanip>

// Vampire headers
#include "vmpi.hpp"
#include "vio.hpp"

// vmpi headers

// local struct for storing dimension ranges
struct dim_t {

   double minx;
   double miny;
   double minz;

   double maxx;
   double maxy;
   double maxz;

};

namespace vmpi{

   // Forward function declaration
   dim_t decompose(int id, int num_blocks, dim_t dimensions, std::string segment_name, std::string segment_type);

   //------------------------------------------------------------------------------
   // Function to subdivide system into cubic blocks with topology awareness
   //------------------------------------------------------------------------------
   void geometric_decomposition(int num_cpus, double system_size[3]){

      // declare complete system size range
      dim_t system_dimensions;

      system_dimensions.minx = 0.0;
      system_dimensions.miny = 0.0;
      system_dimensions.minz = 0.0;

      system_dimensions.maxx = system_size[0];
      system_dimensions.maxy = system_size[1];
      system_dimensions.maxz = system_size[2];

      // declare local range
      dim_t local_dimensions;

      // If no node topology required, decompose per processor
      if(vmpi::ppn == 1){

         local_dimensions = decompose(vmpi::my_rank, vmpi::num_processors, system_dimensions,"System", "processors");

      }
      // Otherwise decompose by ppn followed by num_cpus
      // (This improves the volume/surface ratio for better intranode communication)
      else{

         // temporary for storing node dimensions
         dim_t node_dimensions;

         // determine number of nodes [ == roundup(vmpi::num_processors/vmpi::ppn) ]
         const unsigned int num_nodes = (vmpi::num_processors - 1)/vmpi::ppn + 1;

         // determine node id
         const unsigned int my_node_id = vmpi::my_rank/vmpi::ppn;

         // determine processors per node accounting for possibility of last node having fewer processors
         int num_processors_on_my_node = vmpi::ppn;
         if(my_node_id == num_nodes - 1) num_processors_on_my_node = vmpi::ppn - (vmpi::ppn*num_nodes)%vmpi::num_processors;

         // determine node rank
         const unsigned int node_rank = vmpi::my_rank%vmpi::ppn;

         // decompose by number of nodes
         node_dimensions = decompose(my_node_id, num_nodes, system_dimensions, "System", "nodes");

         /*if(node_rank == 0){
            std::cerr << vmpi::my_rank << " " << my_node << " " << num_nodes << " " << num_processors_on_my_node << " " <<
            node_dimensions.minx << " " << node_dimensions.miny << " " << node_dimensions.minz << " " <<
            node_dimensions.maxx << " " << node_dimensions.maxy << " " << node_dimensions.maxz << "nd" << std::endl;
         }*/

         // sub-decompose by number of processes
         local_dimensions = decompose(node_rank, num_processors_on_my_node, node_dimensions, "Nodes", "processors");

      }

      // set namespace variables for partial system generation
      vmpi::min_dimensions[0]=local_dimensions.minx;
      vmpi::min_dimensions[1]=local_dimensions.miny;
      vmpi::min_dimensions[2]=local_dimensions.minz;
      vmpi::max_dimensions[0]=local_dimensions.maxx;
      vmpi::max_dimensions[1]=local_dimensions.maxy;
      vmpi::max_dimensions[2]=local_dimensions.maxz;

      //std::cerr << vmpi::my_rank << " " << my_node << " "
      //          << vmpi::min_dimensions[0] << " " << vmpi::min_dimensions[1] << " " << vmpi::min_dimensions[2] << " "
      //          << vmpi::max_dimensions[0] << " " << vmpi::max_dimensions[1] << " " << vmpi::max_dimensions[2] << std::endl;

      //------------------------------------------------------------------------
      // Homogenise minima and maxima to be the same on all processors
      //------------------------------------------------------------------------
      // This fixes a persistent bug where the maximum on one processor can be
      // greater than the minimum on another processor due to rounding errors
      //------------------------------------------------------------------------

      #ifdef MPICF // MPI-only region
      // store vectors for easy transportation
      std::vector<double> md(6,0.0);

      md[0] = vmpi::min_dimensions[0];
      md[1] = vmpi::min_dimensions[1];
      md[2] = vmpi::min_dimensions[2];
      md[3] = vmpi::max_dimensions[0];
      md[4] = vmpi::max_dimensions[1];
      md[5] = vmpi::max_dimensions[2];

      // collate dimensions on root process
      std::vector<double> mdg;
      if(vmpi::my_rank==0) mdg.resize(6*vmpi::num_processors,0.0);

      // gather values from all other processes
      MPI_Gather(&md[0], 6, MPI_DOUBLE, &mdg[0], 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      //------------------------------------------------------------------------
      // homogenise (N_proc **2 operation, my take a while for > 1000 CPUs...)
      //------------------------------------------------------------------------
      if(vmpi::my_rank==0){

         for(int i = 0; i< mdg.size(); i++){

            const double value = mdg[i];

            for(int j = 0; j< mdg.size(); j++){

               // calculate absolute numerical diff in value
               const double diff = fabs(mdg[j] - value);

               // if difference is within tolerance, then homogenise
               if(diff < 1.0e-7) mdg[j] = value;

            }

         }

      }

      // output homogenised values to file with very high precision
      /*
      if(vmpi::my_rank==0){
         std::ofstream ofile("mma.txt");
         for(int i=0; i< vmpi::num_processors; i++){
            ofile << i << "\t";
            for(int j=0; j<6; j++) ofile << std::setw(52) << std::setprecision(50) << mdg[6*i+j] << "\t";
            ofile << "\n";
         }
         ofile.close();
      }*/

      // now scatter so everyone has the same data
      MPI_Scatter(&mdg[0], 6, MPI_DOUBLE, &md[0], 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      vmpi::min_dimensions[0] = md[0];
      vmpi::min_dimensions[1] = md[1];
      vmpi::min_dimensions[2] = md[2];
      vmpi::max_dimensions[0] = md[3];
      vmpi::max_dimensions[1] = md[4];
      vmpi::max_dimensions[2] = md[5];

      #endif // End of MPI region

      return;

   }

   //-----------------------------------------------------------------------------------
   // Function to decompose a specified 3D range into n blocks, returning the block
   // dimensions for block #id
   //-----------------------------------------------------------------------------------
   dim_t decompose(int id, int num_blocks, dim_t dimensions, std::string segment_name, std::string segment_type){

      // set local constants
      const int x=num_blocks;
      const double lx = dimensions.maxx - dimensions.minx;
      const double ly = dimensions.maxy - dimensions.miny;
      const double lz = dimensions.maxz - dimensions.minz;

      // local variables
      int nx,ny,nz;  /// Number of blocks in x,y,z
      std::vector<int> factor_array; /// to store the factors of each given num_blocks
      factor_array.reserve(50);
      int counter_factor=0; /// to count the number of factors

      //---------------------------------------------------
      // Determine number of blocks in x,y,z
      //---------------------------------------------------

      // find all the factors of given n_blocks (x)
      for (int i=1;i<x+1;i++){
         if ((x%i)==0){
            factor_array.push_back(i);
            counter_factor++;
         }
      }

      // set the remaining elements of the array as 1 if there are no other factors
      for (size_t i=counter_factor+1;i<factor_array.size();i++){
         factor_array[i]=1;
      }

      double surface_volumn=0.0;
      double compare_sv=100000000.0; /// set a very large number for comparing each surface_volumn to find the minimum

      // determine best decomposition for minimizing surface/volume ratio
      for (int i=0;i<counter_factor;i++){
         for (int j=0;j<counter_factor;j++){
            for (int k=0;k<counter_factor;k++){
               int n1=factor_array[i];
               int n2=factor_array[j];
               int n3=factor_array[k];
               // check for valid solution for x blocks
               if (n1*n2*n3==x){
                  // calculate surface/volume ratio
                  surface_volumn = 2.0*(double(n1)/lx+double(n2)/ly+double(n3)/lz);
                  if (surface_volumn < compare_sv) {
                     compare_sv=surface_volumn;
                     nx=n1;
                     ny=n2;
                     nz=n3;
                  }
               }
            }
         }
      }

      // Output informative message to screen and log file
      if(vmpi::my_rank==0){
         std::cout << segment_name << " decomposed into" << " " << nx << " x " << ny << " x "<< nz << " " << segment_type << " for parallel execution" << std::endl;
         zlog << zTs() << segment_name << " decomposed into" << " " << nx << " x " << ny << " x "<< nz << " "  << segment_type << " for parallel execution" << std::endl;
      }

      //---------------------------------------------------
      // Calculate local mpi_dimensions assuming box
      //---------------------------------------------------
      const int my_rank = id;
      const double dx=lx/double(nx);
      const double dy=ly/double(ny);
      const double dz=lz/double(nz);

      // calculate each rank on x, y, z directions respectively
      const int my_rank_x= int(my_rank%((nx)*(ny))/(ny));
      const int my_rank_y=(my_rank%((nx)*(ny)))%(ny);
      const int my_rank_z= int(my_rank/((nx)*(ny)));

      // determine partial dimensions for block id
      dim_t block_dimensions;

      block_dimensions.minx = double(my_rank_x)*dx + dimensions.minx;
      block_dimensions.miny = double(my_rank_y)*dy + dimensions.miny;
      block_dimensions.minz = double(my_rank_z)*dz + dimensions.minz;

      block_dimensions.maxx = double(my_rank_x)*dx + dx + dimensions.minx;
      block_dimensions.maxy = double(my_rank_y)*dy + dy + dimensions.miny;
      block_dimensions.maxz = double(my_rank_z)*dz + dz + dimensions.minz;

      //std::cerr << vmpi::my_rank << "\t" << id << "\t" << nx << "\t" << ny << "\t" << nz << "\trrrrr" << std::endl;

      return block_dimensions;

   }

} // end of namespace vmpi
