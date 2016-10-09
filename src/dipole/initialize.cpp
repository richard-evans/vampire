//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#include <fstream>
#include <sstream>

// Vampire headers
#include "cells.hpp"
#include "dipole.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

#include "atoms.hpp"

#include <time.h>
#include <fenv.h>
#include <signal.h>

// dipole module headers
#include "internal.hpp"

namespace dipole{

   //----------------------------------------------------------------------------
   // Function to initialize dipole module
   //----------------------------------------------------------------------------
   void initialize(const int cells_num_atoms_in_unit_cell,
                  int cells_num_cells, /// number of macrocells
                  int cells_num_local_cells, /// number of local macrocells
                  const double cells_macro_cell_size,
                  std::vector <int>& cells_local_cell_array,
                  std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                  std::vector < std::vector <int> >& cells_index_atoms_array,

                  const std::vector<double>& cells_volume_array,
                  /*const std::vector<double>& cells_cell_coords_array_x, /// arrays to store cells positions
                  const std::vector<double>& cells_cell_coords_array_y,
                  const std::vector<double>& cells_cell_coords_array_z,*/

                  std::vector<double>& cells_pos_and_mom_array, // array to store positions and moment of cells

                  std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                  std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                  std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,

                  const std::vector<int>& atom_type_array,
                  const std::vector<int>& atom_cell_id_array,

                  const std::vector<double>& atom_coords_x, //atomic coordinates
                  const std::vector<double>& atom_coords_y,
                  const std::vector<double>& atom_coords_z,

                  const int num_atoms
				){

	//-------------------------------------------------------------------------------------
	// Check for dipole calculation enabled, if not do nothing
	//-------------------------------------------------------------------------------------
      if(!dipole::activated) return;

      //if(!internal::enabled) return;
      std::cout << "Initialising demagnetisation field calculation" << std::endl;
		// output informative message
		zlog << zTs() << "Initialising demagnetisation field calculation" << std::endl;

		// check for prior initialisation
		if(dipole::internal::initialised){
      		zlog << zTs() << "Warning:  Demagnetisation field calculation already initialised. Continuing." << std::endl;
      	return;
		}

/*      std::cout << "Initialising demagnetisation field calculation" << std::endl;
      // check for calling of routine
      if(err::check==true){
         terminaltextcolor(RED);
         std::cerr << "demag::set_rij_matrix has been called " << vmpi::my_rank << std::endl;
         terminaltextcolor(WHITE);
      } */


		//-------------------------------------------------------------------------------------
		// Set simulation constants
		//-------------------------------------------------------------------------------------

      //-------------------------------------------------------------------------------------
      // Set const for functions
      //-------------------------------------------------------------------------------------

      dipole::internal::num_atoms                  = num_atoms;
      dipole::internal::atom_type_array            = atom_type_array;
      dipole::internal::atom_cell_id_array         = atom_cell_id_array;

      dipole::internal::cells_num_cells            = cells_num_cells;
      dipole::internal::cells_num_local_cells      = cells_num_local_cells;
      dipole::internal::cells_local_cell_array     = cells_local_cell_array;
      dipole::internal::cells_num_atoms_in_cell    = cells_num_atoms_in_cell;
      dipole::internal::cells_volume_array         = cells_volume_array;

      dipole::internal::cells_pos_and_mom_array    = cells_pos_and_mom_array;

		//-------------------------------------------------------------------------------------
		// Starting calculation of dipolar field
		//-------------------------------------------------------------------------------------

       // timing function
       #ifdef MPICF
          double t1 = MPI_Wtime();
       #else
          time_t t1;
          t1 = time (NULL);
       #endif

         // Check memory requirements and print to screen
         zlog << zTs() << "Fast demagnetisation field calculation has been enabled and requires " << double(dipole::internal::cells_num_cells)*double(dipole::internal::cells_num_local_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;
         std::cout << "Fast demagnetisation field calculation has been enabled and requires " << double(dipole::internal::cells_num_cells)*double(dipole::internal::cells_num_local_cells*6)*8.0/1.0e6 << " MB of RAM" << std::endl;

        // For MPI version, only add local atoms
         #ifdef MPICF
            int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
         #else
            int num_local_atoms = num_atoms;
         #endif

      /*------------------------------------------------------*/
      /* parallelisation section                              */
      /*------------------------------------------------------*/
      std::cout << "\n\nI'm in here Dipole Module\n\n" << std::flush;

      //for(unsigned int i=0; i<cells::index_atoms_array1D.size(); i++){
      //   //std::cout << "cells::index_atoms_array1D[" << i << "] = " << cells::index_atoms_array1D[i] << std::endl << std::flush;
      //   //fprintf(stderr,"cells::index_atoms_array1D[%d] = %d on my_rank %d\n",i,cells::index_atoms_array1D[i],vmpi::my_rank);
      //}
      //std::cout << std::endl << std::flush;
      //MPI::COMM_WORLD.Barrier();


      #ifdef MPICF
         std::vector<double> atom_pos_x(num_local_atoms,0.0);
         std::vector<double> atom_pos_y(num_local_atoms,0.0);
         std::vector<double> atom_pos_z(num_local_atoms,0.0);
         //std::vector<double> atom_mom(num_local_atoms,0.0);
         for(int atom=0; atom<num_local_atoms; atom++){
            atom_pos_x[atom]=atom_coords_x[atom];
            atom_pos_y[atom]=atom_coords_y[atom];
            atom_pos_z[atom]=atom_coords_z[atom];
            //int type = atom_type_array[atom];
            //const double mus  = mp::material[type].mu_s_SI/9.27400915e-24;
            ////////fprintf(stderr,"  atom = %d type = %d mus = %e atoms::type_array[atom] = %d mp::material[type] = %f on my_rank = %d\n",atom,type,mus,atoms::type_array[atom],mp::material[type].mu_s_SI,vmpi::my_rank);
            //atom_mom[atom] = mus;
         }
         //for(int atom=0; atom<num_local_atoms; atom++){
         //   fprintf(stderr," atom %d x = %f y = %f z = %f in cell %d on my_rank = %d\n",atom,atom_pos_x[atom],atom_pos_y[atom],atom_pos_z[atom],dipole::internal::atom_cell_id_array[atom],vmpi::my_rank);
         //}

         //for(int i=0; i<dipole::internal::cells_num_cells; i++){
           // fprintf(stderr,"\n\n  >>>> Print volume array before paralle section <<<<<<\n   !!!!!!cells_volume_array[%d] = %f cells_num_atoms_in_cell[%d] = %d on rank = %d\n",i,dipole::internal::cells_volume_array[i],i,dipole::internal::cells_num_atoms_in_cell[i],vmpi::my_rank);
         //}

         //MPI::COMM_WORLD.Barrier();

         for(int lc=0; lc<dipole::internal::cells_num_cells; lc++){
            //fprintf(stderr,"\t *#*#*#* lc = %d  num_atoms_in_cell[%d] = %d *#*#*#\n",lc,lc,dipole::internal::cells_num_atoms_in_cell[lc]);
            // resize arrays
            cells_atom_in_cell_coords_array_x[lc].resize(dipole::internal::cells_num_atoms_in_cell[lc]);
            cells_atom_in_cell_coords_array_y[lc].resize(dipole::internal::cells_num_atoms_in_cell[lc]);
            cells_atom_in_cell_coords_array_z[lc].resize(dipole::internal::cells_num_atoms_in_cell[lc]);
            cells_index_atoms_array[lc].resize(dipole::internal::cells_num_atoms_in_cell[lc]);
            //fprintf(stderr,"cells_index_atoms_array[%d].size() = %lu cells_index_atoms_array[%d].size() = %lu on my_rank = %d\n",lc,cells_index_atoms_array[lc].size(),lc,cells_atom_in_cell_coords_array_x[lc].size(),vmpi::my_rank);
         }

         //MPI::COMM_WORLD.Barrier();
         // Call paralelisation function
         dipole::send_recv_cells_data(dipole::internal::proc_cell_index_array1D,
                                    cells_atom_in_cell_coords_array_x,
                                    cells_atom_in_cell_coords_array_y,
                                    cells_atom_in_cell_coords_array_z,
                                    cells_index_atoms_array,
                                    dipole::internal::cells_pos_and_mom_array,
                                    dipole::internal::cells_num_atoms_in_cell,
                                    cells::cell_id_array,
                                    dipole::internal::cells_local_cell_array,
                                    dipole::internal::cells_num_local_cells,
                                    dipole::internal::cells_num_cells
                                    );

         //MPI::COMM_WORLD.Barrier();
         //fprintf(stderr,"\n\n >>>>> End of send/recv session for cells ---> cells_num_cells = %d, cells_num_local_cells = %d on rank = %d<<<<<<\n\n",dipole::internal::cells_num_cells,dipole::internal::cells_num_local_cells,vmpi::my_rank);

         int size = cells::cell_id_array.size();
         //for(int i=0; i<dipole::internal::cells_num_local_cells; i++){
         //   //int lc = cells::cell_id_array[i];
         //   int lc = dipole::internal::cells_local_cell_array[i];
         //   ////////fprintf(stderr,"size = %d, cells::cell_id_array[%d] = %d, x = %f y = %f z = %f mu = %e on cpu = %d proc_cell_index_array[%d] = %d\n",size,i,cells::cell_id_array[i],cells::pos_and_mom_array[4*lc+0],cells::pos_and_mom_array[4*lc+1],cells::pos_and_mom_array[4*lc+2],cells::pos_and_mom_array[4*lc+3],vmpi::my_rank,lc,dipole::internal::proc_cell_index_array1D[lc]);
         //   fprintf(stderr,"size = %d, cells_local_cell_array[%d] = %d, x = %f y = %f z = %f mu = %e on cpu = %d proc_cell_index_array[%d] = %d\n",size,i,dipole::internal::cells_local_cell_array[i],cells::pos_and_mom_array[4*lc+0],cells::pos_and_mom_array[4*lc+1],cells::pos_and_mom_array[4*lc+2],cells::pos_and_mom_array[4*lc+3],vmpi::my_rank,lc,dipole::internal::proc_cell_index_array1D[lc]);
         //}
         //for(int lc=0; lc<ceil(dipole::internal::cells_pos_and_mom_array.size()/4.); lc++){
         //   if(dipole::internal::cells_num_atoms_in_cell[lc]>0 || lc>=dipole::internal::cells_num_cells){
         //      fprintf(stderr,"cell=%d, num_atoms_in_cell=%d, x = %f y = %f z = %f on my_rank=%d proc_cell_index_array[%d]=%d\n",lc,dipole::internal::cells_num_atoms_in_cell[lc],dipole::internal::cells_pos_and_mom_array[4*lc+0],dipole::internal::cells_pos_and_mom_array[4*lc+1],dipole::internal::cells_pos_and_mom_array[4*lc+2],vmpi::my_rank,lc,dipole::internal::proc_cell_index_array1D[lc]);
         //   }
         //}

         //MPI::COMM_WORLD.Barrier();
         //fprintf(stderr,"\n\n >>>>> Beginning of send/recv session for atoms ---> cells_num_cells = %d, cells_num_local_cells = %d cells_pos_and_mom_array.size() = %lu on rank = %d<<<<<<\n\n",dipole::internal::cells_num_cells,dipole::internal::cells_num_local_cells,dipole::internal::cells_pos_and_mom_array.size(),vmpi::my_rank);

         dipole::send_recv_atoms_data(dipole::internal::proc_cell_index_array1D,
                                    cells::cell_id_array,
                                    dipole::internal::cells_local_cell_array,
                                    atom_pos_x,
                                    atom_pos_y,
                                    atom_pos_z,
                                    dipole::internal::atom_type_array, // atomic moments (from dipole;:internal::atom_type_array)
                                    cells_atom_in_cell_coords_array_x,
                                    cells_atom_in_cell_coords_array_y,
                                    cells_atom_in_cell_coords_array_z,
                                    cells_index_atoms_array,
                                    dipole::internal::cells_pos_and_mom_array,
                                    dipole::internal::cells_num_atoms_in_cell,
                                    dipole::internal::cells_num_local_cells,
                                    dipole::internal::cells_num_cells,
                                    cells_macro_cell_size
                                    );

         //MPI::COMM_WORLD.Barrier();
         //fprintf(stderr,"\n\n >>>>> End of send/recv session for atoms ---> cells_num_cells = %d, cells_num_local_cells = %d cells_pos_and_mom_array.size() = %lu on rank = %d<<<<<<\n\n",dipole::internal::cells_num_cells,dipole::internal::cells_num_local_cells,dipole::internal::cells_pos_and_mom_array.size(),vmpi::my_rank);
         //MPI::COMM_WORLD.Barrier();
         //for(int lc=0; lc<ceil(dipole::internal::cells_pos_and_mom_array.size()/4.); lc++){
         //   if( (dipole::internal::cells_num_atoms_in_cell[lc]>0 && lc<dipole::internal::cells_num_cells) || (lc>=dipole::internal::cells_num_cells)){
         //      fprintf(stderr,"cell=%d, num_atoms_in_cell=%d, x = %f y = %f z = %f on my_rank=%d proc_cell_index_array[%d]=%d\n",lc,dipole::internal::cells_num_atoms_in_cell[lc],dipole::internal::cells_pos_and_mom_array[4*lc+0],dipole::internal::cells_pos_and_mom_array[4*lc+1],dipole::internal::cells_pos_and_mom_array[4*lc+2],vmpi::my_rank,lc,dipole::internal::proc_cell_index_array1D[lc]);
         //   }
         //}
         //MPI::COMM_WORLD.Barrier();
         //fprintf(stderr,"\n\n >>>>> BEFORE sort_data() func <<<<<<<\n");

      	sort_data(dipole::internal::proc_cell_index_array1D,
                  cells::cell_id_array,
                  cells_atom_in_cell_coords_array_x,
                  cells_atom_in_cell_coords_array_y,
                  cells_atom_in_cell_coords_array_z,
                  cells_index_atoms_array,
                  dipole::internal::cells_pos_and_mom_array,
                  dipole::internal::cells_num_atoms_in_cell,
                  dipole::internal::cells_num_local_cells,
                  dipole::internal::cells_num_cells
                  );
         //MPI::COMM_WORLD.Barrier();
		   //fprintf(stderr,"\n\n =======> after call of sort_data() func ---> cells_num_cells = %d, cells_num_local_cells = %d cells_pos_and_mom_array.size() = %lu on rank = %d <===========\n\n",dipole::internal::cells_num_cells,dipole::internal::cells_num_local_cells,dipole::internal::cells_pos_and_mom_array.size(),vmpi::my_rank);
         //MPI::COMM_WORLD.Barrier();
         //for(int lc=0; lc<ceil(dipole::internal::cells_pos_and_mom_array.size()/4.); lc++){
         //   if( (dipole::internal::cells_num_atoms_in_cell[lc]>0 && lc<dipole::internal::cells_num_cells) || (lc>=dipole::internal::cells_num_cells)){
         //      fprintf(stderr,"cell=%d, num_atoms_in_cell=%d, x = %f y = %f z = %f on my_rank=%d proc_cell_index_array[%d]=%d\n",lc,dipole::internal::cells_num_atoms_in_cell[lc],dipole::internal::cells_pos_and_mom_array[4*lc+0],dipole::internal::cells_pos_and_mom_array[4*lc+1],dipole::internal::cells_pos_and_mom_array[4*lc+2],vmpi::my_rank,lc,dipole::internal::proc_cell_index_array1D[lc]);
         //   }
         //}
         //MPI::COMM_WORLD.Barrier();

         //fprintf(stderr,"\n\n new cells_num_cells = %d cell_id_array.size()=%lu  tmp_cells_size=%d on my_rank = %d\n\n",cells_num_cells,tmp_cells_size,vmpi::my_rank);
         // Clear atom_pos_x,y,z
         atom_pos_x.clear();
         atom_pos_y.clear();
         atom_pos_z.clear();
      #endif
      //MPI::COMM_WORLD.Barrier();
      //fprintf(stderr,"\n\n");

      //fprintf(stderr,"\n !!!! PROBLEMS!!! \n");

		//MPI::COMM_WORLD.Barrier();

      //for(int i=0; i<dipole::internal::cells_num_cells; i++){
      //   //fprintf(stderr,"\n\n  >>>> Print volume array after paralle section <<<<<<\n   !!!!!!cells_volume_array[%d] = %f cells_num_atoms_in_cell[%d] = %d on rank = %d\n",i,dipole::internal::cells_volume_array[i],i,dipole::internal::cells_num_atoms_in_cell[i],vmpi::my_rank);
      //}


      //fprintf(stderr,"\n !!!! PROBLEMS!!!  just after updating dipoolenum_local_cells\n");

      //fprintf(stderr,"\n !!!! PROBLEMS!!! just before resize arrays dipole matrix\n");

      //for(int i=0;i<dipole::internal::cells_num_cells;i++){
      //   if(dipole::internal::cells_num_atoms_in_cell[i]>0){
      //      for(int j=0;j<dipole::internal::cells_num_atoms_in_cell[i];j++){
      //         //std::cout << "\t";
      //      	//std::cout << cells::index_atoms_array[i][j] << "\t";
      //      	//std::cout << i << "\t";
      //      	//std::cout << j << "\t";
      //      	//std::cout << cells_atom_in_cell_coords_array_x[i][j] << "\t";
      //      	//std::cout << cells_atom_in_cell_coords_array_y[i][j] << "\t";
      //      	//std::cout << cells_atom_in_cell_coords_array_z[i][j] << std::endl;
      //         fprintf(stderr," index_atoms_array[%d][%d] = %d atom_in_cell_coords_array_x[%d][%d] = %f atom_in_cell_coords_array_y[%d][%d] = %f atom_in_cell_coords_array_z[%d][%d] = %f on rank = %d\n",i,j,cells::index_atoms_array[i][j],i,j,cells_atom_in_cell_coords_array_x[i][j],i,j,cells_atom_in_cell_coords_array_y[i][j],i,j,cells_atom_in_cell_coords_array_z[i][j],vmpi::my_rank);
      //      }
      //   }
      //}
      //MPI::COMM_WORLD.Barrier();
      //for(int lc=0; lc<ceil(dipole::internal::cells_pos_and_mom_array.size()/4.); lc++){
      //   if( (dipole::internal::cells_num_atoms_in_cell[lc]>0 && lc<dipole::internal::cells_num_cells) || (lc>=dipole::internal::cells_num_cells)){
      //      fprintf(stderr,"cell=%d, num_atoms_in_cell=%d, x = %f y = %f z = %f on my_rank=%d proc_cell_index_array[%d]=%d\n",lc,dipole::internal::cells_num_atoms_in_cell[lc],dipole::internal::cells_pos_and_mom_array[4*lc+0],dipole::internal::cells_pos_and_mom_array[4*lc+1],dipole::internal::cells_pos_and_mom_array[4*lc+2],vmpi::my_rank,lc,dipole::internal::proc_cell_index_array1D[lc]);
      //   }
      //}
      //MPI::COMM_WORLD.Barrier();
      //std::cout << std::endl << std::flush;

       // allocate arrays to store data [nloccell x ncells]
      for(int lc=0;lc<dipole::internal::cells_num_local_cells; lc++){

         dipole::internal::rij_inter_xx.push_back(std::vector<double>());
         dipole::internal::rij_inter_xx[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_inter_xy.push_back(std::vector<double>());
         dipole::internal::rij_inter_xy[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_inter_xz.push_back(std::vector<double>());
         dipole::internal::rij_inter_xz[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_inter_yy.push_back(std::vector<double>());
         dipole::internal::rij_inter_yy[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_inter_yz.push_back(std::vector<double>());
         dipole::internal::rij_inter_yz[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_inter_zz.push_back(std::vector<double>());
         dipole::internal::rij_inter_zz[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_intra_xx.push_back(std::vector<double>());
         dipole::internal::rij_intra_xx[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_intra_xy.push_back(std::vector<double>());
         dipole::internal::rij_intra_xy[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_intra_xz.push_back(std::vector<double>());
         dipole::internal::rij_intra_xz[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_intra_yy.push_back(std::vector<double>());
         dipole::internal::rij_intra_yy[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_intra_yz.push_back(std::vector<double>());
         dipole::internal::rij_intra_yz[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::internal::rij_intra_zz.push_back(std::vector<double>());
         dipole::internal::rij_intra_zz[lc].resize(dipole::internal::cells_num_cells,0.0);

         dipole::cells_field_array_x.resize(dipole::internal::cells_num_cells,0.0);
         dipole::cells_field_array_y.resize(dipole::internal::cells_num_cells,0.0);
         dipole::cells_field_array_z.resize(dipole::internal::cells_num_cells,0.0);

      }

      // calculate matrix prefactors
      zlog << zTs() << "Precalculating rij matrix for demag calculation... " << std::endl;


      //std::cout<< "Number of local cells= "<<dipole::internal::cells_num_local_cells << std::endl;
      //std::cout<< "Number of  cells= "<<dipole::internal::cells_num_cells << std::endl;

      //double cutoff=12.0; //after 12 macrocell of distance, the bare macrocell model gives the same result

      // loop over local cells
      for(int lc=0;lc<dipole::internal::cells_num_local_cells;lc++){

         // reference global cell ID
         //int i = dipole::internal::cells_local_cell_array[lc];
         int i = cells::cell_id_array[lc];
         //fprintf(stderr,"lc = %d, i = %d, x = %f y = %f z = %f M = %e proc_cell_index_array1D[i] = %d on rank = %d\n",lc,i,dipole::internal::cells_pos_and_mom_array[4*i+0],dipole::internal::cells_pos_and_mom_array[4*i+1],dipole::internal::cells_pos_and_mom_array[4*i+2],dipole::internal::cells_pos_and_mom_array[4*i+3], dipole::internal::proc_cell_index_array1D[i],vmpi::my_rank);

         if(dipole::internal::cells_num_atoms_in_cell[i]>0){

         	// Loop over all other cells to calculate contribution to local cell
            for(int j=0;j<dipole::internal::cells_num_cells;j++){

             	if(i!=j && dipole::internal::cells_num_atoms_in_cell[j]>0){
                  //fprintf(stderr," >> inter >> i = %d, x = %f y = %f z = %f\tj = %d, x = %f y = %f z = %f proc_cell_index_array1D[j] = %d on rank = %d\n",i,dipole::internal::cells_pos_and_mom_array[4*i+0],dipole::internal::cells_pos_and_mom_array[4*i+1],dipole::internal::cells_pos_and_mom_array[4*i+2],j,dipole::internal::cells_pos_and_mom_array[4*j+0],dipole::internal::cells_pos_and_mom_array[4*j+1],dipole::internal::cells_pos_and_mom_array[4*j+2], dipole::internal::proc_cell_index_array1D[j],vmpi::my_rank);

                 	double tmp_rij_inter_xx = 0.0;
                 	double tmp_rij_inter_xy = 0.0;
                 	double tmp_rij_inter_xz = 0.0;

                 	double tmp_rij_inter_yy = 0.0;
                 	double tmp_rij_inter_yz = 0.0;
                 	double tmp_rij_inter_zz = 0.0;

                 	/*double rx = cells_cell_coords_array_x[j] - cells_cell_coords_array_x[i]; // Angstroms
                 	double ry = cells_cell_coords_array_y[j] - cells_cell_coords_array_y[i];
                 	double rz = cells_cell_coords_array_z[j] - cells_cell_coords_array_z[i];*/

                  double rx = dipole::internal::cells_pos_and_mom_array[4*j+0] - dipole::internal::cells_pos_and_mom_array[4*i+0];
                  double ry = dipole::internal::cells_pos_and_mom_array[4*j+1] - dipole::internal::cells_pos_and_mom_array[4*i+1];
                  double rz = dipole::internal::cells_pos_and_mom_array[4*j+2] - dipole::internal::cells_pos_and_mom_array[4*i+2];

                 	double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance
                 	double rij_1 = 1.0/rij;
                  //fprintf(stderr," >> i %d j %d rij/mc_size = %f on rank %d\n",i,j,(rij_1)/cells_macro_cell_size,vmpi::my_rank);

                  // If distance between macro-cells > cutoff nm => continuum approach
                  if( (rij_1)/cells_macro_cell_size > dipole::cutoff){
                     //fprintf(stderr,"\n >> inter >> i = %d, x = %f y = %f z = %f\tj = %d, x = %f y = %f z = %f proc_cell_index_array1D[j] = %d rij/mc_size = %f on rank = %d\n",i,dipole::internal::cells_pos_and_mom_array[4*i+0],dipole::internal::cells_pos_and_mom_array[4*i+1],dipole::internal::cells_pos_and_mom_array[4*i+2],j,dipole::internal::cells_pos_and_mom_array[4*j+0],dipole::internal::cells_pos_and_mom_array[4*j+1],dipole::internal::cells_pos_and_mom_array[4*j+2], dipole::internal::proc_cell_index_array1D[j],(rij_1)/cells_macro_cell_size,vmpi::my_rank);

                  	double ex = rx*rij;
                  	double ey = ry*rij;
                  	double ez = rz*rij;

                  	double rij3 = (rij*rij*rij); // Angstroms

                  	dipole::internal::rij_inter_xx[lc][j] = ((3.0*ex*ex - 1.0)*rij3);
                  	dipole::internal::rij_inter_xy[lc][j] = ( 3.0*ex*ey      )*rij3 ;
                  	dipole::internal::rij_inter_xz[lc][j] = ( 3.0*ex*ez      )*rij3 ;

                  	dipole::internal::rij_inter_yy[lc][j] = ((3.0*ey*ey - 1.0)*rij3);
                  	dipole::internal::rij_inter_yz[lc][j] = ( 3.0*ey*ez      )*rij3 ;
                  	dipole::internal::rij_inter_zz[lc][j] = ((3.0*ez*ez - 1.0)*rij3);
                     //fprintf(stderr,"lc = %d j = %d rij3=%f\n ex=%f\tey=%f\tez=%f\n\trij_inter_xx = %f\trij_inter_xy = %f\trij_inter_xz = %f\n\trij_inter_yx = %f\trij_inter_yy = %f\trij_inter_yz = %f\n\trij_inter_zx = %f\trij_inter_zy = %f\trij_inter_zz = %f\n on rank = %d\n",lc,j,rij3,ex,ey,ez,internal::rij_inter_xx[lc][j],internal::rij_inter_xy[lc][j],internal::rij_inter_xz[lc][j],internal::rij_inter_xy[lc][j],internal::rij_inter_yy[lc][j],internal::rij_inter_yz[lc][j],internal::rij_inter_xz[lc][j],internal::rij_inter_yz[lc][j],internal::rij_inter_zz[lc][j],vmpi::my_rank);
                  }
                  else if( (1.0/rij)/cells_macro_cell_size <= dipole::cutoff){
                     for(int pi=0; pi<dipole::internal::cells_num_atoms_in_cell[i]; pi++){
                        for(int qj=0; qj<dipole::internal::cells_num_atoms_in_cell[j]; qj++){

                           rx = cells_atom_in_cell_coords_array_x[j][qj] - cells_atom_in_cell_coords_array_x[i][pi];
                           ry = cells_atom_in_cell_coords_array_y[j][qj] - cells_atom_in_cell_coords_array_y[i][pi];
                           rz = cells_atom_in_cell_coords_array_z[j][qj] - cells_atom_in_cell_coords_array_z[i][pi];

                           rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz);  //Reciprocal of the distance
                           rij_1 = 1.0/rij;

                           if( rij_1==0.0 ) {
                              std::cout << ">>>>> (Inter) Warning: atoms are overlapping in cells i,j " << i << "\t" << j <<"<<<<<" << std::endl;
                              std::cout << "Cell and atomic coordinates used for calculation of Dipolar matrix" << std::endl;
                              /*std::cout << " xj= "  << cells_cell_coords_array_x[j]  << " yj= "  << cells_cell_coords_array_y[j]  << " zj= "  << cells_cell_coords_array_z[j] << std::endl;
                              std::cout << " xi= "  << cells_cell_coords_array_x[i]  << " yi= "  << cells_cell_coords_array_y[i]  << " zi= "  << cells_cell_coords_array_z[i] << std::endl;*/
                              std::cout << " xj= "  << dipole::internal::cells_pos_and_mom_array[4*j+0]  << " yj= "  << dipole::internal::cells_pos_and_mom_array[4*j+1]  << " zj= "  << dipole::internal::cells_pos_and_mom_array[4*j+2] << std::endl;
                              std::cout << " xi= "  << dipole::internal::cells_pos_and_mom_array[4*i+0]  << " yi= "  << dipole::internal::cells_pos_and_mom_array[4*i+1]  << " zi= "  << dipole::internal::cells_pos_and_mom_array[4*i+2] << std::endl;
//                            std::cout << " rx_cell= "  << rx_cell << " ry_cell= " << ry_cell << " rz_cell= " << rz_cell << std::endl;
                              std::cout << " xpi= " << cells_atom_in_cell_coords_array_x[i][pi] << " ypi= " << cells_atom_in_cell_coords_array_y[i][pi] << " zpi= " << cells_atom_in_cell_coords_array_z[i][pi] << std::endl;
                              std::cout << " xqj= " << cells_atom_in_cell_coords_array_x[j][qj] << " yqj= " << cells_atom_in_cell_coords_array_y[j][qj] << " zqj= " << cells_atom_in_cell_coords_array_z[j][qj] << std::endl;
                              std::cout << " rx= "  << rx << " ry= " << ry << " rz= " << rz << std::endl;
                              std::cout << " rij= " << 1.0/rij << " = " << (1.0/rij)/((cells_macro_cell_size*sqrt(3)+0.01)) << " macro-cells" << std::endl;
                              return;
                           }

                           const double ex = rx*rij;
                           const double ey = ry*rij;
                           const double ez = rz*rij;

                           double rij3 = (rij*rij*rij); // Angstroms

                           tmp_rij_inter_xx += ((3.0*ex*ex - 1.0)*rij3);
                           tmp_rij_inter_xy += ((3.0*ex*ey      )*rij3);
                           tmp_rij_inter_xz += ((3.0*ex*ez      )*rij3);

                           tmp_rij_inter_yy += ((3.0*ey*ey - 1.0)*rij3);
                           tmp_rij_inter_yz += ((3.0*ey*ez      )*rij3);
                           tmp_rij_inter_zz += ((3.0*ez*ez - 1.0)*rij3);
                           //fprintf(stderr,"i=%d %f %f %f,j=%d %f %f %f, rank=%d\n\t%f\t%f\t%f\n\t%f\t%f\t%f\n\t%f\t%f\t%f\n",i,cells_atom_in_cell_coords_array_x[i][pi],cells_atom_in_cell_coords_array_y[i][pi],cells_atom_in_cell_coords_array_z[i][pi],j,cells_atom_in_cell_coords_array_x[j][qj],cells_atom_in_cell_coords_array_y[j][qj],cells_atom_in_cell_coords_array_z[j][qj],vmpi::my_rank,tmp_rij_inter_xx,tmp_rij_inter_xy,tmp_rij_inter_xz,tmp_rij_inter_xy,tmp_rij_inter_yy,tmp_rij_inter_yz,tmp_rij_inter_yz,tmp_rij_inter_zz,tmp_rij_inter_xz);
                        }
                     }

                     dipole::internal::rij_inter_xx[lc][j] =  (tmp_rij_inter_xx);
                     dipole::internal::rij_inter_xy[lc][j] =  (tmp_rij_inter_xy);
                     dipole::internal::rij_inter_xz[lc][j] =  (tmp_rij_inter_xz);

                     dipole::internal::rij_inter_yy[lc][j] =  (tmp_rij_inter_yy);
                     dipole::internal::rij_inter_yz[lc][j] =  (tmp_rij_inter_yz);
                     dipole::internal::rij_inter_zz[lc][j] =  (tmp_rij_inter_zz);

                     dipole::internal::rij_inter_xx[lc][j] = dipole::internal::rij_inter_xx[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                     dipole::internal::rij_inter_xy[lc][j] = dipole::internal::rij_inter_xy[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                     dipole::internal::rij_inter_xz[lc][j] = dipole::internal::rij_inter_xz[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));

                     dipole::internal::rij_inter_yy[lc][j] = dipole::internal::rij_inter_yy[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                     dipole::internal::rij_inter_yz[lc][j] = dipole::internal::rij_inter_yz[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                     dipole::internal::rij_inter_zz[lc][j] = dipole::internal::rij_inter_zz[lc][j]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                  }  // End of Inter part calculated atomicstically
                  //fprintf(stderr,"lc = %d j = %d\n\trij_inter_xx = %f\trij_inter_xy = %f\trij_inter_xz = %f\n\trij_inter_yx = %f\trij_inter_yy = %f\trij_inter_yz = %f\n\trij_inter_zx = %f\trij_inter_zy = %f\trij_inter_zz = %f\n on rank = %d\n",lc,j,internal::rij_inter_xx[lc][j],internal::rij_inter_xy[lc][j],internal::rij_inter_xz[lc][j],internal::rij_inter_xy[lc][j],internal::rij_inter_yy[lc][j],internal::rij_inter_yz[lc][j],internal::rij_inter_xz[lc][j],internal::rij_inter_yz[lc][j],internal::rij_inter_zz[lc][j],vmpi::my_rank);
               } // End of Inter part

               // Calculation of INTRA TERM of Dipolar interaction
               //#ifdef MPICF
               //// if the cell is not splitted on more cpus
               ////else if( ((i==j && proc_cell_index_array1D[i]==proc_cell_index_array1D[j]) && dipole::internal::cells_num_atoms_in_cell[j]>0) ){
             	//else if( ( (i==j && dipole::internal::proc_cell_index_array1D[i]==dipole::internal::proc_cell_index_array1D[j]) ||
               //      (dipole::internal::proc_cell_index_array1D[i]!=dipole::internal::proc_cell_index_array1D[j] &&
               //            ((dipole::internal::cells_pos_and_mom_array[4*i+0]==dipole::internal::cells_pos_and_mom_array[4*j+0]) &&
               //             (dipole::internal::cells_pos_and_mom_array[4*i+1]==dipole::internal::cells_pos_and_mom_array[4*j+1]) &&
               //             (dipole::internal::cells_pos_and_mom_array[4*i+2]==dipole::internal::cells_pos_and_mom_array[4*j+2]) )
               //      ) ) && dipole::internal::cells_num_atoms_in_cell[j]>0){
               //#else
               else if( i==j && dipole::internal::cells_num_atoms_in_cell[j]>0){
               //#endif
                  //fprintf(stderr," >> intra >> j = %d, x = %f y = %f z = %f M = %e proc_cell_index_array1D[j] = %d on rank = %d\n",j,dipole::internal::cells_pos_and_mom_array[4*j+0],dipole::internal::cells_pos_and_mom_array[4*j+1],dipole::internal::cells_pos_and_mom_array[4*j+2],dipole::internal::cells_pos_and_mom_array[4*j+3],dipole::internal::proc_cell_index_array1D[j],vmpi::my_rank);

                  double tmp_rij_intra_xx = 0.0;
                  double tmp_rij_intra_xy = 0.0;
                  double tmp_rij_intra_xz = 0.0;

                  double tmp_rij_intra_yy = 0.0;
                  double tmp_rij_intra_yz = 0.0;
                  double tmp_rij_intra_zz = 0.0;

               	for(int pi=0; pi<dipole::internal::cells_num_atoms_in_cell[i]; pi++){
                 		for(int qj=0; qj<dipole::internal::cells_num_atoms_in_cell[i]; qj++){
                        //#ifdef MPICF
                        //if( (dipole::internal::proc_cell_index_array1D[i]==dipole::internal::proc_cell_index_array1D[j] && pi!=qj) ||
                        //    (dipole::internal::proc_cell_index_array1D[i]!=dipole::internal::proc_cell_index_array1D[j]) ){
                        //#else
                  		if( pi!=qj ){
                        //#endif

                   			//double rx = cells_atom_in_cell_coords_array_x[i][qj] - cells_atom_in_cell_coords_array_x[i][pi];
                   			//double ry = cells_atom_in_cell_coords_array_y[i][qj] - cells_atom_in_cell_coords_array_y[i][pi];
                   			//double rz = cells_atom_in_cell_coords_array_z[i][qj] - cells_atom_in_cell_coords_array_z[i][pi];

                   			double rx = cells_atom_in_cell_coords_array_x[j][qj] - cells_atom_in_cell_coords_array_x[i][pi];
                   			double ry = cells_atom_in_cell_coords_array_y[j][qj] - cells_atom_in_cell_coords_array_y[i][pi];
                   			double rz = cells_atom_in_cell_coords_array_z[j][qj] - cells_atom_in_cell_coords_array_z[i][pi];

                   			const double rij = 1.0/sqrt(rx*rx+ry*ry+rz*rz); //Reciprocal of the distance

                   			if( 1.0/rij==0.0 ){
                     			std::cout << ">>>>> (Intra)  Warning: atoms are overlapping in cells i=\t" << i << "\t and j=\t" << j <<"\t<<<<<" << std::endl;
                     			std::cout << "Cell and atomic coordinates used for calculation of Dipolar matrix" << std::endl;
                     			std::cout << " xqj= " << cells_atom_in_cell_coords_array_x[j][qj] << " yqj= " << cells_atom_in_cell_coords_array_y[j][qj] << " zqj= " << cells_atom_in_cell_coords_array_z[j][qj] << std::endl;
                     			std::cout << " xpi= " << cells_atom_in_cell_coords_array_x[i][pi] << " ypi= " << cells_atom_in_cell_coords_array_y[i][pi] << " zpi= " << cells_atom_in_cell_coords_array_z[i][pi] << std::endl;
                              /*std::cout << " xj= "  << cells_cell_coords_array_x[j]  << " yj= "  << cells_cell_coords_array_y[j]  << " zj= "  << cells_cell_coords_array_z[j] << std::endl;
                              std::cout << " xi= "  << cells_cell_coords_array_x[i]  << " yi= "  << cells_cell_coords_array_y[i]  << " zi= "  << cells_cell_coords_array_z[i] << std::endl;*/
                              std::cout << " xj= "  << dipole::internal::cells_pos_and_mom_array[4*j+0]  << " yj= "  << dipole::internal::cells_pos_and_mom_array[4*j+1]  << " zj= "  << dipole::internal::cells_pos_and_mom_array[4*j+2] << std::endl;
                              std::cout << " xi= "  << dipole::internal::cells_pos_and_mom_array[4*i+0]  << " yi= "  << dipole::internal::cells_pos_and_mom_array[4*i+1]  << " zi= "  << dipole::internal::cells_pos_and_mom_array[4*i+2] << std::endl;
                     			std::cout << " rx= "  << rx << " ry= " << ry << " rz= " << rz << std::endl;
                     			std::cout << " rij= " << 1.0/rij << " = " << (1.0/rij)/((cells_macro_cell_size*sqrt(3)+0.01)) << " macro-cells" << std::endl;
                     			return;
                   			}

			                  const double ex = rx*rij;
                   			const double ey = ry*rij;
                   			const double ez = rz*rij;

                   			const double rij3 = (rij*rij*rij); // Angstroms

                   			tmp_rij_intra_xx += ((3.0*ex*ex - 1.0)*rij3);
                   			tmp_rij_intra_xy += ((3.0*ex*ey      )*rij3);
                   			tmp_rij_intra_xz += ((3.0*ex*ez      )*rij3);

                   			tmp_rij_intra_yy += ((3.0*ey*ey - 1.0)*rij3);
                   			tmp_rij_intra_yz += ((3.0*ey*ez      )*rij3);
                   			tmp_rij_intra_zz += ((3.0*ez*ez - 1.0)*rij3);

                 			}
                 		}
                	}

                	dipole::internal::rij_intra_xx[lc][i] =  (tmp_rij_intra_xx);
                	dipole::internal::rij_intra_xy[lc][i] =  (tmp_rij_intra_xy);
                	dipole::internal::rij_intra_xz[lc][i] =  (tmp_rij_intra_xz);

                	dipole::internal::rij_intra_yy[lc][i] =  (tmp_rij_intra_yy);
                	dipole::internal::rij_intra_yz[lc][i] =  (tmp_rij_intra_yz);
                	dipole::internal::rij_intra_zz[lc][i] =  (tmp_rij_intra_zz);

                	dipole::internal::rij_intra_xx[lc][i] = dipole::internal::rij_intra_xx[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                	dipole::internal::rij_intra_xy[lc][i] = dipole::internal::rij_intra_xy[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                	dipole::internal::rij_intra_xz[lc][i] = dipole::internal::rij_intra_xz[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));

                	dipole::internal::rij_intra_yy[lc][i] = dipole::internal::rij_intra_yy[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                	dipole::internal::rij_intra_yz[lc][i] = dipole::internal::rij_intra_yz[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
                	dipole::internal::rij_intra_zz[lc][i] = dipole::internal::rij_intra_zz[lc][i]/(double(dipole::internal::cells_num_atoms_in_cell[i]) * double(dipole::internal::cells_num_atoms_in_cell[j]));
               }
               //fprintf(stderr,"lc = %d j = %d i = %d\n\trij_intra_xx = %f\trij_intra_xy = %f\trij_intra_xz = %f\n\trij_intra_yx = %f\trij_intra_yy = %f\trij_intra_yz = %f\n\trij_intra_zx = %f\trij_intra_zy = %f\trij_intra_zz = %f\n on rank = %d\n",lc,j,i,internal::rij_intra_xx[lc][i],internal::rij_intra_xy[lc][i],internal::rij_intra_xz[lc][i],internal::rij_intra_xy[lc][i],internal::rij_intra_yy[lc][i],internal::rij_intra_yz[lc][i],internal::rij_intra_xz[lc][i],internal::rij_intra_yz[lc][i],internal::rij_intra_zz[lc][i],vmpi::my_rank);
            }
			}
		}

      //MPI::COMM_WORLD.Barrier();
      //fprintf(stderr,"\n\n >>> printing inter part of dipolar matrix <<<< \n\n");
      //MPI::COMM_WORLD.Barrier();
      //// print to check dipolar martrix inter term
      //for(unsigned int i=0; i<dipole::internal::rij_inter_xx.size(); i++){
      //   if(dipole::internal::cells_num_atoms_in_cell[cells::cell_id_array[i]]>0){
      //      for(unsigned int j=0; j<dipole::internal::rij_inter_xx[i].size(); j++){
      //         if(dipole::internal::cells_num_atoms_in_cell[j]>0){
      //            fprintf(stderr,"i = %d j = %d\n\trij_inter_xx = %f\trij_inter_xy = %f\trij_inter_xz = %f\n\trij_inter_yx = %f\trij_inter_yy = %f\trij_inter_yz = %f\n\trij_inter_zx = %f\trij_inter_zy = %f\trij_inter_zz = %f\n on rank = %d\n",i,j,internal::rij_inter_xx[i][j],internal::rij_inter_xy[i][j],internal::rij_inter_xz[i][j],internal::rij_inter_xy[i][j],internal::rij_inter_yy[i][j],internal::rij_inter_yz[i][j],internal::rij_inter_xz[i][j],internal::rij_inter_yz[i][j],internal::rij_inter_zz[i][j],vmpi::my_rank);
      //         }
      //      }
      //   }
      //}
      //MPI::COMM_WORLD.Barrier();
      //fprintf(stderr,"\n\n >>> printing intra part of dipolar matrix <<<< \n\n");
      //MPI::COMM_WORLD.Barrier();
      //// print to check dipolar matrix intra term
      //for(unsigned int i=0; i<dipole::internal::rij_intra_xx.size(); i++){
      //   for(unsigned int j=0; j<dipole::internal::rij_intra_xx[i].size(); j++){
      //      fprintf(stderr,"i = %d j = %d\n\trij_intra_xx = %f\trij_intra_xy = %f\trij_intra_xz = %f\n\trij_intra_yx = %f\trij_intra_yy = %f\trij_intra_yz = %f\n\trij_intra_zx = %f\trij_intra_zy = %f\trij_intra_zz = %f\n on rank = %d\n",i,j,internal::rij_intra_xx[i][j],internal::rij_intra_xy[i][j],internal::rij_intra_xz[i][j],internal::rij_intra_xy[i][j],internal::rij_intra_yy[i][j],internal::rij_intra_yz[i][j],internal::rij_intra_xz[i][j],internal::rij_intra_yz[i][j],internal::rij_intra_zz[i][j],vmpi::my_rank);
      //   }
      //}
      //MPI::COMM_WORLD.Barrier();

	   ////------- CPUs OUTPUT Dij on different fiels ------------//
      //// Set local output filename
      //std::stringstream file_sstr;
      //file_sstr << "Dij-";
      //file_sstr << vmpi::my_rank;
      //// Set CPUID on non-root process
      ////if(vmpi::my_rank!=0){
      //   //file_sstr << std::setfill('0') << std::setw(5) << vmpi::my_rank << "-";
      ////}
      ////file_sstr << std::setfill('0') << std::setw(8) << output_atoms_file_counter;
      ////file_sstr << ".cfg";
      //std::string cfg_file = file_sstr.str();
      //const char* cfg_filec = cfg_file.c_str();
      /
      //// Output informative message to log file
      //zlog << zTs() << "Outputting dipole matrix file " << cfg_file << " to disk" << std::endl;
      /
      //// Declare and open output file
      //std::ofstream cfg_file_ofstr;
      //cfg_file_ofstr.open (cfg_filec);
      /
      //// Every cpus print to check dipolar martrix inter term
      //for(int i=0; i<dipole::internal::cells_num_local_cells; i++){
      //   int lc = cells::cell_id_array[i];
      //   if(dipole::internal::cells_num_atoms_in_cell[cells::cell_id_array[i]]>0){
      //      for(unsigned int j=0; j<dipole::internal::rij_inter_xx[i].size(); j++){
      //         if(dipole::internal::cells_num_atoms_in_cell[j]>0){
      //            //fprintf(cfg_file_ofstr,"i = %d j = %d\n\trij_inter_xx = %f\trij_inter_xy = %f\trij_inter_xz = %f\n\trij_inter_yx = %f\trij_inter_yy = %f\trij_inter_yz = %f\n\trij_inter_zx = %f\trij_inter_zy = %f\trij_inter_zz = %f\n on rank = %d\n",i,j,internal::rij_inter_xx[i][j],internal::rij_inter_xy[i][j],internal::rij_inter_xz[i][j],internal::rij_inter_xy[i][j],internal::rij_inter_yy[i][j],internal::rij_inter_yz[i][j],internal::rij_inter_xz[i][j],internal::rij_inter_yz[i][j],internal::rij_inter_zz[i][j],vmpi::my_rank);
      //            cfg_file_ofstr << "i = " << i << "\tlc = " << lc << "\t";
      //            cfg_file_ofstr << dipole::internal::cells_pos_and_mom_array[4*lc+0] << "\t";
      //            cfg_file_ofstr << dipole::internal::cells_pos_and_mom_array[4*lc+1] << "\t";
      //            cfg_file_ofstr << dipole::internal::cells_pos_and_mom_array[4*lc+2] << "\t";
      //            cfg_file_ofstr << " j = " << j << "\t";
      //            cfg_file_ofstr << dipole::internal::cells_pos_and_mom_array[4*j+0] << "\t";
      //            cfg_file_ofstr << dipole::internal::cells_pos_and_mom_array[4*j+1] << "\t";
      //            cfg_file_ofstr << dipole::internal::cells_pos_and_mom_array[4*j+2] << "\t";
      //            cfg_file_ofstr << "\n";
      //            cfg_file_ofstr << "rij_intra_xx = " << dipole::internal::rij_intra_xx[i][j] << "\trij_intra_xy = " << dipole::internal::rij_intra_xy[i][j] << "\trij_intra_xz = " << dipole::internal::rij_intra_xz[i][j] << "\n";
      //            cfg_file_ofstr << "rij_intra_yx = " << dipole::internal::rij_intra_xy[i][j] << "\trij_intra_yy = " << dipole::internal::rij_intra_yy[i][j] << "\trij_intra_yz = " << dipole::internal::rij_intra_yz[i][j] << "\n";
      //            cfg_file_ofstr << "rij_intra_zx = " << dipole::internal::rij_intra_xz[i][j] << "\trij_intra_zy = " << dipole::internal::rij_intra_yz[i][j] << "\trij_intra_zz = " << dipole::internal::rij_intra_zz[i][j] << "\n";
      //            cfg_file_ofstr << "\n";
      //            cfg_file_ofstr << "rij_inter_xx = " << dipole::internal::rij_inter_xx[i][j] << "\trij_inter_xy = " << dipole::internal::rij_inter_xy[i][j] << "\trij_inter_xz = " << dipole::internal::rij_inter_xz[i][j] << "\n";
      //            cfg_file_ofstr << "rij_inter_yx = " << dipole::internal::rij_inter_xy[i][j] << "\trij_inter_yy = " << dipole::internal::rij_inter_yy[i][j] << "\trij_inter_yz = " << dipole::internal::rij_inter_yz[i][j] << "\n";
      //            cfg_file_ofstr << "rij_inter_zx = " << dipole::internal::rij_inter_xz[i][j] << "\trij_inter_zy = " << dipole::internal::rij_inter_yz[i][j] << "\trij_inter_zz = " << dipole::internal::rij_inter_zz[i][j] << "\n";
      //            cfg_file_ofstr << "\n";
      //         }
      //      }
      //   }
      //}
      /
      //cfg_file_ofstr.close();
	   ////--------------------------------------------------/



      //dipole::internal::cells_pos_and_mom_array = cells_pos_and_mom_array;
      //dipole::internal::cells_num_atoms_in_cell = cells_num_atoms_in_cell;
      //dipole::internal::cells_num_cells = cells_num_cells;

      #ifdef MPICF
         double t2 = MPI_Wtime();
      #else
         time_t t2;
         t2 = time (NULL);
      #endif
      zlog << zTs() << "Precalculation of rij matrix for demag calculation complete. Time taken: " << t2-t1 << "s."<< std::endl;


      // Set initialised flag
      dipole::internal::initialised=true;

      // timing function
      #ifdef MPICF
			t1 = MPI_Wtime();
      #else
         //time_t t1;
         t1 = time (NULL);
      #endif

      //fprintf(stderr,"\n >>>>> PROBLEMS!!! \n");
      // now calculate fields
      dipole::calculate_field();

      // timing function
      #ifdef MPICF
         t2 = MPI_Wtime();
      #else
         //time_t t2;
         t2 = time (NULL);
      #endif

      zlog << zTs() << "Time required for demag update: " << t2-t1 << "s." << std::endl;

    	return;

    }

} // end of dipole namespace
