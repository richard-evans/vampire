//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//   implementation based on the paper Phys. Rev. B 103, 024429, (2021) M.Strungaru, M.O.A. Ellis et al
//------------------------------------------------------------------------------
//

#ifdef MPICF
// Standard Libraries
#include <cmath>
#include <iostream>
#include <vector>

// Vampire Header files
#include "errors.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include "sld.hpp"
#include "atoms.hpp"
#include "create.hpp"
#include "material.hpp"

// Internal header
#include "internal.hpp"







namespace sld{



void suzuki_trotter_parallel_init(std::vector<double> &x, // atomic coordinates
                      std::vector<double> &y,
                      std::vector<double> &z,
                      double min_dim[3], // minimum dimensions on local processor
                      double max_dim[3]){ // maximum dimensions on local processor


     // if (vmpi::my_rank ==0)std::cout<< "inside Suzuki Trotter Parallel init rank0 "<<max_dim[0]<<"\t"<<max_dim[1]<<"\t"<<max_dim[2]<<std::endl;
     // if (vmpi::my_rank ==1)std::cout<< "inside Suzuki Trotter Parallel init rank1"<<max_dim[0]<<"\t"<<max_dim[1]<<"\t"<<max_dim[2]<<std::endl;
     


   // Convenient shorthands
   int catoms = vmpi::num_core_atoms;
   int batoms = vmpi::num_bdry_atoms;
   
   if (vmpi::my_rank ==0) {
   std::cout<< "rank 0 inside  2 Suzuki Trotter Parallel init "<<vmpi::my_rank<<"\t"<<catoms<<"\t"<<batoms<<std::endl;
   }
 
   if (vmpi::my_rank ==1) std::cout<< "rank 1 inside  2 Suzuki Trotter Parallel init "<<vmpi::my_rank<<"\t"<<catoms<<"\t"<<batoms<<std::endl;


   double widthx = max_dim[0] - min_dim[0];
   double widthy = max_dim[1] - min_dim[1];
   double widthz = max_dim[2] - min_dim[2];

   if (widthx<=sld::internal::r_cut_fields ||widthy<=sld::internal::r_cut_fields || widthy<=sld::internal::r_cut_fields ){
     std::cerr << "Error: System size needs to be increased so octants won't interact" << std::endl;
     std::cerr<<"Width of octants is "<<widthx<<"\t"<<widthy<<"\t"<<widthz<< " While interaction cutoff is "<<sld::internal::r_cut_fields<<std::endl;
     err::vexit();
     }

   int octant_num = 0; //Count which octant loop is in

   // Determines which core atoms are in which octant and pushes the index of those
   // atoms into the appropriate octant arrays.
   for(int zoct=0; zoct<2; zoct++){
      for(int yoct=0; yoct<2; yoct++){
         for(int xoct=0; xoct<2; xoct++){
            // Loop through all core atoms
            for (int i=0; i<catoms; i++){
               // Check if current atom is in desired octant
               if (   x[i] >= min_dim[0] + widthx*xoct*0.5 && x[i] < min_dim[0] + widthx*0.5 + widthx*xoct*0.5
                   && y[i] >= min_dim[1] + widthy*yoct*0.5 && y[i] < min_dim[1] + widthy*0.5 + widthy*yoct*0.5
                   && z[i] >= min_dim[2] + widthz*zoct*0.5 && z[i] < min_dim[2] + widthz*0.5 + widthz*zoct*0.5)
               {
                  internal::c_octants[octant_num].push_back(i);
                 // std::cout<<"octant core "<<octant_num<<"\t"<<i<<"\t"<<x[i]<<"\t"<<y[i]<<"\t"<<z[i]<<std::endl;
               }
            }
            octant_num++;
         }
      }
   }
   

   octant_num = 0;
   //Sort boundary atoms into appropriate octant arrays.
   for(int zoct=0; zoct<2; zoct++){
      for(int yoct=0; yoct<2; yoct++){
         for(int xoct=0; xoct<2; xoct++){
            // Loop through all boundary atoms
            for (int i=catoms; i<catoms+batoms; i++){
               // Check if current atom is in desired octant
               if (   x[i] >= min_dim[0] + widthx*xoct*0.5 && x[i] < min_dim[0] + widthx*0.5 + widthx*xoct*0.5
                   && y[i] >= min_dim[1] + widthy*yoct*0.5 && y[i] < min_dim[1] + widthy*0.5 + widthy*yoct*0.5
                   && z[i] >= min_dim[2] + widthz*zoct*0.5 && z[i] < min_dim[2] + widthz*0.5 + widthz*zoct*0.5)
               {
                  internal::b_octants[octant_num].push_back(i);
                 // std::cout<<"octant bound "<<octant_num<<" \t"<<i<<"\t"<<x[i]<<"\t"<<y[i]<<"\t"<<z[i]<<std::endl;

               }
            }
            octant_num++;
         }
      }
   }

   //--------------------------------------------------------------------
   // check that all atoms have been allocated an octant
   //--------------------------------------------------------------------
   // core atoms
   int num_atoms_in_octants = 0;
   for(int i=0; i< 8; i++) {num_atoms_in_octants += internal::c_octants[i].size();
      std::cout<<i<<" num_atoms_in_octants c "<<num_atoms_in_octants<<"\t"<<catoms<<"\t"<<atoms::num_atoms<<std::endl;}

   if(num_atoms_in_octants != catoms){
      std::cerr << "Programmer error: missing atoms in core octants in suzuki-trotter initialisation" << std::endl;
      err::vexit();
   }
   // boundary atoms
   num_atoms_in_octants = 0;
   for(int i=0; i< 8; i++) {
   num_atoms_in_octants += internal::b_octants[i].size();
  // std::cout<<"i="<<i<<"\t"<<internal::b_octants[i].size()<<std::endl;}
      std::cout<<i <<" num_atoms_in_octants b "<<num_atoms_in_octants<<"\t"<<batoms<<"\t"<<atoms::num_atoms<<std::endl;}

   if(num_atoms_in_octants != batoms){
      std::cerr << "Programmer error: missing atoms in boundary octants in suzuki-trotter initialisation" << std::endl;
      err::vexit();
   }



   //create array for suzuki trotter integration index
   //initialise variables
   internal::all_atoms_octant_start_index.push_back(0);
   
   for (int j=0;j<2000;j++){
      
    for(int octant=0; octant< 8; octant++){
    
    int core_at=internal::c_octants[octant].size();
    int bdry_at=internal::b_octants[octant].size();
    
    for (int i=0; i<core_at;i++){
          	int atom = internal::c_octants[octant][i];
          	//std::cout<< " core oct and atom " <<octant<<"\t"<<atom<<"\t"<<x[atom]<<"\t"<<y[atom]<<"\t"<<z[atom]<<std::endl;
            if(x[atom]==sld::x_seq_coord_array[j] && y[atom]==sld::y_seq_coord_array[j]&&z[atom]==sld::z_seq_coord_array[j] ) {
            //std::cout<<"core at "<<i<<"\t"<<j<<std::endl;
            internal::all_atoms_octant.push_back(i);}
    }
    
     for (int i=0; i<bdry_at;i++){
           	int atom = internal::b_octants[octant][i];
           	if(x[atom]==sld::x_seq_coord_array[j] && y[atom]==sld::y_seq_coord_array[j]&&z[atom]==sld::z_seq_coord_array[j] ) {
           	//std::cout<<"bdry at "<<i<<"\t"<<j<<std::endl;
           	internal::all_atoms_octant.push_back(i);}

           	//std::cout<< " bdry oct and atom " <<octant<<"\t"<<atom<<"\t"<<x[atom]<<"\t"<<y[atom]<<"\t"<<z[atom]<<std::endl;

     }
    
    }
    }
                vmpi::barrier();

    
   suzuki_trotter_parallel_initialized = true;


}

//------------------------------------------------------------------------------
// Integrates a Suzuki trotter step in parallel
//------------------------------------------------------------------------------
void suzuki_trotter_step_parallel(std::vector<double> &x_spin_array,
                      std::vector<double> &y_spin_array,
                      std::vector<double> &z_spin_array,
                      std::vector<int> &type_array){
                      
                      
                      

                    
       double cay_dt=-mp::dt/4.0;//-dt4*consts::gyro - mp::dt contains gamma;
       double dt2_m=0.5*mp::dt_SI*1e12/sld::internal::mp[0].mass.get();
       double dt2=0.5*mp::dt_SI*1e12;
       double f_eta=1.0-0.5*sld::internal::mp[0].damp_lat.get()*mp::dt_SI*1e12;
       double lambda=mp::material[0].alpha;
       double spin_noise=mp::material[0].H_th_sigma*sqrt(sim::temperature);
       double velo_noise=sld::internal::mp[0].F_th_sigma.get()*sqrt(sim::temperature);
       
       
         //vectors for thermal noise spin plus lattice
         std::vector <double> Hx_th(atoms::x_spin_array.size());
      	  std::vector <double> Hy_th(atoms::x_spin_array.size());
      	  std::vector <double> Hz_th(atoms::x_spin_array.size());

         generate (Hx_th.begin(),Hx_th.end(), mtrandom::gaussian);
         generate (Hy_th.begin(),Hy_th.end(), mtrandom::gaussian);
         generate (Hz_th.begin(),Hz_th.end(), mtrandom::gaussian);

         //vectors for thermal forces
         std::vector <double> Fx_th(atoms::x_spin_array.size());
         std::vector <double> Fy_th(atoms::x_spin_array.size());
         std::vector <double> Fz_th(atoms::x_spin_array.size());

         generate (Fx_th.begin(),Fx_th.end(), mtrandom::gaussian);
         generate (Fy_th.begin(),Fy_th.end(), mtrandom::gaussian);
         generate (Fz_th.begin(),Fz_th.end(), mtrandom::gaussian);

   int indx_start, indx_end;
   int number_at=0;
   int atom=0;
	// start first octant loop onwards both core and boundary atoms
   for(int octant = 0; octant < 8; octant++) {
      vmpi::mpi_init_halo_swap();

      //std::cout<<"octant  onwards"<<octant<<"\t"<<nspins<<std::endl;

         
         /*indx_start=internal::all_atoms_octant_start_index[octant];
         indx_end=internal::all_atoms_octant_end_index[octant];
         for(int i=indx_start;i<indx_end;i++){
         
         atom=internal::all_atoms_octant[i];
         
         number_at++;*/
         
         int core_at=internal::c_octants[octant].size();
         int bdry_at=internal::b_octants[octant].size();
         
         for (int i=0; i<core_at;i++){
         atom = internal::c_octants[octant][i];
         sld::internal::fields_array_x[atom]=0.0;
         sld::internal::fields_array_y[atom]=0.0;
         sld::internal::fields_array_z[atom]=0.0;
         
         sld::compute_fields(atom, // first atom for exchange interactions to be calculated
                           atom+1, // last +1 atom to be calculated
                           atoms::neighbour_list_start_index,
                           atoms::neighbour_list_end_index,
                           atoms::type_array, // type for atom
                           atoms::neighbour_list_array, // list of interactions between atoms
                           atoms::x_coord_array,
                           atoms::y_coord_array,
                           atoms::z_coord_array,
                           atoms::x_spin_array,
                           atoms::y_spin_array,
                           atoms::z_spin_array,
                           sld::internal::forces_array_x,
                           sld::internal::forces_array_y,
                           sld::internal::forces_array_z,
                           sld::internal::fields_array_x,
                           sld::internal::fields_array_y,
                           sld::internal::fields_array_z);

         sld::internal::add_spin_noise(atom,
                     atom+1,
                     mp::dt_SI*1e12,
                     atoms::type_array, // type for atom
                     atoms::x_spin_array,
                     atoms::y_spin_array,
                     atoms::z_spin_array,
                     sld::internal::fields_array_x,
                     sld::internal::fields_array_y,
                     sld::internal::fields_array_z,
                     Hx_th, //  vectors for fields
                     Hy_th,
                     Hz_th);
       // std::cout<<"fields therm "<< atom<<"\t"<< sld::internal::fields_array_x[atom]<<"\t"<< sld::internal::fields_array_y[atom]<<"\t"<< sld::internal::fields_array_z[atom]<<std::endl;


         sld::internal::cayley_update(atom,
                     atom+1,
                     cay_dt,
                     atoms::x_spin_array,
                     atoms::y_spin_array,
                     atoms::z_spin_array,
                     sld::internal::fields_array_x,
                     sld::internal::fields_array_y,
                     sld::internal::fields_array_z);
                     
                     
       // std::cout<<"fields final suzuki trotter "<< i<<"\t"<<atom<<"\t"<< sld::internal::fields_array_x[atom]<<"\t"<< sld::internal::fields_array_y[atom]<<"\t"<< sld::internal::fields_array_z[atom]<<"\t"<< atoms::x_coord_array[atom]<<"\t"<< atoms::y_coord_array[atom]<<"\t"<< atoms::z_coord_array[atom]<<"\t"<< atoms::x_spin_array[atom]<<"\t"<< atoms::y_spin_array[atom]<<"\t"<< atoms::z_spin_array[atom]<<std::endl;

         } //end spin for loop
         
                     vmpi::mpi_complete_halo_swap();

         
         for (int i=0; i<bdry_at;i++){
         atom = internal::b_octants[octant][i];
         sld::internal::fields_array_x[atom]=0.0;
         sld::internal::fields_array_y[atom]=0.0;
         sld::internal::fields_array_z[atom]=0.0;
         
         sld::compute_fields(atom, // first atom for exchange interactions to be calculated
                           atom+1, // last +1 atom to be calculated
                           atoms::neighbour_list_start_index,
                           atoms::neighbour_list_end_index,
                           atoms::type_array, // type for atom
                           atoms::neighbour_list_array, // list of interactions between atoms
                           atoms::x_coord_array,
                           atoms::y_coord_array,
                           atoms::z_coord_array,
                           atoms::x_spin_array,
                           atoms::y_spin_array,
                           atoms::z_spin_array,
                           sld::internal::forces_array_x,
                           sld::internal::forces_array_y,
                           sld::internal::forces_array_z,
                           sld::internal::fields_array_x,
                           sld::internal::fields_array_y,
                           sld::internal::fields_array_z);

         sld::internal::add_spin_noise(atom,
                     atom+1,
                     mp::dt_SI*1e12,
                     atoms::type_array, // type for atom
                     atoms::x_spin_array,
                     atoms::y_spin_array,
                     atoms::z_spin_array,
                     sld::internal::fields_array_x,
                     sld::internal::fields_array_y,
                     sld::internal::fields_array_z,
                     Hx_th, //  vectors for fields
                     Hy_th,
                     Hz_th);
       // std::cout<<"fields therm "<< atom<<"\t"<< sld::internal::fields_array_x[atom]<<"\t"<< sld::internal::fields_array_y[atom]<<"\t"<< sld::internal::fields_array_z[atom]<<std::endl;


         sld::internal::cayley_update(atom,
                     atom+1,
                     cay_dt,
                     atoms::x_spin_array,
                     atoms::y_spin_array,
                     atoms::z_spin_array,
                     sld::internal::fields_array_x,
                     sld::internal::fields_array_y,
                     sld::internal::fields_array_z);
                     
                     
       // std::cout<<"fields final suzuki trotter "<< i<<"\t"<<atom<<"\t"<< sld::internal::fields_array_x[atom]<<"\t"<< sld::internal::fields_array_y[atom]<<"\t"<< sld::internal::fields_array_z[atom]<<"\t"<< atoms::x_coord_array[atom]<<"\t"<< atoms::y_coord_array[atom]<<"\t"<< atoms::z_coord_array[atom]<<"\t"<< atoms::x_spin_array[atom]<<"\t"<< atoms::y_spin_array[atom]<<"\t"<< atoms::z_spin_array[atom]<<std::endl;

         } //end spin for loop
         
           //    err::vexit();



      vmpi::barrier();


      }	// end first octant loop onwards
      
      
      //first octant loop return    
      for(int octant = 7; octant >= 0; octant--) {
            vmpi::mpi_init_halo_swap();

          

            /*indx_start=internal::all_atoms_octant_start_index[octant];
            indx_end=internal::all_atoms_octant_end_index[octant];
            for(int i=indx_end-1;i<=indx_start;i--){
            
            atom=internal::all_atoms_octant[i];*/
            int core_at=internal::c_octants[octant].size();
            int bdry_at=internal::b_octants[octant].size();
            
            for (int i=core_at-1;i>=0;i--){
            atom = internal::c_octants[octant][i];
            
            sld::internal::fields_array_x[atom]=0.0;
            sld::internal::fields_array_y[atom]=0.0;
            sld::internal::fields_array_z[atom]=0.0;
            
            sld::compute_fields(atom, // first atom for exchange interactions to be calculated
                              atom+1, // last +1 atom to be calculated
                              atoms::neighbour_list_start_index,
                              atoms::neighbour_list_end_index,
                              atoms::type_array, // type for atom
                              atoms::neighbour_list_array, // list of interactions between atoms
                              atoms::x_coord_array,
                              atoms::y_coord_array,
                              atoms::z_coord_array,
                              atoms::x_spin_array,
                              atoms::y_spin_array,
                              atoms::z_spin_array,
                              sld::internal::forces_array_x,
                              sld::internal::forces_array_y,
                              sld::internal::forces_array_z,
                              sld::internal::fields_array_x,
                              sld::internal::fields_array_y,
                              sld::internal::fields_array_z);

            sld::internal::add_spin_noise(atom,
                        atom+1,
                        mp::dt_SI*1e12,
                        atoms::type_array, // type for atom
                        atoms::x_spin_array,
                        atoms::y_spin_array,
                        atoms::z_spin_array,
                        sld::internal::fields_array_x,
                        sld::internal::fields_array_y,
                        sld::internal::fields_array_z,
                        Hx_th, //  vectors for fields
                        Hy_th,
                        Hz_th);


            sld::internal::cayley_update(atom,
                        atom+1,
                        cay_dt,
                        atoms::x_spin_array,
                        atoms::y_spin_array,
                        atoms::z_spin_array,
                        sld::internal::fields_array_x,
                        sld::internal::fields_array_y,
                        sld::internal::fields_array_z);

            }//end spin loop
            
            vmpi::mpi_complete_halo_swap();

            
            for (int i=bdry_at-1;i>=0;i--){
            atom = internal::b_octants[octant][i];
            
            sld::internal::fields_array_x[atom]=0.0;
            sld::internal::fields_array_y[atom]=0.0;
            sld::internal::fields_array_z[atom]=0.0;
            
            sld::compute_fields(atom, // first atom for exchange interactions to be calculated
                              atom+1, // last +1 atom to be calculated
                              atoms::neighbour_list_start_index,
                              atoms::neighbour_list_end_index,
                              atoms::type_array, // type for atom
                              atoms::neighbour_list_array, // list of interactions between atoms
                              atoms::x_coord_array,
                              atoms::y_coord_array,
                              atoms::z_coord_array,
                              atoms::x_spin_array,
                              atoms::y_spin_array,
                              atoms::z_spin_array,
                              sld::internal::forces_array_x,
                              sld::internal::forces_array_y,
                              sld::internal::forces_array_z,
                              sld::internal::fields_array_x,
                              sld::internal::fields_array_y,
                              sld::internal::fields_array_z);

            sld::internal::add_spin_noise(atom,
                        atom+1,
                        mp::dt_SI*1e12,
                        atoms::type_array, // type for atom
                        atoms::x_spin_array,
                        atoms::y_spin_array,
                        atoms::z_spin_array,
                        sld::internal::fields_array_x,
                        sld::internal::fields_array_y,
                        sld::internal::fields_array_z,
                        Hx_th, //  vectors for fields
                        Hy_th,
                        Hz_th);


            sld::internal::cayley_update(atom,
                        atom+1,
                        cay_dt,
                        atoms::x_spin_array,
                        atoms::y_spin_array,
                        atoms::z_spin_array,
                        sld::internal::fields_array_x,
                        sld::internal::fields_array_y,
                        sld::internal::fields_array_z);

            }//end spin loop
             vmpi::barrier();


         }    // end first octant loop return




       
       for(int octant = 0; octant < 8; octant++) {

          //std::cout<<"octant  onwards"<<octant<<"\t"<<nspins<<std::endl;

             vmpi::mpi_init_halo_swap();
             int core_at=internal::c_octants[octant].size();
             int bdry_at=internal::b_octants[octant].size();


            /*
             indx_start=internal::all_atoms_octant_start_index[octant];
             indx_end=internal::all_atoms_octant_end_index[octant];
             for(int i=indx_start;i<indx_end;i++){
             
             atom=internal::all_atoms_octant[i];*/
             
             for (int i=0; i<core_at;i++){
             atom = internal::c_octants[octant][i];
             
             sld::internal::fields_array_x[atom]=0.0;
             sld::internal::fields_array_y[atom]=0.0;
             sld::internal::fields_array_z[atom]=0.0;
             
             sld::compute_fields(atom, // first atom for exchange interactions to be calculated
                               atom+1, // last +1 atom to be calculated
                               atoms::neighbour_list_start_index,
                               atoms::neighbour_list_end_index,
                               atoms::type_array, // type for atom
                               atoms::neighbour_list_array, // list of interactions between atoms
                               atoms::x_coord_array,
                               atoms::y_coord_array,
                               atoms::z_coord_array,
                               atoms::x_spin_array,
                               atoms::y_spin_array,
                               atoms::z_spin_array,
                               sld::internal::forces_array_x,
                               sld::internal::forces_array_y,
                               sld::internal::forces_array_z,
                               sld::internal::fields_array_x,
                               sld::internal::fields_array_y,
                               sld::internal::fields_array_z);

             sld::internal::add_spin_noise(atom,
                         atom+1,
                         mp::dt_SI*1e12,
                         atoms::type_array, // type for atom
                         atoms::x_spin_array,
                         atoms::y_spin_array,
                         atoms::z_spin_array,
                         sld::internal::fields_array_x,
                         sld::internal::fields_array_y,
                         sld::internal::fields_array_z,
                         Hx_th, //  vectors for fields
                         Hy_th,
                         Hz_th);
           // std::cout<<"fields therm "<< atom<<"\t"<< sld::internal::fields_array_x[atom]<<"\t"<< sld::internal::fields_array_y[atom]<<"\t"<< sld::internal::fields_array_z[atom]<<std::endl;


             sld::internal::cayley_update(atom,
                         atom+1,
                         cay_dt,
                         atoms::x_spin_array,
                         atoms::y_spin_array,
                         atoms::z_spin_array,
                         sld::internal::fields_array_x,
                         sld::internal::fields_array_y,
                         sld::internal::fields_array_z);
                         
                         
           // std::cout<<"fields final suzuki trotter "<< atom<<"\t"<< sld::internal::fields_array_x[atom]<<"\t"<< sld::internal::fields_array_y[atom]<<"\t"<< sld::internal::fields_array_z[atom]<<"\t"<< atoms::x_coord_array[atom]<<"\t"<< atoms::y_coord_array[atom]<<"\t"<< atoms::z_coord_array[atom]<<std::endl;

             } //end spin for loop
        
        vmpi::mpi_complete_halo_swap();

             
        for (int i=0; i<bdry_at;i++){
        atom = internal::b_octants[octant][i];
        
        sld::internal::fields_array_x[atom]=0.0;
        sld::internal::fields_array_y[atom]=0.0;
        sld::internal::fields_array_z[atom]=0.0;
        
        sld::compute_fields(atom, // first atom for exchange interactions to be calculated
                          atom+1, // last +1 atom to be calculated
                          atoms::neighbour_list_start_index,
                          atoms::neighbour_list_end_index,
                          atoms::type_array, // type for atom
                          atoms::neighbour_list_array, // list of interactions between atoms
                          atoms::x_coord_array,
                          atoms::y_coord_array,
                          atoms::z_coord_array,
                          atoms::x_spin_array,
                          atoms::y_spin_array,
                          atoms::z_spin_array,
                          sld::internal::forces_array_x,
                          sld::internal::forces_array_y,
                          sld::internal::forces_array_z,
                          sld::internal::fields_array_x,
                          sld::internal::fields_array_y,
                          sld::internal::fields_array_z);

        sld::internal::add_spin_noise(atom,
                    atom+1,
                    mp::dt_SI*1e12,
                    atoms::type_array, // type for atom
                    atoms::x_spin_array,
                    atoms::y_spin_array,
                    atoms::z_spin_array,
                    sld::internal::fields_array_x,
                    sld::internal::fields_array_y,
                    sld::internal::fields_array_z,
                    Hx_th, //  vectors for fields
                    Hy_th,
                    Hz_th);
      // std::cout<<"fields therm "<< atom<<"\t"<< sld::internal::fields_array_x[atom]<<"\t"<< sld::internal::fields_array_y[atom]<<"\t"<< sld::internal::fields_array_z[atom]<<std::endl;


        sld::internal::cayley_update(atom,
                    atom+1,
                    cay_dt,
                    atoms::x_spin_array,
                    atoms::y_spin_array,
                    atoms::z_spin_array,
                    sld::internal::fields_array_x,
                    sld::internal::fields_array_y,
                    sld::internal::fields_array_z);
                    
                    
      // std::cout<<"fields final suzuki trotter "<< atom<<"\t"<< sld::internal::fields_array_x[atom]<<"\t"<< sld::internal::fields_array_y[atom]<<"\t"<< sld::internal::fields_array_z[atom]<<"\t"<< atoms::x_coord_array[atom]<<"\t"<< atoms::y_coord_array[atom]<<"\t"<< atoms::z_coord_array[atom]<<std::endl;

        } //end spin for loop

    vmpi::barrier();


          }	// end second octant loop onwards

          
          
          //first octant loop return    
          for(int octant = 7; octant >= 0; octant--) {
          
                 vmpi::mpi_init_halo_swap();
                 int core_at=internal::c_octants[octant].size();
                 int bdry_at=internal::b_octants[octant].size();

              
                /*
                indx_start=internal::all_atoms_octant_start_index[octant];
                indx_end=internal::all_atoms_octant_end_index[octant];
                for(int i=indx_end-1;i<=indx_start;i--){
                
                atom=internal::all_atoms_octant[i];*/
                
                for (int i=core_at-1;i>=0;i--){
                atom = internal::c_octants[octant][i];
                
                
                sld::internal::fields_array_x[atom]=0.0;
                sld::internal::fields_array_y[atom]=0.0;
                sld::internal::fields_array_z[atom]=0.0;
                sld::compute_fields(atom, // first atom for exchange interactions to be calculated
                                  atom+1, // last +1 atom to be calculated
                                  atoms::neighbour_list_start_index,
                                  atoms::neighbour_list_end_index,
                                  atoms::type_array, // type for atom
                                  atoms::neighbour_list_array, // list of interactions between atoms
                                  atoms::x_coord_array,
                                  atoms::y_coord_array,
                                  atoms::z_coord_array,
                                  atoms::x_spin_array,
                                  atoms::y_spin_array,
                                  atoms::z_spin_array,
                                  sld::internal::forces_array_x,
                                  sld::internal::forces_array_y,
                                  sld::internal::forces_array_z,
                                  sld::internal::fields_array_x,
                                  sld::internal::fields_array_y,
                                  sld::internal::fields_array_z);

                sld::internal::add_spin_noise(atom,
                            atom+1,
                            mp::dt_SI*1e12,
                            atoms::type_array, // type for atom
                            atoms::x_spin_array,
                            atoms::y_spin_array,
                            atoms::z_spin_array,
                            sld::internal::fields_array_x,
                            sld::internal::fields_array_y,
                            sld::internal::fields_array_z,
                            Hx_th, //  vectors for fields
                            Hy_th,
                            Hz_th);


                sld::internal::cayley_update(atom,
                            atom+1,
                            cay_dt,
                            atoms::x_spin_array,
                            atoms::y_spin_array,
                            atoms::z_spin_array,
                            sld::internal::fields_array_x,
                            sld::internal::fields_array_y,
                            sld::internal::fields_array_z);

                }//end spin loop
                
                vmpi::mpi_complete_halo_swap();

                for (int i=bdry_at-1;i>=0;i--){
                atom = internal::b_octants[octant][i];
                
                
                sld::internal::fields_array_x[atom]=0.0;
                sld::internal::fields_array_y[atom]=0.0;
                sld::internal::fields_array_z[atom]=0.0;
                sld::compute_fields(atom, // first atom for exchange interactions to be calculated
                                  atom+1, // last +1 atom to be calculated
                                  atoms::neighbour_list_start_index,
                                  atoms::neighbour_list_end_index,
                                  atoms::type_array, // type for atom
                                  atoms::neighbour_list_array, // list of interactions between atoms
                                  atoms::x_coord_array,
                                  atoms::y_coord_array,
                                  atoms::z_coord_array,
                                  atoms::x_spin_array,
                                  atoms::y_spin_array,
                                  atoms::z_spin_array,
                                  sld::internal::forces_array_x,
                                  sld::internal::forces_array_y,
                                  sld::internal::forces_array_z,
                                  sld::internal::fields_array_x,
                                  sld::internal::fields_array_y,
                                  sld::internal::fields_array_z);

                sld::internal::add_spin_noise(atom,
                            atom+1,
                            mp::dt_SI*1e12,
                            atoms::type_array, // type for atom
                            atoms::x_spin_array,
                            atoms::y_spin_array,
                            atoms::z_spin_array,
                            sld::internal::fields_array_x,
                            sld::internal::fields_array_y,
                            sld::internal::fields_array_z,
                            Hx_th, //  vectors for fields
                            Hy_th,
                            Hz_th);


                sld::internal::cayley_update(atom,
                            atom+1,
                            cay_dt,
                            atoms::x_spin_array,
                            atoms::y_spin_array,
                            atoms::z_spin_array,
                            sld::internal::fields_array_x,
                            sld::internal::fields_array_y,
                            sld::internal::fields_array_z);

                }//end spin loop

             vmpi::barrier();

             }      //end second octant loop return 
             
             

    vmpi::barrier();
           
 // Swap timers compute -> wait
 vmpi::TotalComputeTime+=vmpi::SwapTimer(vmpi::ComputeTime, vmpi::WaitTime);

 // Wait for other processors
 vmpi::barrier();

 // Swap timers wait -> compute
 vmpi::TotalWaitTime += vmpi::SwapTimer(vmpi::WaitTime, vmpi::ComputeTime);
    
   
   return;

}

} // End of namespace sld
#endif
