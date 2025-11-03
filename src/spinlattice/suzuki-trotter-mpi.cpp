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
#include "atoms.hpp"
#include "create.hpp"
#include "errors.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "material.hpp"
#include "vmpi.hpp"

// Internal header
#include "internal.hpp"







namespace sld{



void suzuki_trotter_parallel_init(std::vector<double> &x, // atomic coordinates
                      std::vector<double> &y,
                      std::vector<double> &z,
                      double min_dim[3], // minimum dimensions on local processor
                      double max_dim[3]){ // maximum dimensions on local processor


   // Convenient shorthands
   int catoms = vmpi::num_core_atoms;
   int batoms = vmpi::num_bdry_atoms;



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
   for(int i=0; i< 8; i++) {num_atoms_in_octants += internal::c_octants[i].size();}

   if(num_atoms_in_octants != catoms){
      std::cerr << "Programmer error: missing atoms in core octants in suzuki-trotter initialisation" << std::endl;
      err::vexit();
   }
   // boundary atoms
   num_atoms_in_octants = 0;
   for(int i=0; i< 8; i++) {
   num_atoms_in_octants += internal::b_octants[i].size();}

   if(num_atoms_in_octants != batoms){
      std::cout<< "num_atoms "<<num_atoms_in_octants<<"\t"<<batoms<<std::endl;
      std::cerr << "Programmer error: missing atoms in boundary octants in suzuki-trotter initialisation" << std::endl;
      err::vexit();
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
       double dt2=0.5*mp::dt_SI*1e12;



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

   //int indx_start, indx_end;
   //int number_at=0;
   int atom=0;


/*for(int atom=0;atom<=atoms::num_atoms-1;atom++){

   if (atom==0) {
   std::cout<<std::setprecision(15)<<std::endl;

   std::cout <<" b i=atom="<<atom<<"\txyz "<<atoms::x_coord_array[atom]<<"\t"<<atoms::y_coord_array[atom]<<"\t"<<atoms::z_coord_array[atom]<<std::endl;
   std::cout <<" b i=atom="<<atom<<"\tvel "<<atoms::x_velo_array[atom]<<"\t"<<atoms::y_velo_array[atom]<<"\t"<<atoms::z_velo_array[atom]<<std::endl;


   }}*/

   //for velocity and position update
   const int pre_comm_si = 0;
   const int pre_comm_ei = vmpi::num_core_atoms;
   const int post_comm_si = vmpi::num_core_atoms;
   const int post_comm_ei = vmpi::num_core_atoms+vmpi::num_bdry_atoms;




      //----------------------------------------
      // Initialise storage array (all)
      //----------------------------------------


     for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
 sld::internal::x_coord_storage_array[atom]=atoms::x_coord_array[atom];
 sld::internal::y_coord_storage_array[atom]=atoms::y_coord_array[atom];
 sld::internal::z_coord_storage_array[atom]=atoms::z_coord_array[atom];
 }

      vmpi::barrier();
	// start first octant loop onwards both core and boundary atoms
	//all fields set to 0
	std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
    std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
    std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);


    for(int octant = 0; octant < 8; octant++) {
      vmpi::mpi_init_halo_swap();


         int core_at=internal::c_octants[octant].size();
         int bdry_at=internal::b_octants[octant].size();

         for (int i=0; i<core_at;i++){
         atom = internal::c_octants[octant][i];


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



         } //end spin for loop

         vmpi::mpi_complete_halo_swap();
         vmpi::barrier();


         for (int i=0; i<bdry_at;i++){
         atom = internal::b_octants[octant][i];



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



         } //end spin for loop




      vmpi::barrier();


      }	// end first octant loop onwards


      //all fields set to 0
  	  std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
      std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
      std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);


      //first octant loop return
      for(int octant = 7; octant >= 0; octant--) {
            vmpi::mpi_init_halo_swap();


            int core_at=internal::c_octants[octant].size();
            int bdry_at=internal::b_octants[octant].size();

            for (int i=core_at-1;i>=0;i--){
            atom = internal::c_octants[octant][i];


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
            vmpi::barrier();


            for (int i=bdry_at-1;i>=0;i--){
            atom = internal::b_octants[octant][i];

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


     //vrv update
     //CORE atoms first
     //all forces are set to 0 for computation
     std::fill(sld::internal::forces_array_x.begin(), sld::internal::forces_array_x.end(), 0.0);
     std::fill(sld::internal::forces_array_y.begin(), sld::internal::forces_array_y.end(), 0.0);
     std::fill(sld::internal::forces_array_z.begin(), sld::internal::forces_array_z.end(), 0.0);

     vmpi::mpi_init_halo_swap_coords();




     sld::compute_fields(pre_comm_si, // first atom for exchange interactions to be calculated
                       pre_comm_ei, // last +1 atom to be calculated
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




     sld::compute_forces(pre_comm_si, // first atom for exchange interactions to be calculated
                       pre_comm_ei, // last +1 atom to be calculated
                       atoms::neighbour_list_start_index,
                       atoms::neighbour_list_end_index,
                       atoms::type_array, // type for atom
                       atoms::neighbour_list_array, // list of interactions between atoms
                       sld::internal::x0_coord_array, // list of isotropic exchange constants
                       sld::internal::y0_coord_array, // list of vectorial exchange constants
                       sld::internal::z0_coord_array, // list of tensorial exchange constants
                       atoms::x_coord_array,
                       atoms::y_coord_array,
                       atoms::z_coord_array,
                       sld::internal::forces_array_x,
                       sld::internal::forces_array_y,
                       sld::internal::forces_array_z,
                       sld::internal::potential_eng);

     //update position, Velocity
     		for(int atom=pre_comm_si;atom<pre_comm_ei;atom++){



     		     const unsigned int imat = atoms::type_array[atom];
 		         double f_eta=1.0-0.5*sld::internal::mp[imat].damp_lat.get()*mp::dt_SI*1e12;
                 double velo_noise=sld::internal::mp[imat].F_th_sigma.get()*sqrt(sim::temperature);
                 double dt2_m=0.5*mp::dt_SI*1e12/sld::internal::mp[imat].mass.get();

                //if during equilibration:
                if (sim::time < sim::equilibration_time) {
                      f_eta=1.0-0.5*sld::internal::mp[imat].eq_damp_lat.get()*mp::dt_SI*1e12;
                      velo_noise=sld::internal::mp[imat].F_th_sigma_eq.get()*sqrt(sim::temperature);
                }

                 atoms::x_velo_array[atom] =  f_eta*atoms::x_velo_array[atom] + dt2_m * sld::internal::forces_array_x[atom]+dt2*velo_noise*Fx_th[atom];
                 atoms::y_velo_array[atom] =  f_eta*atoms::y_velo_array[atom] + dt2_m * sld::internal::forces_array_y[atom]+dt2*velo_noise*Fy_th[atom];
                 atoms::z_velo_array[atom] =  f_eta*atoms::z_velo_array[atom] + dt2_m * sld::internal::forces_array_z[atom]+dt2*velo_noise*Fz_th[atom];

                 sld::internal::x_coord_storage_array[atom] +=  mp::dt_SI*1e12 * atoms::x_velo_array[atom];
                 sld::internal::y_coord_storage_array[atom] +=  mp::dt_SI*1e12 * atoms::y_velo_array[atom];
                 sld::internal::z_coord_storage_array[atom] +=  mp::dt_SI*1e12 * atoms::z_velo_array[atom];




   }


   vmpi::mpi_complete_halo_swap_coords();
   vmpi::barrier();

   sld::compute_fields(post_comm_si, // first atom for exchange interactions to be calculated
                      post_comm_ei, // last +1 atom to be calculated
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

   sld::compute_forces(post_comm_si, // first atom for exchange interactions to be calculated
                     post_comm_ei, // last +1 atom to be calculated
                     atoms::neighbour_list_start_index,
                     atoms::neighbour_list_end_index,
                     atoms::type_array, // type for atom
                     atoms::neighbour_list_array, // list of interactions between atoms
                     sld::internal::x0_coord_array, // list of isotropic exchange constants
                     sld::internal::y0_coord_array, // list of vectorial exchange constants
                     sld::internal::z0_coord_array, // list of tensorial exchange constants
                     atoms::x_coord_array,
                     atoms::y_coord_array,
                     atoms::z_coord_array,
                     sld::internal::forces_array_x,
                     sld::internal::forces_array_y,
                     sld::internal::forces_array_z,
                     sld::internal::potential_eng);







      //update position, Velocity for boundary atoms
      		for(int atom=post_comm_si;atom<post_comm_ei;atom++){


     		     const unsigned int imat = atoms::type_array[atom];
 		         double f_eta=1.0-0.5*sld::internal::mp[imat].damp_lat.get()*mp::dt_SI*1e12;
                 double velo_noise=sld::internal::mp[imat].F_th_sigma.get()*sqrt(sim::temperature);
                 double dt2_m=0.5*mp::dt_SI*1e12/sld::internal::mp[imat].mass.get();

                   //if during equilibration:
                   if (sim::time < sim::equilibration_time) {
                         f_eta=1.0-0.5*sld::internal::mp[imat].eq_damp_lat.get()*mp::dt_SI*1e12;
                         velo_noise=sld::internal::mp[imat].F_th_sigma_eq.get()*sqrt(sim::temperature);
                   }

                   atoms::x_velo_array[atom] =  f_eta*atoms::x_velo_array[atom] + dt2_m * sld::internal::forces_array_x[atom]+dt2*velo_noise*Fx_th[atom];
                   atoms::y_velo_array[atom] =  f_eta*atoms::y_velo_array[atom] + dt2_m * sld::internal::forces_array_y[atom]+dt2*velo_noise*Fy_th[atom];
                   atoms::z_velo_array[atom] =  f_eta*atoms::z_velo_array[atom] + dt2_m * sld::internal::forces_array_z[atom]+dt2*velo_noise*Fz_th[atom];


                  sld::internal::x_coord_storage_array[atom] +=  mp::dt_SI*1e12 * atoms::x_velo_array[atom];
                  sld::internal::y_coord_storage_array[atom] +=  mp::dt_SI*1e12 * atoms::y_velo_array[atom];
                  sld::internal::z_coord_storage_array[atom] +=  mp::dt_SI*1e12 * atoms::z_velo_array[atom];


             }

             vmpi::barrier();


            //----------------------------------------
    		// Copy new positions array (all)
    		//----------------------------------------
    		for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
    			atoms::x_coord_array[atom]=sld::internal::x_coord_storage_array[atom];
    			atoms::y_coord_array[atom]=sld::internal::y_coord_storage_array[atom];
    			atoms::z_coord_array[atom]=sld::internal::z_coord_storage_array[atom];
    		}

            //reset forces to 0 for v integration
            vmpi::mpi_init_halo_swap_coords();



            std::fill(sld::internal::forces_array_x.begin(), sld::internal::forces_array_x.end(), 0.0);
            std::fill(sld::internal::forces_array_y.begin(), sld::internal::forces_array_y.end(), 0.0);
            std::fill(sld::internal::forces_array_z.begin(), sld::internal::forces_array_z.end(), 0.0);



            sld::compute_fields(pre_comm_si, // first atom for exchange interactions to be calculated
                               pre_comm_ei, // last +1 atom to be calculated
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



            sld::compute_forces(pre_comm_si, // first atom for exchange interactions to be calculated
                               pre_comm_ei, // last +1 atom to be calculated
                               atoms::neighbour_list_start_index,
                               atoms::neighbour_list_end_index,
                               atoms::type_array, // type for atom
                               atoms::neighbour_list_array, // list of interactions between atoms
                               sld::internal::x0_coord_array, // list of isotropic exchange constants
                               sld::internal::y0_coord_array, // list of vectorial exchange constants
                               sld::internal::z0_coord_array, // list of tensorial exchange constants
                               atoms::x_coord_array,
                               atoms::y_coord_array,
                               atoms::z_coord_array,
                               sld::internal::forces_array_x,
                               sld::internal::forces_array_y,
                               sld::internal::forces_array_z,
                               sld::internal::potential_eng);



           for(int atom=pre_comm_si;atom<pre_comm_ei;atom++){


     		     const unsigned int imat = atoms::type_array[atom];
 		         double f_eta=1.0-0.5*sld::internal::mp[imat].damp_lat.get()*mp::dt_SI*1e12;
                 double velo_noise=sld::internal::mp[imat].F_th_sigma.get()*sqrt(sim::temperature);
                 double dt2_m=0.5*mp::dt_SI*1e12/sld::internal::mp[imat].mass.get();

                 //if during equilibration:
                 if (sim::time < sim::equilibration_time) {
                       f_eta=1.0-0.5*sld::internal::mp[imat].eq_damp_lat.get()*mp::dt_SI*1e12;
                       velo_noise=sld::internal::mp[imat].F_th_sigma_eq.get()*sqrt(sim::temperature);
                 }

              atoms::x_velo_array[atom] =  f_eta*atoms::x_velo_array[atom] + dt2_m * sld::internal::forces_array_x[atom]+dt2*velo_noise*Fx_th[atom];
              atoms::y_velo_array[atom] =  f_eta*atoms::y_velo_array[atom] + dt2_m * sld::internal::forces_array_y[atom]+dt2*velo_noise*Fy_th[atom];
              atoms::z_velo_array[atom] =  f_eta*atoms::z_velo_array[atom] + dt2_m * sld::internal::forces_array_z[atom]+dt2*velo_noise*Fz_th[atom];


           }


    vmpi::mpi_complete_halo_swap_coords();
    vmpi::barrier();


           sld::compute_fields(post_comm_si, // first atom for exchange interactions to be calculated
                               post_comm_ei, // last +1 atom to be calculated
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

            sld::compute_forces(post_comm_si, // first atom for exchange interactions to be calculated
                              post_comm_ei, // last +1 atom to be calculated
                              atoms::neighbour_list_start_index,
                              atoms::neighbour_list_end_index,
                              atoms::type_array, // type for atom
                              atoms::neighbour_list_array, // list of interactions between atoms
                              sld::internal::x0_coord_array, // list of isotropic exchange constants
                              sld::internal::y0_coord_array, // list of vectorial exchange constants
                              sld::internal::z0_coord_array, // list of tensorial exchange constants
                              atoms::x_coord_array,
                              atoms::y_coord_array,
                              atoms::z_coord_array,
                              sld::internal::forces_array_x,
                              sld::internal::forces_array_y,
                              sld::internal::forces_array_z,
                              sld::internal::potential_eng);


        //

           for(int atom=post_comm_si;atom<post_comm_ei;atom++){


     		     const unsigned int imat = atoms::type_array[atom];
 		         double f_eta=1.0-0.5*sld::internal::mp[imat].damp_lat.get()*mp::dt_SI*1e12;
                 double velo_noise=sld::internal::mp[imat].F_th_sigma.get()*sqrt(sim::temperature);
                 double dt2_m=0.5*mp::dt_SI*1e12/sld::internal::mp[imat].mass.get();

                 //if during equilibration:
                 if (sim::time < sim::equilibration_time) {
                       f_eta=1.0-0.5*sld::internal::mp[imat].eq_damp_lat.get()*mp::dt_SI*1e12;
                       velo_noise=sld::internal::mp[imat].F_th_sigma_eq.get()*sqrt(sim::temperature);
                       }

              atoms::x_velo_array[atom] =  f_eta*atoms::x_velo_array[atom] + dt2_m * sld::internal::forces_array_x[atom]+dt2*velo_noise*Fx_th[atom];
              atoms::y_velo_array[atom] =  f_eta*atoms::y_velo_array[atom] + dt2_m * sld::internal::forces_array_y[atom]+dt2*velo_noise*Fy_th[atom];
              atoms::z_velo_array[atom] =  f_eta*atoms::z_velo_array[atom] + dt2_m * sld::internal::forces_array_z[atom]+dt2*velo_noise*Fz_th[atom];


            }


       		vmpi::barrier();

      //spin update
      //all fields set to 0
  	  std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
      std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
      std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);


       for(int octant = 0; octant < 8; octant++) {

          //std::cout<<"octant  onwards"<<octant<<"\t"<<nspins<<std::endl;

             vmpi::mpi_init_halo_swap();
             int core_at=internal::c_octants[octant].size();
             int bdry_at=internal::b_octants[octant].size();


             for (int i=0; i<core_at;i++){
             atom = internal::c_octants[octant][i];


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



             } //end spin for loop

        vmpi::mpi_complete_halo_swap();
        vmpi::barrier();



        for (int i=0; i<bdry_at;i++){
        atom = internal::b_octants[octant][i];


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



        } //end spin for loop

    vmpi::barrier();


          }	// end second octant loop onwards


  	  std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
      std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
      std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);


          //first octant loop return
          for(int octant = 7; octant >= 0; octant--) {

                 vmpi::mpi_init_halo_swap();
                 int core_at=internal::c_octants[octant].size();
                 int bdry_at=internal::b_octants[octant].size();



                for (int i=core_at-1;i>=0;i--){
                atom = internal::c_octants[octant][i];

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
                vmpi::barrier();


                for (int i=bdry_at-1;i>=0;i--){
                atom = internal::b_octants[octant][i];


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



/*for(int atom=0;atom<=atoms::num_atoms-1;atom++){

   if ( (atom==0)) {
   std::cout<<std::setprecision(15)<<std::endl;

   std::cout <<" end i=atom="<<atom<<"\txyz "<<atoms::x_coord_array[atom]<<"\t"<<atoms::y_coord_array[atom]<<"\t"<<atoms::z_coord_array[atom]<<std::endl;
   std::cout <<" e i=atom="<<atom<<"\tvel "<<atoms::x_velo_array[atom]<<"\t"<<atoms::y_velo_array[atom]<<"\t"<<atoms::z_velo_array[atom]<<std::endl;


   }}*/

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
