//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2023. All rights reserved.
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
#include "atoms.hpp"
#include "create.hpp"
#include "material.hpp"
#include "internal.hpp"









namespace sim{


void STDspin_parallel_init(std::vector<double> &x, // atomic coordinates
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

  /* if (widthx<=sld::internal::r_cut_fields ||widthy<=sld::internal::r_cut_fields || widthy<=sld::internal::r_cut_fields ){
     std::cerr << "Error: System size needs to be increased so octants won't interact" << std::endl;
     std::cerr<<"Width of octants is "<<widthx<<"\t"<<widthy<<"\t"<<widthz<< " While interaction cutoff is "<<sld::internal::r_cut_fields<<std::endl;
     err::vexit();
     }*/

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
                  sim::c_octants[octant_num].push_back(i);
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
                  sim::b_octants[octant_num].push_back(i);

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
   for(int i=0; i< 8; i++) {num_atoms_in_octants += sim::c_octants[i].size();}

   if(num_atoms_in_octants != catoms){
      std::cerr << "Programmer error: missing atoms in core octants in suzuki-trotter initialisation" << std::endl;
      err::vexit();
   }
   // boundary atoms
   num_atoms_in_octants = 0;
   for(int i=0; i< 8; i++) {
   num_atoms_in_octants += sim::b_octants[i].size();}

   if(num_atoms_in_octants != batoms){
      std::cerr << "Programmer error: missing atoms in boundary octants in suzuki-trotter initialisation" << std::endl;
      err::vexit();
   }


    vmpi::barrier();

    
   STDspin_parallel_initialized = true;


}

//------------------------------------------------------------------------------
// Integrates a Suzuki trotter step in parallel
//------------------------------------------------------------------------------
void STDspin_step_parallel(std::vector<double> &x_spin_array,
                      std::vector<double> &y_spin_array,
                      std::vector<double> &z_spin_array,
                      std::vector<int> &type_array){
                      
                      
                      

                    
       double cay_dt=-mp::dt/4.0;//-dt4*consts::gyro - mp::dt contains gamma;
       double dt2=0.5*mp::dt_SI*1e12;
       double lambda;
       double spin_noise;
       

       
         //vectors for thermal noise spin plus lattice
         std::vector <double> Hx_th(atoms::x_spin_array.size());
      	  std::vector <double> Hy_th(atoms::x_spin_array.size());
      	  std::vector <double> Hz_th(atoms::x_spin_array.size());

         generate (Hx_th.begin(),Hx_th.end(), mtrandom::gaussian);
         generate (Hy_th.begin(),Hy_th.end(), mtrandom::gaussian);
         generate (Hz_th.begin(),Hz_th.end(), mtrandom::gaussian);

         

   int indx_start, indx_end;
   int number_at=0;
   int atom=0;
   
   //for velocity and position update
   const int pre_comm_si = 0;
   const int pre_comm_ei = vmpi::num_core_atoms;
   const int post_comm_si = vmpi::num_core_atoms;
   const int post_comm_ei = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
   

	
    for(int octant = 0; octant < 8; octant++) {
      vmpi::mpi_init_halo_swap();

         
         int core_at=sim::c_octants[octant].size();
         int bdry_at=sim::b_octants[octant].size();
         
         for (int i=0; i<core_at;i++){
         atom = sim::c_octants[octant][i];
         const unsigned int imat = atoms::type_array[atom];
         lambda=mp::material[imat].alpha;
         spin_noise=mp::material[imat].H_th_sigma*sqrt(sim::temperature);
         
         
         calculate_spin_fields(atom, atom+1);



         add_spin_noise(atom,
                     atom+1,
                     mp::dt_SI*1e12,
                     atoms::type_array, // type for atom
                     atoms::x_spin_array,
                     atoms::y_spin_array,
                     atoms::z_spin_array,
                     atoms::x_total_spin_field_array,
                     atoms::y_total_spin_field_array,
                     atoms::z_total_spin_field_array,
                     Hx_th, //  vectors for fields
                     Hy_th,
                     Hz_th);
                     



        cayley_update(atom,
                     atom+1,
                     cay_dt,
                     atoms::x_spin_array,
                     atoms::y_spin_array,
                     atoms::z_spin_array,
                     atoms::x_total_spin_field_array,
                     atoms::y_total_spin_field_array,
                     atoms::z_total_spin_field_array);
          

         } //end spin for loop
         
         vmpi::mpi_complete_halo_swap();
         vmpi::barrier();

         
         for (int i=0; i<bdry_at;i++){
         atom = sim::b_octants[octant][i];
         const unsigned int imat = atoms::type_array[atom];
         lambda=mp::material[imat].alpha;
         spin_noise=mp::material[imat].H_th_sigma*sqrt(sim::temperature);
         
         calculate_spin_fields(atom, atom+1);

         add_spin_noise(atom,
                     atom+1,
                     mp::dt_SI*1e12,
                     atoms::type_array, // type for atom
                     atoms::x_spin_array,
                     atoms::y_spin_array,
                     atoms::z_spin_array,
                     atoms::x_total_spin_field_array,
                     atoms::y_total_spin_field_array,
                     atoms::z_total_spin_field_array,
                     Hx_th, //  vectors for fields
                     Hy_th,
                     Hz_th);
                     



        cayley_update(atom,
                     atom+1,
                     cay_dt,
                     atoms::x_spin_array,
                     atoms::y_spin_array,
                     atoms::z_spin_array,
                     atoms::x_total_spin_field_array,
                     atoms::y_total_spin_field_array,
                     atoms::z_total_spin_field_array);
          
                     

         } //end spin for loop
         



      vmpi::barrier();


      }	// end first octant loop onwards
      
  	
      
      //first octant loop return    
      for(int octant = 7; octant >= 0; octant--) {
            vmpi::mpi_init_halo_swap();

          
            int core_at=sim::c_octants[octant].size();
            int bdry_at=sim::b_octants[octant].size();
            
            for (int i=core_at-1;i>=0;i--){
            atom = sim::c_octants[octant][i];
            const unsigned int imat = atoms::type_array[atom];
            lambda=mp::material[imat].alpha;
            spin_noise=mp::material[imat].H_th_sigma*sqrt(sim::temperature);
            
            
            calculate_spin_fields(atom, atom+1);

            add_spin_noise(atom,
                        atom+1,
                        mp::dt_SI*1e12,
                        atoms::type_array, // type for atom
                        atoms::x_spin_array,
                        atoms::y_spin_array,
                        atoms::z_spin_array,
                        atoms::x_total_spin_field_array,
                        atoms::y_total_spin_field_array,
                        atoms::z_total_spin_field_array,
                        Hx_th, //  vectors for fields
                        Hy_th,
                        Hz_th);
                        



           cayley_update(atom,
                        atom+1,
                        cay_dt,
                        atoms::x_spin_array,
                        atoms::y_spin_array,
                        atoms::z_spin_array,
                        atoms::x_total_spin_field_array,
                        atoms::y_total_spin_field_array,
                        atoms::z_total_spin_field_array);
             
          

            }//end spin loop
            
            vmpi::mpi_complete_halo_swap();
            vmpi::barrier();

            
            for (int i=bdry_at-1;i>=0;i--){
            atom = sim::b_octants[octant][i];
            
            const unsigned int imat = atoms::type_array[atom];
            lambda=mp::material[imat].alpha;
            spin_noise=mp::material[imat].H_th_sigma*sqrt(sim::temperature);
            
            calculate_spin_fields(atom, atom+1);

            add_spin_noise(atom,
                        atom+1,
                        mp::dt_SI*1e12,
                        atoms::type_array, // type for atom
                        atoms::x_spin_array,
                        atoms::y_spin_array,
                        atoms::z_spin_array,
                        atoms::x_total_spin_field_array,
                        atoms::y_total_spin_field_array,
                        atoms::z_total_spin_field_array,
                        Hx_th, //  vectors for fields
                        Hy_th,
                        Hz_th);
                        



           cayley_update(atom,
                        atom+1,
                        cay_dt,
                        atoms::x_spin_array,
                        atoms::y_spin_array,
                        atoms::z_spin_array,
                        atoms::x_total_spin_field_array,
                        atoms::y_total_spin_field_array,
                        atoms::z_total_spin_field_array);
             

            }//end spin loop
             vmpi::barrier();


         }    // end first octant loop return


  	
       
       for(int octant = 0; octant < 8; octant++) {

          //std::cout<<"octant  onwards"<<octant<<"\t"<<nspins<<std::endl;

             vmpi::mpi_init_halo_swap();
             int core_at=sim::c_octants[octant].size();
             int bdry_at=sim::b_octants[octant].size();
           
             
             for (int i=0; i<core_at;i++){
             atom = sim::c_octants[octant][i];
             const unsigned int imat = atoms::type_array[atom];
             lambda=mp::material[imat].alpha;
              spin_noise=mp::material[imat].H_th_sigma*sqrt(sim::temperature);
             
             
             calculate_spin_fields(atom, atom+1);
         	  

             add_spin_noise(atom,
                         atom+1,
                         mp::dt_SI*1e12,
                         atoms::type_array, // type for atom
                         atoms::x_spin_array,
                         atoms::y_spin_array,
                         atoms::z_spin_array,
                         atoms::x_total_spin_field_array,
                         atoms::y_total_spin_field_array,
                         atoms::z_total_spin_field_array,
                         Hx_th, //  vectors for fields
                         Hy_th,
                         Hz_th);
                         



            cayley_update(atom,
                         atom+1,
                         cay_dt,
                         atoms::x_spin_array,
                         atoms::y_spin_array,
                         atoms::z_spin_array,
                         atoms::x_total_spin_field_array,
                         atoms::y_total_spin_field_array,
                         atoms::z_total_spin_field_array);
              
                         

             } //end spin for loop
        
        vmpi::mpi_complete_halo_swap();
        vmpi::barrier();


             
        for (int i=0; i<bdry_at;i++){
        atom = sim::b_octants[octant][i];
        const unsigned int imat = atoms::type_array[atom];
        lambda=mp::material[imat].alpha;
        spin_noise=mp::material[imat].H_th_sigma*sqrt(sim::temperature);
        
        
        calculate_spin_fields(atom, atom+1);
    	

        add_spin_noise(atom,
                    atom+1,
                    mp::dt_SI*1e12,
                    atoms::type_array, // type for atom
                    atoms::x_spin_array,
                    atoms::y_spin_array,
                    atoms::z_spin_array,
                    atoms::x_total_spin_field_array,
                    atoms::y_total_spin_field_array,
                    atoms::z_total_spin_field_array,
                    Hx_th, //  vectors for fields
                    Hy_th,
                    Hz_th);
                    



       cayley_update(atom,
                    atom+1,
                    cay_dt,
                    atoms::x_spin_array,
                    atoms::y_spin_array,
                    atoms::z_spin_array,
                    atoms::x_total_spin_field_array,
                    atoms::y_total_spin_field_array,
                    atoms::z_total_spin_field_array);
         

        } //end spin for loop

    vmpi::barrier();


          }	// end second octant loop onwards

          
  	
          
          //first octant loop return    
          for(int octant = 7; octant >= 0; octant--) {
          
                 vmpi::mpi_init_halo_swap();
                 int core_at=sim::c_octants[octant].size();
                 int bdry_at=sim::b_octants[octant].size();

             
                
                for (int i=core_at-1;i>=0;i--){
                atom = sim::c_octants[octant][i];
                const unsigned int imat = atoms::type_array[atom];
                lambda=mp::material[imat].alpha;
                spin_noise=mp::material[imat].H_th_sigma*sqrt(sim::temperature);
                
                calculate_spin_fields(atom, atom+1);
            	
                add_spin_noise(atom,
                            atom+1,
                            mp::dt_SI*1e12,
                            atoms::type_array, // type for atom
                            atoms::x_spin_array,
                            atoms::y_spin_array,
                            atoms::z_spin_array,
                            atoms::x_total_spin_field_array,
                            atoms::y_total_spin_field_array,
                            atoms::z_total_spin_field_array,
                            Hx_th, //  vectors for fields
                            Hy_th,
                            Hz_th);
                            



               cayley_update(atom,
                            atom+1,
                            cay_dt,
                            atoms::x_spin_array,
                            atoms::y_spin_array,
                            atoms::z_spin_array,
                            atoms::x_total_spin_field_array,
                            atoms::y_total_spin_field_array,
                            atoms::z_total_spin_field_array);
                 

                }//end spin loop
                
                vmpi::mpi_complete_halo_swap();
                vmpi::barrier();


                for (int i=bdry_at-1;i>=0;i--){
                atom = sim::b_octants[octant][i];
                
                const unsigned int imat = atoms::type_array[atom];
                lambda=mp::material[imat].alpha;
                spin_noise=mp::material[imat].H_th_sigma*sqrt(sim::temperature);
                
               calculate_spin_fields(atom, atom+1);


               add_spin_noise(atom,
                           atom+1,
                           mp::dt_SI*1e12,
                           atoms::type_array, // type for atom
                           atoms::x_spin_array,
                           atoms::y_spin_array,
                           atoms::z_spin_array,
                           atoms::x_total_spin_field_array,
                           atoms::y_total_spin_field_array,
                           atoms::z_total_spin_field_array,
                           Hx_th, //  vectors for fields
                           Hy_th,
                           Hz_th);
                           



              cayley_update(atom,
                           atom+1,
                           cay_dt,
                           atoms::x_spin_array,
                           atoms::y_spin_array,
                           atoms::z_spin_array,
                           atoms::x_total_spin_field_array,
                           atoms::y_total_spin_field_array,
                           atoms::z_total_spin_field_array);
                

                }//end spin loop

             vmpi::barrier();

             }      //end second octant loop return 
             
          

           
 // Swap timers compute -> wait
 vmpi::TotalComputeTime+=vmpi::SwapTimer(vmpi::ComputeTime, vmpi::WaitTime);

 // Wait for other processors
 vmpi::barrier();

 // Swap timers wait -> compute
 vmpi::TotalWaitTime += vmpi::SwapTimer(vmpi::WaitTime, vmpi::ComputeTime);
  
   
   return;

}

} // End of namespace sim
#endif
