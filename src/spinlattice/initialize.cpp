//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
// Vampire headers
#include "sld.hpp"
#include "atoms.hpp"
#include "iostream"
#include "neighbours.hpp"
#include "material.hpp"
#include "constants.hpp"



// sld module headers
#include "internal.hpp"


namespace sld{


   //----------------------------------------------------------------------------
   // Function to initialize sld module
   //----------------------------------------------------------------------------
   void initialize(){
      std::cout<<"Input parameters for Spin-lattice dynamics simulations:"<<std::endl;
      std::cout<<"*******************************************************"<<std::endl;

      std::cout<<"Potential Cutoff: "<<sld::internal::r_cut_pot<<std::endl;
      std::cout<<"Fields Cutoff: "<<sld::internal::r_cut_fields<<std::endl;
      std::cout<<"Mass: "<<sld::internal::mp[0].mass.get()<<std::endl;
      std::cout<<"Lattice damping: "<<sld::internal::mp[0].damp_lat.get()<<std::endl;
      std::cout<<"Coupling C0 "<<sld::internal::mp[0].C0.get()<<std::endl;
      if(sld::internal::harmonic)std::cout<<"Harmonic potential is used of potential well depth V0="<<sld::internal::mp[0].V0.get()<<std::endl;
      if(sld::internal::pseudodipolar)std::cout<<"Pseudodipolar coupling is used of strength C0="<<sld::internal::mp[0].C0.get()<<std::endl;
      std::cout<<"*******************************************************"<<std::endl;
    


     sld::internal::initialise_positions(sld::internal::x0_coord_array, // coord vectors for atoms
                 sld::internal::y0_coord_array,
                 sld::internal::z0_coord_array,
                 atoms::x_coord_array, // coord vectors for atoms
                 atoms::y_coord_array,
                 atoms::z_coord_array,
                 sld::internal::dr_init);



    //initialise exchange, coupling parameters
     sld::internal::initialise_sld_parameters();
     

     
     sld::internal::thermal_velocity(atoms::x_velo_array, atoms::y_velo_array,atoms::z_velo_array);
     
     
     //initialise for Parallel simulations
     //Initialize parallel mc variables
     suzuki_trotter_parallel_initialized = false;
     internal::c_octants.resize(8);
     internal::b_octants.resize(8);
     
     internal::all_atoms_octant_start_index.reserve(8); 
     internal::all_atoms_octant_end_index.reserve(8); 
     internal::all_atoms_octant.reserve(atoms::num_atoms);
     
    
     
    // sld::tests();
    
    

      return;


   }


   namespace internal{

   void initialise_positions(std::vector<double>& x0_coord_array, // coord vectors for atoms
               std::vector<double>& y0_coord_array,
               std::vector<double>& z0_coord_array,
               std::vector<double>& x_coord_array, // coord vectors for atoms
               std::vector<double>& y_coord_array,
               std::vector<double>& z_coord_array,
               const double dr){


   x0_coord_array.resize(atoms::num_atoms,0);
   y0_coord_array.resize(atoms::num_atoms,0);
   z0_coord_array.resize(atoms::num_atoms,0);
   
 
 
   x_coord_storage_array.resize(atoms::num_atoms,0);
   y_coord_storage_array.resize(atoms::num_atoms,0);
   z_coord_storage_array.resize(atoms::num_atoms,0);

   for( int i = 0; i < atoms::num_atoms; i++){

      x0_coord_array[i]=x_coord_array[i];
      y0_coord_array[i]=y_coord_array[i];
      z0_coord_array[i]=z_coord_array[i];
      
     x_coord_storage_array[i]=x_coord_array[i];
     y_coord_storage_array[i]=y_coord_array[i];
     z_coord_storage_array[i]=z_coord_array[i];
            
      

   }
   //std::cout<<"random positions  "<<dr<<std::endl;
   /*
    for( int i = 0; i < atoms::num_atoms; i++)
    {

             x_coord_array[i] += dr* (2.0*rand()/double(RAND_MAX) -1.0);
             y_coord_array[i] += dr* (2.0*rand()/double(RAND_MAX) -1.0);
             z_coord_array[i] += dr* (2.0*rand()/double(RAND_MAX) -1.0); 

        
    }*/
    
           std::cout<<"pos after"<<x_coord_array[0]<<"\t"<<z_coord_array[10]<<std::endl;

   return;
   
   
   }//end of initialise initialise_positions
   
   
   void thermal_velocity(std::vector<double>& x_velo_array, // velocities vectors
                         std::vector<double>& y_velo_array,
                         std::vector<double>& z_velo_array)
                         
   {    sld::internal::th_velo=10;
       const double v_therm = sqrt( 3.0 * constants::kB_eV * sld::internal::th_velo / sld::internal::mp[0].mass.get());
       int N=atoms::num_atoms;
       double rx,ry,s;
       double net_v[] = {0.0, 0.0, 0.0};
       
       std::cout<<"velocities before "<<x_velo_array[0]<<"\t"<<x_velo_array[10]<<"\t"<<th_velo<<std::endl;

       for( int i = 0; i < N; i++)
       {
           /*s = 2.0;
           while( s >= 1.0)
           {
               rx = 2.0;//2.0* (rand()/double(RAND_MAX)) - 1.0;
               ry = 2.0;//2.0* (rand()/double(RAND_MAX)) - 1.0;

               s = rx*rx + ry*ry;
               
           }*/
           
            s=0.5;
            rx=2.0;
            ry=2.0;

           x_velo_array[i] = v_therm * 2.0 * rx * sqrt( 1.0 - s);
           y_velo_array[i] = v_therm * 2.0 * ry * sqrt( 1.0 - s);
           z_velo_array[i] = v_therm * ( 1.0 - 2.0 * s);

          // net_v[0] += x_velo_array[i];
          // net_v[1] += y_velo_array[i];
          // net_v[2] += z_velo_array[i];
       }

       //net_v[0] /= double(N);
       //net_v[1] /= double(N);
       //net_v[2] /= double(N);

      /* for( int i = 0; i < N; i++)
       {
           x_velo_array[i] -= net_v[0];
           y_velo_array[i] -= net_v[1];
           z_velo_array[i] -= net_v[2];
       }
       */
       std::cout<<"velocities after"<<x_velo_array[0]<<"\t"<<y_velo_array[0]<<"\t"<<std::endl;
       return;

   }
   
   



   void initialise_sld_parameters(){

      sld::internal::forces_array_x.resize(atoms::num_atoms,0);
      sld::internal::forces_array_y.resize(atoms::num_atoms,0);
      sld::internal::forces_array_z.resize(atoms::num_atoms,0);

      sld::internal::fields_array_x.resize(atoms::num_atoms,0);
      sld::internal::fields_array_y.resize(atoms::num_atoms,0);
      sld::internal::fields_array_z.resize(atoms::num_atoms,0);

      sld::internal::velo_array_x.resize(atoms::num_atoms,0);
      sld::internal::velo_array_y.resize(atoms::num_atoms,0);
      sld::internal::velo_array_z.resize(atoms::num_atoms,0);

      sld::internal::potential_eng.resize(atoms::num_atoms,0);
      sld::internal::sumJ.resize(atoms::num_atoms,0);
      sld::internal::sumC.resize(atoms::num_atoms,0);
      sld::internal::exch_eng.resize(atoms::num_atoms,0);
      sld::internal::coupl_eng.resize(atoms::num_atoms,0);





      //sqrt( 2.0 * eta * consts::kB * T / (mass *dt));
      for(int mat=0;mat<mp::num_materials; mat++){

         double C=sld::internal::mp[mat].C0.get();
         sld::internal::mp[mat].C0.set(sld::internal::mp[mat].J0.get()*C);
         sld::internal::mp[mat].J0_ms.set(sld::internal::mp[mat].J0.get()/mp::material[mat].mu_s_SI);
         sld::internal::mp[mat].C0_ms.set(sld::internal::mp[mat].C0.get()/mp::material[mat].mu_s_SI);
         sld::internal::mp[mat].J0_prime.set(3.0*sld::internal::mp[mat].J0.get()/sld::internal::r_cut_fields);
         sld::internal::mp[mat].F_th_sigma.set(sqrt(2.0*sld::internal::mp[0].damp_lat.get()*constants::kB_eV / (sld::internal::mp[mat].mass.get()*mp::dt_SI*1e12)));

         //now change to ev the following
         sld::internal::mp[0].V0.set(sld::internal::mp[0].V0.get()*6.242e18);
         
         //initialise mass for the rest of the simulations
         mp::material[mat].mass=sld::internal::mp[0].mass.get();
       }


       //initialise mass array
       
       
       for(int atom=0;atom<atoms::num_atoms;atom++){

   		int mat=atoms::type_array[atom];
        atoms::mass_spin_array[atom]=mp::material[mat].mass;
        }
   	
       //test cayley_update
      /* sld::internal::cayley_update(0,
                   0,
                   -mp::dt/4.0,
                   atoms::x_spin_array,
                   atoms::y_spin_array,
                   atoms::z_spin_array,
                   sld::internal::fields_array_x,
                   sld::internal::fields_array_y,
                   sld::internal::fields_array_z);
*/


   return;
   }

} //end of internal
} // end of sld namespace
