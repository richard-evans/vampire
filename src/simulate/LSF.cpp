// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <random>
#include <fstream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "exchange.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "constants.hpp"
#include "random.hpp"
#include "LSF.hpp"
#include "internal.hpp"

// Field calculation functions
int calculate_spin_fields(const int,const int);
int calculate_external_fields(const int,const int);

// TESTING***
std::vector <int> mod_S_round;
int counter;
int loc_counter;
double mod_S_pt;
double mod_S_fe;
double exchsum;
double ptdir[3];
double k_B;

namespace LSF_arrays{

   // Local arrays for LSF integration
   std::vector <double> x_lsf_array;
   std::vector <double> y_lsf_array;
   std::vector <double> z_lsf_array;

   std::vector <double> x_euler_array;
   std::vector <double> y_euler_array;
   std::vector <double> z_euler_array;

   std::vector <double> x_heun_array;
   std::vector <double> y_heun_array;
   std::vector <double> z_heun_array;

   std::vector <double> x_spin_storage_array;
   std::vector <double> y_spin_storage_array;
   std::vector <double> z_spin_storage_array;

   std::vector <double> x_initial_spin_array;
   std::vector <double> y_initial_spin_array;
   std::vector <double> z_initial_spin_array;

   // Flag to define state of LSF arrays (initialised/uninitialised)
   bool LSF_set=false;

   // *** TESTING ***
   std::ofstream MyFile("spinlength");
   double simtemp=1;
   std::vector <double> mod_S;

}
namespace sim{

int LSFinit(){
   // Check calling of routine if error checking is activated
   if(err::check==true){std::cout << "sim:LSF_init has been called" << std::endl;}

   using namespace LSF_arrays;

   x_lsf_array.resize(atoms::num_atoms,0.0);
   y_lsf_array.resize(atoms::num_atoms,0.0);
   z_lsf_array.resize(atoms::num_atoms,0.0);

   x_spin_storage_array.resize(atoms::num_atoms,0.0);
	y_spin_storage_array.resize(atoms::num_atoms,0.0);
	z_spin_storage_array.resize(atoms::num_atoms,0.0);

	x_initial_spin_array.resize(atoms::num_atoms,0.0);
	y_initial_spin_array.resize(atoms::num_atoms,0.0);
	z_initial_spin_array.resize(atoms::num_atoms,0.0);

	x_euler_array.resize(atoms::num_atoms,0.0);
	y_euler_array.resize(atoms::num_atoms,0.0);
	z_euler_array.resize(atoms::num_atoms,0.0);

	x_heun_array.resize(atoms::num_atoms,0.0);
	y_heun_array.resize(atoms::num_atoms,0.0);
	z_heun_array.resize(atoms::num_atoms,0.0);

	LSF_set=true;

   sim::hamiltonian_simulation_flags[3]=0;

   mod_S.resize(atoms::num_atoms, 1.0);
   mod_S_round.resize(401, 0);
   counter=0;
   loc_counter=0;
   mod_S_pt=0;
   mod_S_fe=0;
   exchsum = 0.0;
   k_B=1.38064e-23;
   for(int i=0; i<3; i++){
      ptdir[i]=0;
   }

  	return EXIT_SUCCESS;

}

// LSF magnetic field function
void calculate_lsf_magnetic_field(const int start_index, const int end_index){

   // Check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_lsf_magnetic_field has been called" << std::endl;}
      
   // LSF Hamiltonian calculation
   for(int atom=start_index;atom<end_index;atom++){
         
      const int imaterial = atoms::type_array[atom];

      const double sx = atoms::x_spin_array[atom];
      const double sy = atoms::y_spin_array[atom];
      const double sz = atoms::z_spin_array[atom];

      const double imu_S = -1.0/mp::material[imaterial].mu_s_SI;

      const double L2 = 2.0*sim::internal::lsf_second_order_coefficient[imaterial]*imu_S;
      const double L4 = 4.0*sim::internal::lsf_fourth_order_coefficient[imaterial]*imu_S;
      const double L6 = 6.0*sim::internal::lsf_sixth_order_coefficient[imaterial]*imu_S;

      const double ss2 = sx*sx + sy*sy + sz*sz;   

      LSF_arrays::x_lsf_array[atom] = L2*sx + L4*sx*ss2 + L6*sx*ss2*ss2;
      LSF_arrays::y_lsf_array[atom] = L2*sy + L4*sy*ss2 + L6*sy*ss2*ss2;
      LSF_arrays::z_lsf_array[atom] = L2*sz + L4*sz*ss2 + L6*sz*ss2*ss2;
      
   } 

}

namespace internal{

void lsf_step(){

   // Check calling of routine if error checking is activated
   if(err::check==true){std::cout << "sim::LSF has been called" << std::endl;}

   using namespace LSF_arrays;

   // Check for initialisation of LSF integration arrays
   if(LSF_set==false) sim::LSFinit();

   // Local variables for system integration
   const int num_atoms=atoms::num_atoms;
   double xyz[3]; // Local delta spin components
   double S_new[3]; // New local spin moment

   // Calculate fields
   calculate_spin_fields(0, num_atoms);
   calculate_lsf_magnetic_field(0, num_atoms);
   calculate_external_fields(0, num_atoms); 

   // Store initial spin positions
   for (int atom = 0; atom < num_atoms; atom++){
      x_initial_spin_array[atom]=atoms::x_spin_array[atom];
		y_initial_spin_array[atom]=atoms::y_spin_array[atom];
		z_initial_spin_array[atom]=atoms::z_spin_array[atom];
   }

   // *** TESTING ***
   double spinlength[2];
   for (int atom = 0; atom < num_atoms; atom++){

      const int imaterial = atoms::type_array[atom];

      const double sx = atoms::x_spin_array[atom];
      const double sy = atoms::y_spin_array[atom];
      const double sz = atoms::z_spin_array[atom];

      spinlength[imaterial] = sqrt(sx*sx + sy*sy + sz*sz);

      if(imaterial==1){
         ptdir[0]+=sx/spinlength[imaterial];
         ptdir[1]+=sy/spinlength[imaterial];
         ptdir[2]+=sz/spinlength[imaterial];
      }

      if(imaterial==0){
         mod_S_fe += spinlength[imaterial];
      }

      if(imaterial==1){
         mod_S_pt += spinlength[imaterial];
      }

      int spinround = floor((spinlength[imaterial]*100)+0.5);
      mod_S_round[spinround]+= 1;
      
      //double exch = exchange::single_spin_energy(atom, atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]);
      //exchsum += exch;

   }
   if(simtemp!=sim::temperature){
      if(loc_counter==0) loc_counter=1;
      MyFile << sim::temperature << " " << mod_S_fe/((atoms::num_atoms/mp::num_materials)*loc_counter) << " " << mod_S_pt/((atoms::num_atoms/mp::num_materials)*loc_counter) << std::endl;
      //MyFile << sim::temperature << " " << exchsum/(loc_counter*atoms::num_atoms) << std::endl;
      loc_counter=0;
      exchsum=0;
      mod_S_pt=0;
      mod_S_fe=0;
      simtemp=sim::temperature;
   }
   spinlength[0]=0;
   spinlength[1]=0;

   std::vector <double> tx(num_atoms), ty(num_atoms), tz(num_atoms);

   //generate (tx.begin(),tx.begin()+num_atoms, mtrandom::gaussian);
   //generate (ty.begin(),ty.begin()+num_atoms, mtrandom::gaussian);
   //generate (tz.begin(),tz.begin()+num_atoms, mtrandom::gaussian);

   double sigma = ( sqrt( ( 2.0 * k_B * sim::temperature * mp::gamma_SI) / ( mp::dt_SI ) ) );

   std::mt19937 temperature(counter);

   for (int atom = 0; atom < num_atoms; atom++){

      //std::normal_distribution<double> randomTemp(0, (sigma * sqrt(mp::material[atoms::type_array[atom]].alpha / mp::material[atoms::type_array[atom]].mu_s_SI)) );
      std::normal_distribution<double> randomTemp(0, 1.0); 

      tx[atom]=randomTemp(temperature);
      ty[atom]=randomTemp(temperature);
      tz[atom]=randomTemp(temperature);

   }

   // Calculate first GSE step
   for (int atom = 0; atom < num_atoms; atom++){

      const int imaterial = atoms::type_array[atom];

		// Store local spin in S and local field in H
		const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};

		const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom]+LSF_arrays::x_lsf_array[atom],
									atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom]+LSF_arrays::y_lsf_array[atom],
									atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]+LSF_arrays::z_lsf_array[atom]};

      // Calculate Delta S
      xyz[0]=(-mp::gamma_SI*(S[1]*H[2]-S[2]*H[1]))+(mp::gamma_SI*mp::material[imaterial].alpha*H[0])+(tx[atom]*(sigma * sqrt(mp::material[atoms::type_array[atom]].alpha / mp::material[atoms::type_array[atom]].mu_s_SI)));
      xyz[1]=(-mp::gamma_SI*(S[2]*H[0]-S[0]*H[2]))+(mp::gamma_SI*mp::material[imaterial].alpha*H[1])+(ty[atom]*(sigma * sqrt(mp::material[atoms::type_array[atom]].alpha / mp::material[atoms::type_array[atom]].mu_s_SI)));
      xyz[2]=(-mp::gamma_SI*(S[0]*H[1]-S[1]*H[0]))+(mp::gamma_SI*mp::material[imaterial].alpha*H[2])+(tz[atom]*(sigma * sqrt(mp::material[atoms::type_array[atom]].alpha / mp::material[atoms::type_array[atom]].mu_s_SI)));

      // Store dS in euler array
		x_euler_array[atom]=xyz[0];
		y_euler_array[atom]=xyz[1];
		z_euler_array[atom]=xyz[2];

      // Calculate Euler Step
		S_new[0]=S[0]+xyz[0]*mp::dt_SI;
		S_new[1]=S[1]+xyz[1]*mp::dt_SI;
		S_new[2]=S[2]+xyz[2]*mp::dt_SI;

      //Writing of Spin Values to Storage Array
		x_spin_storage_array[atom]=S_new[0];
		y_spin_storage_array[atom]=S_new[1];
		z_spin_storage_array[atom]=S_new[2];

   }
   // Copy new spins to spin array
	for(int atom=0;atom<num_atoms;atom++){
		atoms::x_spin_array[atom]=x_spin_storage_array[atom];
		atoms::y_spin_array[atom]=y_spin_storage_array[atom];
		atoms::z_spin_array[atom]=z_spin_storage_array[atom];

	}

   // Recalculate spin dependent fields
   calculate_spin_fields(0, num_atoms);
   calculate_lsf_magnetic_field(0, num_atoms);

   // Calculate second GSE step
	for(int atom=0;atom<num_atoms;atom++){

		const int imaterial=atoms::type_array[atom];

		// Store local spin in S and local field in H
		const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};

		const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom]+LSF_arrays::x_lsf_array[atom],
									atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom]+LSF_arrays::y_lsf_array[atom],
									atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]+LSF_arrays::z_lsf_array[atom]};

      // Calculate Delta S
      xyz[0]=(-mp::gamma_SI*(S[1]*H[2]-S[2]*H[1]))+(mp::gamma_SI*mp::material[imaterial].alpha*H[0])+(tx[atom]*(sigma * sqrt(mp::material[atoms::type_array[atom]].alpha / mp::material[atoms::type_array[atom]].mu_s_SI)));
      xyz[1]=(-mp::gamma_SI*(S[2]*H[0]-S[0]*H[2]))+(mp::gamma_SI*mp::material[imaterial].alpha*H[1])+(ty[atom]*(sigma * sqrt(mp::material[atoms::type_array[atom]].alpha / mp::material[atoms::type_array[atom]].mu_s_SI)));
      xyz[2]=(-mp::gamma_SI*(S[0]*H[1]-S[1]*H[0]))+(mp::gamma_SI*mp::material[imaterial].alpha*H[2])+(tz[atom]*(sigma * sqrt(mp::material[atoms::type_array[atom]].alpha / mp::material[atoms::type_array[atom]].mu_s_SI)));

      // Store dS in Heun array
		x_heun_array[atom]=xyz[0];
		y_heun_array[atom]=xyz[1];
		z_heun_array[atom]=xyz[2];

   }
   // Calculate Heun Step
	for(int atom=0;atom<num_atoms;atom++){
		S_new[0]=x_initial_spin_array[atom]+(0.5*mp::dt_SI*(x_euler_array[atom]+x_heun_array[atom]));
		S_new[1]=y_initial_spin_array[atom]+(0.5*mp::dt_SI*(y_euler_array[atom]+y_heun_array[atom]));
		S_new[2]=z_initial_spin_array[atom]+(0.5*mp::dt_SI*(z_euler_array[atom]+z_heun_array[atom]));

      // Copy new spins to spin array
		atoms::x_spin_array[atom]=S_new[0];
		atoms::y_spin_array[atom]=S_new[1];
		atoms::z_spin_array[atom]=S_new[2];
   }

   // Store spin length data
   for(int atom=0;atom<num_atoms;atom++){
      const double sx = atoms::x_spin_array[atom];
      const double sy = atoms::y_spin_array[atom];
      const double sz = atoms::z_spin_array[atom];

      mod_S[atom] = sqrt(sx*sx + sy*sy + sz*sz);
   }

   counter++;
   loc_counter++;

   
   // Write spin length sampling data
   if(counter==sim::total_time){
      for (int i=0; i<401; i++){
         //MyFile << i << " " << mod_S_round[i] << std::endl;  
      }
      //std::cout << ptdir[0]/(counter*2662) << " " << ptdir[1]/(counter*2662) << " " << ptdir[2]/(counter*2662) << std::endl;
   }
   

   return;

}

} // End of indental namespace

} // End of sim namespace
