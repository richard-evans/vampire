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
#include "vio.hpp"

// Field calculation functions
int calculate_spin_fields(const int,const int);
int calculate_external_fields(const int,const int);

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

   std::vector <double> mod_S;

   std::vector <double> tx;
   std::vector <double> ty;
   std::vector <double> tz;

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

   mod_S.resize(atoms::num_atoms,1.0);

   tx.resize(atoms::num_atoms,0.0);
   ty.resize(atoms::num_atoms,0.0);
   tz.resize(atoms::num_atoms,0.0);

   // Disable external thermal field calculations
   sim::hamiltonian_simulation_flags[3]=0;

   /*
   // Function to check if Landau coefficients are reasonably set
   std::vector<double> spin_length_init_container;
   std::vector<int> sample_atoms;
   bool complete_flag = false;
   int current_material = 0;
   for(int atom=0; complete_flag==false; atom++){ // Loop over all atoms until a sample atom is found for all materials
      if(atoms::type_array[atom]==current_material){
         sample_atoms.push_back(atom);
         current_material++;
      }
      if(sample_atoms.size()==mp::num_materials || atom==atoms::num_atoms-1) complete_flag = true;
   }
   spin_length_init_container.resize(500,0.0);
   for(int imaterial=0; imaterial<mp::num_materials; imaterial++){
      const double A = sim::internal::lsf_second_order_coefficient[imaterial];
      const double B = sim::internal::lsf_fourth_order_coefficient[imaterial];
      const double C = sim::internal::lsf_sixth_order_coefficient[imaterial];
      const double J = exchange::single_spin_energy(sample_atoms[imaterial],0.0,0.0,1.0); // Assume S=(0,0,1) when calculating initial exchange energy

      // Iterate spin length values between |S|= 0 - 5
      for(double sl=0.0; sl<500.0; sl++){
         const double mods = sl/100;
         const double landau_energy = A*mods*mods + B*mods*mods*mods*mods + C*mods*mods*mods*mods*mods*mods + J*mp::material[imaterial].mu_s_SI*mods;
         spin_length_init_container[sl] = landau_energy;
      }

      // Index of lowest energy value
      const int index = std::distance(std::begin(spin_length_init_container), std::min_element(std::begin(spin_length_init_container), std::end(spin_length_init_container)));
      if(index<=70 || index >=130){ // If |S| is between 0.7 - 1.3, accept coefficients. Otherwise, project error message
         terminaltextcolor(RED);
         std::cerr << "Error in LSF integration! - Landau coefficients set for material " << imaterial+1 << " initialise spin length too far from |S|=1!" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
   }
   */

   LSF_set = true;

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
   const double kB = 1.3806503e-23;

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

   // Thermal noise based on Gaussian function
   double sigma = ( sqrt( ( 2.0 * kB * sim::temperature * mp::gamma_SI) / ( mp::dt_SI ) ) );
   generate (tx.begin(),tx.begin()+num_atoms, mtrandom::gaussian);
   generate (ty.begin(),ty.begin()+num_atoms, mtrandom::gaussian);
   generate (tz.begin(),tz.begin()+num_atoms, mtrandom::gaussian);

   /*
   std::mt19937 temperature(seeding);
   for (int atom = 0; atom < num_atoms; atom++){

      std::normal_distribution<double> randomTemp(0, 1.0); 

      // Generate Gaussian
      tx[atom]=randomTemp(temperature);
      ty[atom]=randomTemp(temperature);
      tz[atom]=randomTemp(temperature);

   }
   */

   // Calculate first GSE step
   for (int atom = 0; atom < num_atoms; atom++){

      const int imaterial = atoms::type_array[atom];
      const double alpha = mp::material[atoms::type_array[atom]].alpha;
      const double mu = mp::material[atoms::type_array[atom]].mu_s_SI;

		// Store local spin in S and local field in H
		const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};

		const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom]+LSF_arrays::x_lsf_array[atom],
									atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom]+LSF_arrays::y_lsf_array[atom],
									atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]+LSF_arrays::z_lsf_array[atom]};

      // Calculate Delta S
      xyz[0]=(-mp::gamma_SI*(S[1]*H[2]-S[2]*H[1]))+(mp::gamma_SI*alpha*H[0])+(tx[atom]*(sigma * sqrt(alpha / mu)));
      xyz[1]=(-mp::gamma_SI*(S[2]*H[0]-S[0]*H[2]))+(mp::gamma_SI*alpha*H[1])+(ty[atom]*(sigma * sqrt(alpha / mu)));
      xyz[2]=(-mp::gamma_SI*(S[0]*H[1]-S[1]*H[0]))+(mp::gamma_SI*alpha*H[2])+(tz[atom]*(sigma * sqrt(alpha / mu)));

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
      const double alpha = mp::material[atoms::type_array[atom]].alpha;
      const double mu = mp::material[atoms::type_array[atom]].mu_s_SI;

		// Store local spin in S and local field in H
		const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};

		const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom]+LSF_arrays::x_lsf_array[atom],
									atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom]+LSF_arrays::y_lsf_array[atom],
									atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]+LSF_arrays::z_lsf_array[atom]};

      // Calculate Delta S
      xyz[0]=(-mp::gamma_SI*(S[1]*H[2]-S[2]*H[1]))+(mp::gamma_SI*alpha*H[0])+(tx[atom]*(sigma * sqrt(alpha / mu)));
      xyz[1]=(-mp::gamma_SI*(S[2]*H[0]-S[0]*H[2]))+(mp::gamma_SI*alpha*H[1])+(ty[atom]*(sigma * sqrt(alpha / mu)));
      xyz[2]=(-mp::gamma_SI*(S[0]*H[1]-S[1]*H[0]))+(mp::gamma_SI*alpha*H[2])+(tz[atom]*(sigma * sqrt(alpha / mu)));

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

   return;

}

} // End of indental namespace

} // End of sim namespace
