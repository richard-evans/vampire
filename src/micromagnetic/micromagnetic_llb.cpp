//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

#include "random.hpp"
#include "errors.hpp"
#include "atoms.hpp"
#include "cells.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include "dipole.hpp"
#include "vio.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>

namespace mm = micromagnetic::internal;
using namespace std;

//function declaration for the serial LLB
int LLB_serial_heun( std::vector <int> local_cell_array,
							int num_steps,
                     int num_cells,
						   int num_local_cells,
                     double temperature,
                     std::vector<double>& x_mag_array,
                     std::vector<double>& y_mag_array,
                     std::vector<double>& z_mag_array,
                     double Hx,
                     double Hy,
                     double Hz,
                     double H,
                     double dt,
                     std::vector <double> volume_array);

namespace micromagnetic_arrays{

	// Local arrays for LLG integration
	std::vector <double> x_euler_array;
	std::vector <double> y_euler_array;
	std::vector <double> z_euler_array;

	std::vector <double> x_array;
	std::vector <double> y_array;
	std::vector <double> z_array;

	std::vector <double> x_heun_array;
	std::vector <double> y_heun_array;
	std::vector <double> z_heun_array;

	std::vector <double> x_spin_storage_array;
	std::vector <double> y_spin_storage_array;
	std::vector <double> z_spin_storage_array;

	std::vector <double> x_initial_spin_array;
	std::vector <double> y_initial_spin_array;
	std::vector <double> z_initial_spin_array;

	std::vector <double> GW1x;
	std::vector <double> GW1y;
	std::vector <double> GW1z;
	std::vector <double> GW2x;
	std::vector <double> GW2y;
	std::vector <double> GW2z;

	bool LLG_set=false; ///< Flag to define state of LLG arrays (initialised/uninitialised)

}


namespace micromagnetic{

//---------------------------------------------------------------------------
// persistant file scope variables
//---------------------------------------------------------------------------
const double kB = 1.3806503e-23;
//arrays to store the magnetisation and the field
std::vector<double> m(3,0.0);
std::vector<double> spin_field(3,0.0);

//---------------------------------------------------------------------------
// Function to initialise micromagnetic integration arrays
//---------------------------------------------------------------------------
int micromagnetic_init(int num_cells,
                       std::vector<double>& x_mag_array,
                       std::vector<double>& y_mag_array,
                       std::vector<double>& z_mag_array){

   using namespace micromagnetic_arrays;

   x_spin_storage_array.resize(num_cells,0.0);
   y_spin_storage_array.resize(num_cells,0.0);
   z_spin_storage_array.resize(num_cells,0.0);

   x_array.resize(num_cells,0.0);
   y_array.resize(num_cells,0.0);
   z_array.resize(num_cells,0.0);

   x_initial_spin_array.resize(num_cells,0.0);
   y_initial_spin_array.resize(num_cells,0.0);
   z_initial_spin_array.resize(num_cells,0.0);

   x_euler_array.resize(num_cells,0.0);
   y_euler_array.resize(num_cells,0.0);
   z_euler_array.resize(num_cells,0.0);

   x_heun_array.resize(num_cells,0.0);
   y_heun_array.resize(num_cells,0.0);
   z_heun_array.resize(num_cells,0.0);

   GW1x.resize(num_cells,0.0);
   GW1y.resize(num_cells,0.0);
   GW1z.resize(num_cells,0.0);
   GW2x.resize(num_cells,0.0);
   GW2y.resize(num_cells,0.0);
   GW2z.resize(num_cells,0.0);

   LLG_set=true;

   return EXIT_SUCCESS;

}

//------------------------------------------------------------------------------
// Function to integrate LLB for micromagnetic cells
//------------------------------------------------------------------------------
int LLB( std::vector <int>& local_cell_array,
         int num_steps,
         int num_cells,
         int num_local_cells,
         double temperature,
         std::vector<double>& x_mag_array,
         std::vector<double>& y_mag_array,
         std::vector<double>& z_mag_array,
         double Hx,
         double Hy,
         double Hz,
         double H,
         double dt,
         std::vector <double>& volume_array){

   // import micromagnetic arrays for brevity
   using namespace micromagnetic_arrays;

//std::cout << num_cells <<std::endl;
	// Check for initialisation of LLG integration arrays
	if(LLG_set== false) micromagnetic::micromagnetic_init(num_cells, x_mag_array, y_mag_array, z_mag_array);

   // set cell temperatures
   for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
      int cell = list_of_micromagnetic_cells[lc];
      mm::T[cell] = temperature;
   }

   // Output cell parameters for debugging
   //mm::output_system_parameters(cells::pos_and_mom_array, x_mag_array, y_mag_array, z_mag_array); std::cin.get();

   // compute parallel and perpendicular susceptibilities
   mm::calculate_chi_para(number_of_micromagnetic_cells, list_of_micromagnetic_cells, mm::one_o_chi_para, mm::T, mm::Tc);
   mm::calculate_chi_perp(number_of_micromagnetic_cells, list_of_micromagnetic_cells, mm::one_o_chi_perp, mm::T, mm::Tc);

   //initialise the x_array to Mx/Ms
   for (int cell = 0; cell < num_cells; cell++){
		double ims =1.0/mm::ms[cell];
      x_array[cell] = x_mag_array[cell]*ims;
      y_array[cell] = y_mag_array[cell]*ims;
      z_array[cell] = z_mag_array[cell]*ims;
   }

   //save this new m as the initial value, so it can be saved and used in the final equation.
   for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
      int cell = list_of_micromagnetic_cells[lc];
      //converting units from J/T to M/Ms which has a range from 0-1.
      x_initial_spin_array[cell] = x_array[cell];
      y_initial_spin_array[cell] = y_array[cell];
      z_initial_spin_array[cell] = z_array[cell];
   }

   //The external fields equal the length of the field times the applied field vector.
   mm::ext_field[0] = H*Hx;
   mm::ext_field[1] = H*Hy;
   mm::ext_field[2] = H*Hz;

   //6 arrays of gaussian random numbers to store the stochastic noise terms for x,y,z parallel and perperdicular
   //fill the noise terms
   for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
      int cell = list_of_micromagnetic_cells[lc];
      GW1x[cell] = mtrandom::gaussian();
      GW1y[cell] = mtrandom::gaussian();
      GW1z[cell] = mtrandom::gaussian();
      GW2x[cell] = mtrandom::gaussian();
      GW2y[cell] = mtrandom::gaussian();
      GW2z[cell] = mtrandom::gaussian();
   }

   double xyz[3] = {0.0,0.0,0.0};

   //calculate the euler gradient (for cells on this processor)
   for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
      int cell = list_of_micromagnetic_cells[lc];

      //save to m for easy access
      m[0] = x_array[cell];
      m[1] = y_array[cell];
      m[2] = z_array[cell];

      //calculate spin fields
      spin_field = mm::calculate_llb_fields(m, temperature, num_cells, cell, x_array,y_array,z_array);

      //calculates the stochatic parallel and perpendicular terms
      double sigma_para = sqrt(2*kB*temperature*mm::alpha_para[cell]/(mm::ms[cell]*dt));
      double sigma_perp = sqrt(2*kB*temperature*(mm::alpha_perp[cell]-mm::alpha_para[cell])/(dt*mm::ms[cell]*mm::alpha_perp[cell]*mm::alpha_perp[cell]));
      const double H[3] = {spin_field[0], spin_field[1], spin_field[2]};


      //saves the noise terms to an array
      const double GW2t[3] = {GW2x[cell],GW2y[cell],GW2z[cell]};
      const double one_o_m_squared = 1.0/(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      const double SdotH = m[0]*H[0] + m[1]*H[1] + m[2]*H[2];

		//calculates delta terms
      xyz[0]=  - (m[1]*H[2]-m[2]*H[1])
               + mm::alpha_para[cell]*m[0]*SdotH*one_o_m_squared
               - mm::alpha_perp[cell]*(m[1]*(m[0]*H[1]-m[1]*H[0])-m[2]*(m[2]*H[0]-m[0]*H[2]))*one_o_m_squared
               + GW1x[cell]*sigma_para
               - mm::alpha_perp[cell]*(m[1]*(m[0]*GW2t[1]-m[1]*GW2t[0])-m[2]*(m[2]*GW2t[0]-m[0]*GW2t[2]))*one_o_m_squared*sigma_perp;

      xyz[1]=  - (m[2]*H[0]-m[0]*H[2])
               + mm::alpha_para[cell]*m[1]*SdotH*one_o_m_squared
               - mm::alpha_perp[cell]*(m[2]*(m[1]*H[2]-m[2]*H[1])-m[0]*(m[0]*H[1]-m[1]*H[0]))*one_o_m_squared
               + GW1y[cell]*sigma_para
               - mm::alpha_perp[cell]*(m[2]*(m[1]*GW2t[2]-m[2]*GW2t[1])-m[0]*(m[0]*GW2t[1]-m[1]*GW2t[0]))*one_o_m_squared*sigma_perp;

      xyz[2]=  - (m[0]*H[1]-m[1]*H[0])
               + mm::alpha_para[cell]*m[2]*SdotH*one_o_m_squared
               - mm::alpha_perp[cell]*(m[0]*(m[2]*H[0]-m[0]*H[2])-m[1]*(m[1]*H[2]-m[2]*H[1]))*one_o_m_squared
               + GW1z[cell]*sigma_para
               - mm::alpha_perp[cell]*(m[0]*(m[2]*GW2t[0]-m[0]*GW2t[2])-m[1]*(m[1]*GW2t[2]-m[2]*GW2t[1]))*one_o_m_squared*sigma_perp;

      x_euler_array[cell] = xyz[0];
      y_euler_array[cell] = xyz[1];
      z_euler_array[cell] = xyz[2];

   }


	// For parallel version set arrays to zero to allow parallel reduction to work
   #ifdef MPICF
      const int data_size = x_mag_array.size();
      for(int cell=0; cell < data_size; cell++){
         x_spin_storage_array[cell] = 0.0;
         y_spin_storage_array[cell] = 0.0;
         z_spin_storage_array[cell] = 0.0;
      }
   #endif

   //these new x postiion are stored in an array (store)
   //x = x+step*dt
   for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
      int cell = list_of_micromagnetic_cells[lc];
      x_spin_storage_array[cell] = x_array[cell] + x_euler_array[cell]*dt;
      y_spin_storage_array[cell] = y_array[cell] + y_euler_array[cell]*dt;
      z_spin_storage_array[cell] = z_array[cell] + z_euler_array[cell]*dt;
   }

	// Reduce cell magnetizations on all processors to enable correct exchange field calculations
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &x_spin_storage_array[0], data_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &y_spin_storage_array[0], data_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &z_spin_storage_array[0], data_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

   //calculates the heun gradient
   for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
   	int cell = list_of_micromagnetic_cells[lc];
   	m[0] = x_spin_storage_array[cell];
   	m[1] = y_spin_storage_array[cell];
   	m[2] = z_spin_storage_array[cell];

   	spin_field = mm::calculate_llb_fields(m, temperature, num_cells, cell, x_spin_storage_array,y_spin_storage_array,z_spin_storage_array);

   	//fill the noise terms
   	//calculates the stochatic parallel and perpendicular terms
   	double sigma_para = sqrt(2*kB*temperature*mm::alpha_para[cell]/(mm::ms[cell]*dt)); //why 1e-27
   	double sigma_perp = sqrt(2*kB*temperature*(mm::alpha_perp[cell]-mm::alpha_para[cell])/(dt*mm::ms[cell]*mm::alpha_perp[cell]*mm::alpha_perp[cell]));

   	const double H[3] = {spin_field[0], spin_field[1], spin_field[2]};

   	//saves the noise terms to an array
   	const double GW2t[3] = {GW2x[cell],GW2y[cell],GW2z[cell]};
   	const double one_o_m_squared = 1.0/(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
   	const double SdotH = m[0]*H[0] + m[1]*H[1] + m[2]*H[2];

   	double xyz[3] = {0.0,0.0,0.0};
   	//calculates the LLB deltas
   	xyz[0]=  - (m[1]*H[2]-m[2]*H[1])
               + mm::alpha_para[cell]*m[0]*SdotH*one_o_m_squared
               - mm::alpha_perp[cell]*(m[1]*(m[0]*H[1]-m[1]*H[0])-m[2]*(m[2]*H[0]-m[0]*H[2]))*one_o_m_squared
               + GW1x[cell]*sigma_para
               - mm::alpha_perp[cell]*(m[1]*(m[0]*GW2t[1]-m[1]*GW2t[0])-m[2]*(m[2]*GW2t[0]-m[0]*GW2t[2]))*one_o_m_squared*sigma_perp;

   	xyz[1]=  - (m[2]*H[0]-m[0]*H[2])
               + mm::alpha_para[cell]*m[1]*SdotH*one_o_m_squared
               - mm::alpha_perp[cell]*(m[2]*(m[1]*H[2]-m[2]*H[1])-m[0]*(m[0]*H[1]-m[1]*H[0]))*one_o_m_squared
               + GW1y[cell]*sigma_para
               - mm::alpha_perp[cell]*(m[2]*(m[1]*GW2t[2]-m[2]*GW2t[1])-m[0]*(m[0]*GW2t[1]-m[1]*GW2t[0]))*one_o_m_squared*sigma_perp;

   	xyz[2]=	- (m[0]*H[1]-m[1]*H[0])
               + mm::alpha_para[cell]*m[2]*SdotH*one_o_m_squared
               - mm::alpha_perp[cell]*(m[0]*(m[2]*H[0]-m[0]*H[2])-m[1]*(m[1]*H[2]-m[2]*H[1]))*one_o_m_squared
               + GW1z[cell]*sigma_para
               - mm::alpha_perp[cell]*(m[0]*(m[2]*GW2t[0]-m[0]*GW2t[2])-m[1]*(m[1]*GW2t[2]-m[2]*GW2t[1]))*one_o_m_squared*sigma_perp;

      x_heun_array[cell] = xyz[0];
   	y_heun_array[cell] = xyz[1];
   	z_heun_array[cell] = xyz[2];

   }

	// For parallel version set arrays to zero to allow parallel reduction to work
	#ifdef MPICF
		for(int cell=0; cell< data_size; cell++){
			x_array[cell] = 0.0;
			y_array[cell] = 0.0;
			z_array[cell] = 0.0;
			cells::mag_array_x[cell] = 0.0;
			cells::mag_array_y[cell] = 0.0;
			cells::mag_array_z[cell] = 0.0;
		}
	#endif

   //update spin arrays
   for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){
   	int cell = list_of_micromagnetic_cells[lc];

   	 //m = initial + 1/2 dt (euler + heun)
   	 x_array[cell] = x_initial_spin_array[cell] + 0.5*dt*(x_euler_array[cell] + x_heun_array[cell]);
   	 y_array[cell] = y_initial_spin_array[cell] + 0.5*dt*(y_euler_array[cell] + y_heun_array[cell]);
   	 z_array[cell] = z_initial_spin_array[cell] + 0.5*dt*(z_euler_array[cell] + z_heun_array[cell]);

   	//convert from M/Ms to J/T
   	cells::mag_array_x[cell] = x_array[cell]*mm::ms[cell];
   	cells::mag_array_y[cell] = y_array[cell]*mm::ms[cell];
   	cells::mag_array_z[cell] = z_array[cell]*mm::ms[cell];
   }

	// Reduce unit vectors and moments to all processors
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_x[0], data_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_y[0], data_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_z[0], data_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &x_array[0],            data_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &y_array[0],            data_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &z_array[0],            data_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   #endif

   //update atom positions
   if (discretisation_type == 2 || sim::time%vout::output_rate -1){
      for(int atom_list=0;atom_list<number_of_none_atomistic_atoms;atom_list++){
         int atom = list_of_none_atomistic_atoms[atom_list];
         int cell = cells::atom_cell_id_array[atom];
         atoms::x_spin_array[atom] = x_array[cell];
         atoms::y_spin_array[atom] = y_array[cell];
         atoms::z_spin_array[atom] = z_array[cell];
      }
    }

   if (sim::time%100000 && enable_resistance && mm::resistance_layer_2 != mm::resistance_layer_1 ){
      micromagnetic::MR_resistance = mm::calculate_resistance();
   }

	return 0;

}


} // end of micromagnetic namespace
