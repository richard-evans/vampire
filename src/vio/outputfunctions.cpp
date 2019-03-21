//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and Rory Pond 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk and rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
// Headers
#include "vio.hpp"
#include "grains.hpp"
#include "material.hpp"
#include "montecarlo.hpp"
#include "stats.hpp"
#include "sim.hpp"
#include "spintransport.hpp"

// vio module headers
#include "internal.hpp"

namespace vout{
	// Output Function 0
	void time(std::ostream& stream){
		stream << sim::time << "\t";
	}

	// Output Function 1
	void real_time(std::ostream& stream){
		stream << sim::time*mp::dt_SI << "\t";
	}

	// Output Function 2
	void temperature(std::ostream& stream){
		stream << sim::temperature << "\t";
	}

	// Output Function 3
	void Happ(std::ostream& stream){
		stream << sim::H_applied << "\t";
	}

	// Output Function 4
	void Hvec(std::ostream& stream){
		stream << sim::H_vec[0] << "\t"<< sim::H_vec[1] << "\t"<< sim::H_vec[2] << "\t";
	}

	// Output Function 5
	void mvec(std::ostream& stream){
		stream << stats::system_magnetization.output_normalized_magnetization();
	}

	// Output Function 6
	void magm(std::ostream& stream){
		stream << stats::system_magnetization.output_normalized_magnetization_length() << "\t";
	}

	// Output Function 7
	void mean_magm(std::ostream& stream){
		stream << stats::system_magnetization.output_normalized_mean_magnetization_length();
	}

	// Output Function 8
	void mat_mvec(std::ostream& stream){
		stream << stats::material_magnetization.output_normalized_magnetization();
	}

	// Output Function 9
	void mat_mean_magm(std::ostream& stream){
		stream << stats::material_magnetization.output_normalized_mean_magnetization_length();
	}

	// Output Function 10
	void grain_mvec(std::ostream& stream){

		unsigned int id=0; // grain id (excluding grains with zero atoms)

		// loop over all grains
		for(int grain=0;grain<grains::num_grains;grain++){
			// check for grains with zero atoms
			if(grains::grain_size_array[grain]!=0){
				stream << grains::x_mag_array[grain] << "\t";
				stream << grains::y_mag_array[grain] << "\t";
				stream << grains::z_mag_array[grain] << "\t";
				stream << grains::mag_m_array[grain] << "\t";
				id++;
			}
		}
	}

	// Output Function 11
	void grain_magm(std::ostream& stream){

		unsigned int id=0; // grain id (excluding grains with zero atoms)

		// loop over all grains
		for(int grain=0;grain<grains::num_grains;grain++){
			// check for grains with zero atoms
			if(grains::grain_size_array[grain]!=0){
				stream << grains::mag_m_array[grain] << "\t";
				id++;
			}
		}
	}

	// Output Function 12
	void mdoth(std::ostream& stream){
		// initialise vector of H
		std::vector<double> H(&sim::H_vec[0], &sim::H_vec[0]+3);
		stream << stats::system_magnetization.output_normalized_magnetization_dot_product(H);
	}

	// Output Function 13
	void grain_mat_mvec(std::ostream& stream){

		grains::output_mat_mag(stream);

	}

	// Output Function 14
	void systorque(std::ostream& stream){
		stream << stats::total_system_torque[0] << "\t";
		stream << stats::total_system_torque[1] << "\t";
		stream << stats::total_system_torque[2] << "\t";
	}

	// Output Function 15
	void mean_systorque(std::ostream& stream){
		stream << stats::total_mean_system_torque[0]/stats::torque_data_counter << "\t";
		stream << stats::total_mean_system_torque[1]/stats::torque_data_counter << "\t";
		stream << stats::total_mean_system_torque[2]/stats::torque_data_counter << "\t";
	}

	// Output Function 16
	void constraint_phi(std::ostream& stream){
		stream << sim::constraint_phi << "\t";
	}

	// Output Function 17
	void constraint_theta(std::ostream& stream){
		stream << sim::constraint_theta << "\t";
	}

	// Output Function 18
	void material_constraint_phi(std::ostream& stream){
		for(int mat=0;mat<mp::num_materials;mat++){
			stream << montecarlo::cmc::cmc_mat[mat].constraint_phi << "\t";
		}
	}

	// Output Function 19
	void material_constraint_theta(std::ostream& stream){
		for(int mat=0;mat<mp::num_materials;mat++){
			stream << montecarlo::cmc::cmc_mat[mat].constraint_theta << "\t";
		}
	}

	// Output Function 20
	void material_mean_systorque(std::ostream& stream){
		for(int mat=0;mat<mp::num_materials;mat++){
			stream << stats::sublattice_mean_torque_x_array[mat]/stats::torque_data_counter << "\t";
			stream << stats::sublattice_mean_torque_y_array[mat]/stats::torque_data_counter << "\t";
			stream << stats::sublattice_mean_torque_z_array[mat]/stats::torque_data_counter << "\t";
		}
	}

	// Output Function 21
	void mean_system_susceptibility(std::ostream& stream){
		stream << stats::system_susceptibility.output_mean_susceptibility(sim::temperature);
	}

	// Output Function 22
	void phonon_temperature(std::ostream& stream){
		stream << sim::TTTp << "\t";
	}

	// Output Function 23
	void material_temperature(std::ostream& stream){
		for(int mat=0;mat<mp::material.size();mat++){
			stream << mp::material[mat].temperature << "\t";
		}
	}

	// Output Function 24
	void material_applied_field_strength(std::ostream& stream){
		for(int mat=0;mat<mp::material.size();mat++){
			stream << mp::material[mat].applied_field_strength << "\t";
		}
	}

	// Output Function 25
	void material_fmr_field_strength(std::ostream& stream){
		const double real_time=sim::time*mp::dt_SI;

		for(int mat=0;mat<mp::material.size();mat++){
			const double Hsinwt_local=mp::material[mat].fmr_field_strength*sin(2.0*M_PI*real_time*mp::material[mat].fmr_field_frequency);
			stream << Hsinwt_local << "\t";
		}
	}

	// Output Function 26
	void mat_mdoth(std::ostream& stream){
		// initialise vector of H
		std::vector<double> H(&sim::H_vec[0], &sim::H_vec[0]+3);
		stream << stats::material_magnetization.output_normalized_magnetization_dot_product(H);
	}

	// Output Function 27
	void total_energy(std::ostream& stream){
      stream << stats::system_energy.output_energy(stats::total);
	}

	// Output Function 28
	void mean_total_energy(std::ostream& stream){
      stream << stats::system_energy.output_mean_energy(stats::total);
	}

	// Output Function 29
	void total_anisotropy_energy(std::ostream& stream){
      stream << stats::system_energy.output_energy(stats::anisotropy);
	}

	// Output Function 30
	void mean_total_anisotropy_energy(std::ostream& stream){
      stream << stats::system_energy.output_mean_energy(stats::anisotropy);
	}

	// Output Function 31
	/*void total_cubic_anisotropy_energy(std::ostream& stream){
		stats::output_energy(stream, stats::cubic_anisotropy, stats::total);
	}

	// Output Function 32
	void mean_total_cubic_anisotropy_energy(std::ostream& stream){
		stats::output_energy(stream, stats::cubic_anisotropy, stats::mean);
	}

	// Output Function 33
	void total_surface_anisotropy_energy(std::ostream& stream){
		stats::output_energy(stream, stats::surface_anisotropy, stats::total);
	}

	// Output Function 34
	void mean_total_surface_anisotropy_energy(std::ostream& stream){
		stats::output_energy(stream, stats::surface_anisotropy, stats::mean);
	}*/

	// Output Function 35
	void total_exchange_energy(std::ostream& stream){
      stream << stats::system_energy.output_energy(stats::exchange);
	}

	// Output Function 36
	void mean_total_exchange_energy(std::ostream& stream){
      stream << stats::system_energy.output_mean_energy(stats::exchange);
	}

	// Output Function 37
	void total_applied_field_energy(std::ostream& stream){
      stream << stats::system_energy.output_energy(stats::applied_field);
	}

	// Output Function 38
	void mean_total_applied_field_energy(std::ostream& stream){
      stream << stats::system_energy.output_mean_energy(stats::applied_field);
	}

	// Output Function 39
	void total_magnetostatic_energy(std::ostream& stream){
      stream << stats::system_energy.output_energy(stats::magnetostatic);
	}

	// Output Function 40
	void mean_total_magnetostatic_energy(std::ostream& stream){
      stream << stats::system_energy.output_mean_energy(stats::magnetostatic);
	}

	// Output Function 41
	/*void total_so_anisotropy_energy(std::ostream& stream){
		stats::output_energy(stream, stats::second_order_anisotropy, stats::total);
	}

	// Output Function 42
	void mean_total_so_anisotropy_energy(std::ostream& stream){
		stats::output_energy(stream, stats::second_order_anisotropy, stats::mean);
	}*/

	// Output Function 43
	void height_mvec(std::ostream& stream){
		stream << stats::height_magnetization.output_normalized_magnetization();
	}

	// Output Function 44
	void material_height_mvec(std::ostream& stream){
		stream << stats::material_height_magnetization.output_normalized_magnetization();
	}

	// Output Function 45
	void height_mvec_actual(std::ostream& stream){
		stream << stats::height_magnetization.output_magnetization();
	}

	// Output Function 46
	void material_height_mvec_actual(std::ostream& stream){
		stream << stats::material_height_magnetization.output_magnetization();
	}

	// Output Function 47
	void fmr_field_strength(std::ostream& stream){
		stream << sim::fmr_field << "\t";
	}

   // Output Function 48
	void mean_mvec(std::ostream& stream){
		stream << stats::system_magnetization.output_normalized_mean_magnetization();
	}

   // Output Function 49
	void mat_mean_mvec(std::ostream& stream){
		stream << stats::material_magnetization.output_normalized_mean_magnetization();
	}

   // Output Function 50
   void mean_material_susceptibility(std::ostream& stream){
		stream << stats::material_susceptibility.output_mean_susceptibility(sim::temperature);
	}

	// Output Function 51
   void mean_height_magnetisation_length(std::ostream& stream){
		stream << stats::height_magnetization.output_mean_magnetization_length();
	}

	// Output Function 51
   void mean_height_magnetisation(std::ostream& stream){
		stream << stats::height_magnetization.output_mean_magnetization();
	}

	// Output Function 60
	void MPITimings(std::ostream& stream){
		stream << vmpi::AverageComputeTime+vmpi::AverageWaitTime << "\t" << vmpi::AverageComputeTime << "\t" << vmpi::AverageWaitTime;
		stream << "\t" << vmpi::MaximumComputeTime << "\t" << vmpi::MaximumWaitTime << "\t";
	}

   // Output Function 61
   void mean_system_specific_heat(std::ostream& stream){
      stream << stats::system_specific_heat.output_mean_specific_heat(sim::temperature);
   }

   // Output Function 62
   void mean_material_specific_heat(std::ostream& stream){
      stream << stats::material_specific_heat.output_mean_specific_heat(sim::temperature);
   }

   // Output Function 63
   void material_total_energy(std::ostream& stream){
      stream << stats::material_energy.output_energy(stats::total);
   }

   // Output Function 64
   void material_mean_total_energy(std::ostream& stream){
      stream << stats::material_energy.output_mean_energy(stats::total);
   }

   // Output Function 65
   void resistance(std::ostream& stream){
      stream << spin_transport::total_resistance << "\t";
   }

   // Output Function 66
   void current(std::ostream& stream){
      stream << spin_transport::total_current << "\t";
   }

   // Output 67 reserved for voltage
}
