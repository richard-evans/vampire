//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrew J Naden, Richard F L Evans and Rory Pond 2016-2019.
//       All rights reserved.
//
//   Email: richard.evans@york.ac.uk and rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "vio.hpp"
#include "material.hpp"
#include "montecarlo.hpp"
#include "program.hpp"
#include "stats.hpp"
#include "sim.hpp"
#include "micromagnetic.hpp"
#include "spintransport.hpp"

// vio module headers
#include "internal.hpp"

namespace vout{
	// Output Function 0
   std::string generic_output_int(std::string str,uint64_t i, bool header){
      std::ostringstream res;
      vout::fixed_width_output result(res, vout::fw_size_int);
      if(header){
         result << str;
      }
      else{
         result << i;
      }
      return result.str();
    }
    //--------------------------------------------------------------------------
    // Function to format a standard double variable including a file header
    //--------------------------------------------------------------------------
    std::string generic_output_double(const std::string str, const double d, const bool header){

      // resuting output stringstream
      std::ostringstream res;

      // set up fixed width if defined
      vout::fixed_width_output result(res,vout::fw_size);

      // switch between header text and actual variable
      if(header) result << str;
      else result << d;

      // return  resulting sstream as string
      return result.str();

    }

    // why don't we do the fixed width stuff here? Yep - TBD
    void time(std::ostream& stream, bool header){
       stream << generic_output_int("Time_steps", sim::time,header);
   }

   // Output Function 1 - with Header
   void real_time(std::ostream& stream, bool header){
      stream << generic_output_double("Real_time",sim::time*mp::dt_SI,header);
   }

   // Output Function 2 - with Header
   void temperature(std::ostream& stream, bool header){
      stream << generic_output_double("Temperature" ,sim::temperature,header);
   }

   // Output Function 3 - with Header
   void Happ(std::ostream& stream, bool header){
      stream << generic_output_double("B_applied" ,sim::H_applied,header);
   }

   // Output Function 4 - with Header
   void Hvec(std::ostream& stream, bool header){
      std::ostringstream res;
      vout::fixed_width_output result(res,vout::fw_size);
      if(header) result << "B_vector_x" << "B_vector_y" << "B_vector_z";
      else result << sim::H_vec[0] << sim::H_vec[1] << sim::H_vec[2];
      stream << result.str();
   }

   // Output Function 5 - with Header
   void mvec(std::ostream& stream, bool header){
      // header included in function
      stream << stats::system_magnetization.output_normalized_magnetization(header);
   }

   // Output Function 6 - with Header
   void magm(std::ostream& stream, bool header){
      // header included in function
      stream << stats::system_magnetization.output_normalized_magnetization_length(header);
   }

   // Output Function 7 - with Header
   void mean_magm(std::ostream& stream, bool header){
      // header included in function
      stream << stats::system_magnetization.output_normalized_mean_magnetization_length(header);
   }

   // Output Function 8 - with Header
   void mat_mvec(std::ostream& stream, bool header){
      // header included in function
      stream << stats::material_magnetization.output_normalized_magnetization(header);
   }

   // Output Function 9 - with Header
   void mat_mean_magm(std::ostream& stream, bool header){
      // header included in function
      stream << stats::material_magnetization.output_normalized_mean_magnetization_length(header);
   }

   // Output Function 12 - with Header
   void mdoth(std::ostream& stream, bool header){
      // initialise vector of H
      std::vector<double> H(&sim::H_vec[0], &sim::H_vec[0]+3);
      stream << stats::system_magnetization.output_normalized_magnetization_dot_product(H,header);
   }

   // Output Function 14 - with Header
   void systorque(std::ostream& stream, bool header){
      stream << stats::system_torque.output_torque(header);
   }

   // Output Function 15 - with Header
   void mean_systorque(std::ostream& stream, bool header){
      stream << stats::system_torque.output_mean_torque(header);
   }

   // Output Function 16 - with Header
   void constraint_phi(std::ostream& stream, bool header){
      stream << generic_output_double("Con_phi",sim::constraint_phi,header);
   }

   // Output Function 17 - with Header
   void constraint_theta(std::ostream& stream, bool header){
      stream << generic_output_double("Con_theta",sim::constraint_theta,header);
   }

   // Output Function 18 - with Header
   void material_constraint_phi(std::ostream& stream, bool header){
      std::ostringstream res;
      vout::fixed_width_output result(res,vout::fw_size);
      for(int mat=0;mat<mp::num_materials;mat++){
         if(header){
            result << "ID" + std::to_string(mat) + "_Con_phi";
         }
         else{
            result << montecarlo::cmc::cmc_mat[mat].constraint_phi;
         }
      }
      stream << result.str();
   }

   // Output Function 19 - with Header
   void material_constraint_theta(std::ostream& stream, bool header){
      std::ostringstream res;
      vout::fixed_width_output result(res,vout::fw_size);
      for(int mat=0;mat<mp::num_materials;mat++){
         if(header){
            result << "ID" + std::to_string(mat) + "_Con_theta";
         }
         else{
            result << montecarlo::cmc::cmc_mat[mat].constraint_theta;
         }
      }
      stream << result.str();
   }

   // Output Function 20
   void material_mean_systorque(std::ostream& stream, bool header){
      stream << stats::material_torque.output_mean_torque(header);
   }

   void material_torque(std::ostream& stream, bool header){
      stream << stats::material_torque.output_torque(header);
   }

   // Output Function 21 - with Header
   void mean_system_susceptibility(std::ostream& stream, bool header){
      stream << stats::system_susceptibility.output_mean_susceptibility(sim::temperature,header);
   }

    // Output Function 998 - with Header
   void system_binder_cumulant(std::ostream& stream, bool header){
      stream << stats::system_binder_cumulant.output_binder_cumulant(header);
   }

   // Output Function 999 - with Header
   void standard_deviation(std::ostream& stream, bool header){
      stream << stats::material_standard_deviation.output_standard_deviation(header);
   }
   // Output Function 22
   void phonon_temperature(std::ostream& stream, bool header){
     stream << generic_output_double("Phonon_temp",sim::TTTp,header);
   }

   // Output Function 23 - with Header
   void material_temperature(std::ostream& stream, bool header){
      std::ostringstream res;
      vout::fixed_width_output result(res,vout::fw_size);
      for(size_t mat=0;mat<mp::material.size();mat++){
         if(header){
            result << "ID" + std::to_string(mat) + "_Temp";
         }
         else{
            result << mp::material[mat].temperature << "\t";
         }
      }
      stream << result.str();
   }

   // Output Function 24 - with Header
   void material_applied_field_strength(std::ostream& stream, bool header){
      std::ostringstream res;
      vout::fixed_width_output result(res,vout::fw_size);
      for(size_t mat=0;mat<mp::material.size();mat++){
         if(header){
            result << "ID" + std::to_string(mat) + "_H";
         }
         else{
            result << mp::material[mat].applied_field_strength;
         }
      }
      stream << result.str();
   }

   // Output Function 25 - with Header
   void material_fmr_field_strength(std::ostream& stream, bool header){
      std::ostringstream res;
      vout::fixed_width_output result(res,vout::fw_size);

      const double real_time=sim::time*mp::dt_SI;

      for(size_t mat=0;mat<mp::material.size();mat++){
         if(header){
            result << "ID" + std::to_string(mat) + "_fmr_H";
         }
         else{
            const double Hsinwt_local=mp::material[mat].fmr_field_strength*sin(2.0*M_PI*real_time*mp::material[mat].fmr_field_frequency);
            result << Hsinwt_local;
         }
      }
      stream << result.str();
   }

	// Output Function 26 - with Header
	void mat_mdoth(std::ostream& stream, bool header){
		// initialise vector of H
		std::vector<double> H(&sim::H_vec[0], &sim::H_vec[0]+3);
		stream << stats::material_magnetization.output_normalized_magnetization_dot_product(H,header);
	}

	// Output Function 27 - with Header
	void total_energy(std::ostream& stream, bool header){
      stream << stats::system_energy.output_energy(stats::total,header);
	}

	// Output Function 28 - with Header
	void mean_total_energy(std::ostream& stream, bool header){
      stream << stats::system_energy.output_mean_energy(stats::total,header);
	}

	// Output Function 29 - with Header
	void total_anisotropy_energy(std::ostream& stream, bool header){
      stream << stats::system_energy.output_energy(stats::anisotropy,header);
	}

	// Output Function 30 - with Header
	void mean_total_anisotropy_energy(std::ostream& stream, bool header){
      stream << stats::system_energy.output_mean_energy(stats::anisotropy,header);
	}

	// Output Function 31
	/*void total_cubic_anisotropy_energy(std::ostream& stream,bool header){
		stats::output_energy(stream, stats::cubic_anisotropy, stats::total);
	}

	// Output Function 32
	void mean_total_cubic_anisotropy_energy(std::ostream& stream,bool header){
		stats::output_energy(stream, stats::cubic_anisotropy, stats::mean);
	}

	// Output Function 33
	void total_surface_anisotropy_energy(std::ostream& stream,bool header){
		stats::output_energy(stream, stats::surface_anisotropy, stats::total);
	}

	// Output Function 34
	void mean_total_surface_anisotropy_energy(std::ostream& stream,bool header){
		stats::output_energy(stream, stats::surface_anisotropy, stats::mean);
	}*/

	// Output Function 35 - with Header
	void total_exchange_energy(std::ostream& stream, bool header){
      stream << stats::system_energy.output_energy(stats::exchange,header);
	}

	// Output Function 36 - with Header
	void mean_total_exchange_energy(std::ostream& stream, bool header){
      stream << stats::system_energy.output_mean_energy(stats::exchange,header);
	}

	// Output Function 37 - with Header
	void total_applied_field_energy(std::ostream& stream, bool header){
      stream << stats::system_energy.output_energy(stats::applied_field,header);
	}

	// Output Function 38 - with Header
	void mean_total_applied_field_energy(std::ostream& stream, bool header){
      stream << stats::system_energy.output_mean_energy(stats::applied_field,header);
	}

	// Output Function 39 - with Header
	void total_magnetostatic_energy(std::ostream& stream, bool header){
      stream << stats::system_energy.output_energy(stats::magnetostatic,header);
	}

	// Output Function 40 - with Header
	void mean_total_magnetostatic_energy(std::ostream& stream, bool header){
      stream << stats::system_energy.output_mean_energy(stats::magnetostatic,header);
	}

	// Output Function 41
	/*void total_so_anisotropy_energy(std::ostream& stream,bool header){
		stats::output_energy(stream, stats::second_order_anisotropy, stats::total);
	}

	// Output Function 42
	void mean_total_so_anisotropy_energy(std::ostream& stream,bool header){
		stats::output_energy(stream, stats::second_order_anisotropy, stats::mean);
	}*/

	// Output Function 43 - with Header
	void height_mvec(std::ostream& stream, bool header){
		stream << stats::height_magnetization.output_normalized_magnetization(header);
	}

	// Output Function 44 - with Header
	void material_height_mvec(std::ostream& stream, bool header){
		stream << stats::material_height_magnetization.output_normalized_magnetization(header);
	}

	// Output Function 45 - with Header
	void height_mvec_actual(std::ostream& stream, bool header){
		stream << stats::height_magnetization.output_magnetization(header);
	}

	// Output Function 46 - with Header
	void material_height_mvec_actual(std::ostream& stream, bool header){
		stream << stats::material_height_magnetization.output_magnetization(header);
	}

	// Output Function 47 - with Header
	void fmr_field_strength(std::ostream& stream, bool header){
		stream << sim::fmr_field << "\t";
	}

   // Output Function 48 - with Header
	void mean_mvec(std::ostream& stream, bool header){
		stream << stats::system_magnetization.output_normalized_mean_magnetization(header);
	}

   // Output Function 49 - with Header
	void mat_mean_mvec(std::ostream& stream, bool header){
		stream << stats::material_magnetization.output_normalized_mean_magnetization(header);
	}

   // Output Function 50 - with Header
   void mean_material_susceptibility(std::ostream& stream, bool header){
		stream << stats::material_susceptibility.output_mean_susceptibility(sim::temperature,header);
	}

   // Output Function 997 - with Header
   void material_binder_cumulant(std::ostream& stream, bool header){
		stream << stats::material_binder_cumulant.output_binder_cumulant(header);
   }

	// Output Function 51 - with Header
   void mean_height_magnetisation_length(std::ostream& stream, bool header){
		stream << stats::height_magnetization.output_mean_magnetization_length(header);
	}

	// Output Function 51 - with Header
   void mean_height_magnetisation(std::ostream& stream, bool header){
		stream << stats::height_magnetization.output_mean_magnetization(header);
	}

	// Output Function 60
	void MPITimings(std::ostream& stream, bool header){
		stream << vmpi::AverageComputeTime+vmpi::AverageWaitTime << "\t" << vmpi::AverageComputeTime << "\t" << vmpi::AverageWaitTime;
		stream << "\t" << vmpi::MaximumComputeTime << "\t" << vmpi::MaximumWaitTime << "\t";
	}

   // Output Function 61 - with Header
   void mean_system_specific_heat(std::ostream& stream, bool header){
      stream << stats::system_specific_heat.output_mean_specific_heat(sim::temperature,header);
   }

   // Output Function 62 - with Header
   void mean_material_specific_heat(std::ostream& stream, bool header){
      stream << stats::material_specific_heat.output_mean_specific_heat(sim::temperature,header);
   }

   // Output Function 63 - with Header
   void material_total_energy(std::ostream& stream, bool header){
      stream << stats::material_energy.output_energy(stats::total,header);
   }

   // Output Function 64 - with Header
   void material_mean_total_energy(std::ostream& stream, bool header){
      stream << stats::material_energy.output_mean_energy(stats::total,header);
   }

   // Output Function 65
   void resistance(std::ostream& stream, bool header){
      stream << generic_output_double("resistance", spin_transport::total_resistance, header);
   }

   // Output Function 66
   void current(std::ostream& stream, bool header){
      stream << generic_output_double("current", spin_transport::total_current, header);
   }

   // Output Function 67
   void domain_wall_position(std::ostream& stream, bool header){
      stream << sim::domain_wall_centre;
   }

   // Output Function 68
   void MRresistance(std::ostream& stream, bool header){
      if(header){
         stream << "MR" << "\t";
      }
      else{
         stream << micromagnetic::MR_resistance << "\t";
      }
   }

   // Output Function 69
   void lfa_ms(std::ostream& stream, bool header){
      if(header){
         stream << "MS" << "\t";
      }
      else{
         stream << sim::Ms << "\t";
      }
   }

   // Output Function 70
   void x_track_pos(std::ostream& stream, bool header){
      if(header){
         stream << "x-pos" << "\t";
      }
      else{
         stream << sim::track_pos_x << "\t";
      }
   }

   // Output Function 71
   void z_track_pos(std::ostream& stream, bool header){
      if(header){
         stream << "z pos" << "\t";
      }
      else{
         stream << sim::track_pos_z << "\t";
      }
   }

   // Output Function 72
   void fractional_electric_field_strength(std::ostream& stream, bool header){
      stream << generic_output_double("frac_voltage", program::fractional_electric_field_strength, header);
   }

}
