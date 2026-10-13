//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrew Naden, Richard F L Evans and Rory Pond 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk and rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <sstream>

// Vampire headers
#include "config.hpp"
#include "errors.hpp"
#include "info.hpp"
#include "gpu.hpp"
#include "grains.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "vio.hpp"
#include "vutil.hpp"
#include "micromagnetic.hpp"

// vio module headers
#include "internal.hpp"

///-------------------------------------------------------
/// Function to write header information about simulation
///-------------------------------------------------------
void write_output_file_header(std::ofstream& ofile, std::vector<vout::output_t>& file_output_list){

	//------------------------------------
	// Determine current time
	//------------------------------------
	time_t seconds;

	// get time now
	seconds = time (NULL);
	struct tm * timeinfo;
	char oftime [80];

	timeinfo = localtime ( &seconds );
	// format time string
	strftime (oftime,80,"%Y-%m-%d %X ",timeinfo);

	//------------------------------------
	// Determine current directory
	//------------------------------------
	char directory [256];

	#ifdef WIN_COMPILE
		if(_getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog." << std::endl;
        }
	#else
		if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog." << std::endl;
        }
	#endif

	//------------------------------------
	// Output output file header
	//------------------------------------
	ofile << "#----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	ofile << "# " << "Output file for vampire simulation" << std::endl;
	ofile << "# " << "  time       : " << oftime << "    process id : " << vout::zLogPid << std::endl;
	ofile << "# " << "  hostname   : " << vout::zLogHostName << std::endl;
	ofile << "# " << "  path       : " << directory << std::endl;
   ofile << "# " << "  version    : " << vinfo::version() << std::endl;
   ofile << "# " << "  githash    : " << vinfo::githash() << std::endl;
	ofile << "#----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	//ofile << "# time" << "\t" << "temperature" << "\t" <<  "|m|" << "\t" << "..." << std::endl; // to be concluded...
    if(vout::header_option){
        vout::write_out(ofile,file_output_list);
    }
	return;

}

namespace vout{

   void output_switch(std::ostream& stream,output_t id,bool header){
      //stream.precision(vout::precision);
      switch(id){
      	case output_t::time_steps:
      		vout::time(stream,header);
      		break;
      	case output_t::real_time:
      		vout::real_time(stream,header);
      		break;
      	case output_t::temperature:
      		vout::temperature(stream,header);
      		break;
      	case output_t::applied_field_strength:
      		vout::Happ(stream,header);
      		break;
      	case output_t::applied_field_unit_vector:
      		vout::Hvec(stream,header);
      		break;
      	case output_t::magnetisation:
      		vout::mvec(stream,header);
      		break;
      	case output_t::magnetisation_length:
      		vout::magm(stream,header);
      		break;
      	case output_t::mean_magnetisation_length:
      		vout::mean_magm(stream,header);
      		break;
      	case output_t::material_magnetisation:
      		vout::mat_mvec(stream,header);
      		break;
      	case output_t::material_mean_magnetisation_length:
      		vout::mat_mean_magm(stream,header);
      		break;
      	case output_t::applied_field_alignment:
      		vout::mdoth(stream,header);
      		break;
      	case output_t::total_torque:
      		vout::systorque(stream, header);
      		break;
      	case output_t::mean_total_torque:
      		vout::mean_systorque(stream, header);
      		break;
      	case output_t::constraint_phi:
      		vout::constraint_phi(stream,header);
      		break;
      	case output_t::constraint_theta:
      		vout::constraint_theta(stream,header);
      		break;
      	case output_t::material_constraint_phi:
      		vout::material_constraint_phi(stream,header);
      		break;
      	case output_t::material_constraint_theta:
      		vout::material_constraint_theta(stream,header);
      		break;
      	case output_t::material_mean_torque:
      		vout::material_mean_systorque(stream, header);
      		break;
      	case output_t::mean_susceptibility:
      		vout::mean_system_susceptibility(stream,header);
      		break;
      	case output_t::phonon_temperature:
      		vout::phonon_temperature(stream,header);
      		break;
      	case output_t::material_temperature:
      		vout::material_temperature(stream,header);
      		break;
      	case output_t::material_applied_field_strength:
      		vout::material_applied_field_strength(stream,header);
      		break;
      	case output_t::material_fmr_field_strength:
      		vout::material_fmr_field_strength(stream,header);
      		break;
      	case output_t::material_applied_field_alignment:
      		vout::mat_mdoth(stream,header);
      		break;
      	case output_t::total_energy:
      		vout::total_energy(stream,header);
      		break;
      	case output_t::mean_total_energy:
      		vout::mean_total_energy(stream,header);
      		break;
      	case output_t::anisotropy_energy:
      		vout::total_anisotropy_energy(stream,header);
      		break;
      	case output_t::mean_anisotropy_energy:
      		vout::mean_total_anisotropy_energy(stream,header);
      		break;
      	case output_t::material_torque:
      		vout::material_torque(stream, header);
      		break;
      	case output_t::exchange_energy:
      		vout::total_exchange_energy(stream,header);
      		break;
      	case output_t::mean_exchange_energy:
      		vout::mean_total_exchange_energy(stream,header);
      		break;
      	case output_t::applied_field_energy:
      		vout::total_applied_field_energy(stream,header);
      		break;
      	case output_t::mean_applied_field_energy:
      		vout::mean_total_applied_field_energy(stream,header);
      		break;
      	case output_t::magnetostatic_energy:
      		vout::total_magnetostatic_energy(stream,header);
      		break;
      	case output_t::mean_magnetostatic_energy:
      		vout::mean_total_magnetostatic_energy(stream,header);
      		break;
      	case output_t::dmi_energy:
      		vout::total_dmi_energy(stream,header);
      		break;
      	case output_t::mean_dmi_energy:
      		vout::mean_total_dmi_energy(stream,header);
      		break;
      	case output_t::height_magnetisation_normalised:
      		vout::height_mvec(stream,header);
      		break;
      	case output_t::material_height_magnetisation_normalised:
      		vout::material_height_mvec(stream,header);
      		break;
      	case output_t::height_magnetisation:
      		vout::height_mvec_actual(stream,header);
      		break;
      	case output_t::material_height_magnetisation:
      		vout::material_height_mvec_actual(stream,header);
      		break;
      	case output_t::fmr_field_strength:
      		vout::fmr_field_strength(stream,header);
      		break;
         case output_t::mean_magnetisation:
      		vout::mean_mvec(stream,header);
      		break;
         case output_t::material_mean_magnetisation:
      		vout::mat_mean_mvec(stream,header);
      		break;
         case output_t::material_mean_susceptibility:
      		vout::mean_material_susceptibility(stream,header);
      		break;
      	case output_t::mean_height_magnetisation_length:
      		vout::mean_height_magnetisation_length(stream,header);
      		break;
      	case output_t::mean_height_magnetisation:
      		vout::mean_height_magnetisation(stream,header);
      		break;
      	case output_t::mpi_timings:
      		vout::MPITimings(stream,header);
      		break;
      	case output_t::mean_specific_heat:
      		vout::mean_system_specific_heat(stream,header);
      		break;
      	case output_t::material_mean_specific_heat:
      		vout::mean_material_specific_heat(stream,header);
      		break;
			case output_t::material_total_energy:
				vout::material_total_energy(stream,header);
				break;
			case output_t::material_mean_total_energy:
				vout::material_mean_total_energy(stream,header);
				break;
			case output_t::material_dmi_energy:
				vout::material_dmi_energy(stream,header);
				break;
			case output_t::material_mean_dmi_energy:
				vout::material_mean_dmi_energy(stream,header);
				break;
			case output_t::resistance:
				vout::resistance(stream, header);
			  	break;
			case output_t::current:
			  	vout::current(stream, header);
			  	break;
			case output_t::domain_wall_centre:
				vout::domain_wall_position(stream,header);
				break;
			case output_t::magneto_resistance:
				vout::MRresistance(stream,header);
				break;
			case output_t::lfa_ms:
				vout::lfa_ms(stream,header);
				break;
			case output_t::track_pos_x:
				vout::x_track_pos(stream,header);
				break;
			case output_t::track_pos_z:
				vout::z_track_pos(stream,header);
				break;
			case output_t::fractional_electric_field_strength:
			   vout::fractional_electric_field_strength(stream, header);
   			break;
			case output_t::spin_temperature:
				vout::sysspintemp(stream, header);
				break;
			case output_t::mean_spin_temperature:
				vout::mean_sysspintemp(stream, header);
				break;
			case output_t::material_mean_spin_temperature:
				vout::material_mean_sysspintemp(stream, header);
				break;
			case output_t::material_spin_temperature:
				vout::material_spin_temp(stream, header);
				break;
			case output_t::mean_spin_length:
				vout::mean_system_spin_length(stream,header);
				break;
			case output_t::material_mean_spin_length:
				vout::mean_material_spin_length(stream,header);
				break;
			case output_t::mean_height_spin_length:
				vout::mean_height_spin_length(stream,header);
				break;
			case output_t::potential_energy:
				vout::potential_energy(stream, header);
				break;
			case output_t::kinetic_energy:
				vout::kinetic_energy(stream, header);
				break;
			case output_t::sld_exchange_energy:
				vout::sld_exchange_energy(stream, header);
				break;
			case output_t::sld_coupling_energy:
				vout::sld_coupling_energy(stream, header);
				break;
			case output_t::sld_total_energy:
				vout::sld_total_energy(stream, header);
				break;
			case output_t::lattice_temperature:
				vout::syslatticetemp(stream, header);
				break;
			case output_t::mean_lattice_temperature:
				vout::mean_syslatticetemp(stream, header);
				break;
			case output_t::material_lattice_temperature:
				vout::material_lattice_temp(stream, header);
				break;
			case output_t::material_mean_lattice_temperature:
			   vout::material_mean_syslatticetemp(stream, header);
			   break;
			case output_t::mean_coupling_field_vec:
				vout::mean_coupling_field_vec(stream,header);
				break;
			case output_t::material_binder_cumulant:
				vout::material_binder_cumulant(stream,header);
				break;
   		case output_t::binder_cumulant:
      		vout::system_binder_cumulant(stream,header);
      		break;
			case output_t::material_standard_deviation:
				vout::standard_deviation(stream,header);
				break;

      }

      return;

   }

   void write_out(std::ostream& stream,std::vector<output_t>& list){
      static bool header = true;
      // if header is false then header_option is never checked.
      if(header && !vout::header_option){
         header = false;
      };
      // Output data to output
      if(vmpi::my_rank==0){

         // For gpu acceleration get statistics from device
         if(gpu::acceleration) gpu::stats::get();

         for(unsigned int item=0;item<list.size();item++){
            output_switch(stream,list[item],header);
         }
         // Carriage return
         if(list.size()>0) stream << std::endl;

      } // end of code for rank 0 only
      header = false;
   }

   //-------------------------------------------
	// Data output wrapper function
   //-------------------------------------------
	void data(){


      //if(micromagnetic::discretisation_type != 1){
		// check calling of routine if error checking is activated
		if(err::check==true){std::cout << "vout::data has been called" << std::endl;}

		// Calculate MPI Timings since last data output
		#ifdef MPICF
		if(vmpi::DetailedMPITiming){

			// Calculate Average times
			MPI_Reduce (&vmpi::TotalComputeTime,&vmpi::AverageComputeTime,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce (&vmpi::TotalWaitTime,&vmpi::AverageWaitTime,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			vmpi::AverageComputeTime/=double(vmpi::num_processors);
			vmpi::AverageWaitTime/=double(vmpi::num_processors);

			// Calculate Maximum times
			MPI_Reduce (&vmpi::TotalComputeTime,&vmpi::MaximumComputeTime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
			MPI_Reduce (&vmpi::TotalWaitTime,&vmpi::MaximumWaitTime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

			// Save times for timing matrix
			vmpi::ComputeTimeArray.push_back(vmpi::TotalComputeTime);
			vmpi::WaitTimeArray.push_back(vmpi::TotalWaitTime);

			// reset until next data output
			vmpi::TotalComputeTime=0.0;
			vmpi::TotalWaitTime=0.0;
		}
		#endif

      // check for open ofstream on root process only
      if(vmpi::my_rank == 0){
         if(!zmag.is_open()){
            // check for checkpoint continue and append data
            if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag) zmag.open(vout::output_file_name,std::ofstream::app);
            // otherwise overwrite file
            else{
               zmag.open(vout::output_file_name,std::ofstream::trunc);
               // write file header information
               write_output_file_header(zmag, file_output_list);
            }
         }
      }

      // Only output 1/output_rate time steps// This is all serialised inside the write_output fn - AJN
      if(sim::time%vout::output_rate==0){
         write_out(zmag,file_output_list);
      } // end of if statement for output rate

      if(sim::time%vout::output_rate==0){ // needs to be altered to separate variable at some point
         write_out(std::cout,screen_output_list);
      } // End of if statement to output data to screen

		vout::write_grain_file();

		// Output configuration files to disk
		config::output();

		// optionally save checkpoint file
		if(sim::save_checkpoint_flag==true && sim::save_checkpoint_continuous_flag==true && sim::time%sim::save_checkpoint_rate==0) save_checkpoint();
     // }
      if (micromagnetic::discretisation_type ==1){
         micromagnetic::outputs();
      }
      return;

   } // end of data()

} // end of namespace vout
