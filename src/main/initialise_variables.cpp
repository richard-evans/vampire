///
/// @file
/// @brief This is the brief (one line only) description of the function of this file. 
///
/// @details This is the detailed description of the funtion of this file
///
/// @section notes Implementation Notes
/// This is a list of other notes, not related to functionality but rather to implementation. 
/// Also include references, formulae and other notes here.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    11/01/2010
/// @internal
///	Created:		11/01/2010
///	Revision:	  ---
///=====================================================================================
///

// Headers
#include "errors.hpp"
#include "demag.hpp"
#include "voronoi.hpp"
#include "material.hpp"
//#include "multilayer.hpp"
#include "sim.hpp"
#include "random.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

#include <cmath>
#include <iostream>
//==========================================================
// Namespace material_parameters
//==========================================================
namespace mp{
	//----------------------------------
	// Material Container
	//----------------------------------

	//const int max_materials=100;

	int num_materials=1;


	std::vector <materials_t> material(1);


	
	//----------------------------------
	//Input Integration parameters
	//----------------------------------
	double dt_SI;
	double gamma_SI = 1.76E11;
	
	//----------------------------------
	//Derived Integration parameters
	//----------------------------------
	double dt;
	double half_dt;
	
	// Unrolled material parameters for speed
	std::vector <double> MaterialMuSSIArray(0);
	std::vector <zkval_t> MaterialScalarAnisotropyArray(0);
	std::vector <zkten_t> MaterialTensorAnisotropyArray(0);
	std::vector <double> MaterialCubicAnisotropyArray(0);

///
/// @brief Function to initialise program variables prior to system creation.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    19/01/2010
///
/// @param[in] infile Main input file name for system initialisation
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		19/01/2010
///	Revision:	  ---
///=====================================================================================
///
int initialise(std::string const infile){
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "initialise_variables has been called" << std::endl;}

	if(vmpi::my_rank==0){
		std::cout << "================================================================================" << std::endl;
		std::cout << " " << std::endl;
		std::cout << "Initialising system variables" << std::endl;
	}
	
	// Setup default system settings
	mp::default_system();
	
	// Read values from input files
	int iostat = vin::read(infile);
	if(iostat==EXIT_FAILURE){
		std::cerr << "Error - input file \'" << infile << "\' not found, exiting" << std::endl;
		err::vexit();
	}
	
	// Print out material properties
	//mp::material[0].print();

	// Check for keyword parameter overide
	if(cs::single_spin==true){
		mp::single_spin_system();
	}
	
	// Set derived system parameters
	mp::set_derived_parameters();
	
	// Return
	return EXIT_SUCCESS;
}

int default_system(){

	// Initialise system creation flags to zero
	for (int i=0;i<10;i++){
		cs::system_creation_flags[i] = 0;
		sim::hamiltonian_simulation_flags[i] = 0; 
	}
	
	// Set system dimensions !Angstroms
	cs::unit_cell_size[0] = 3.0;
	cs::unit_cell_size[1] = 3.0;
	cs::unit_cell_size[2] = 3.0;

	cs::system_dimensions[0] = 100.0;
	cs::system_dimensions[1] = 100.0;
	cs::system_dimensions[2] = 100.0;

	cs::particle_scale   = 50.0;
	cs::particle_spacing = 10.0;
	
	cs::particle_creation_parity=0;
	cs::crystal_structure = "sc";

	// Voronoi Variables
	create_voronoi::voronoi_sd=0.1;
	create_voronoi::parity=0;
	
	// Setup Hamiltonian Flags
	sim::hamiltonian_simulation_flags[0] = 1;	// Exchange
	sim::hamiltonian_simulation_flags[1] = 1;	// Anisotropy
	sim::hamiltonian_simulation_flags[2] = 1;	// Applied
	sim::hamiltonian_simulation_flags[3] = 1;	// Thermal
	sim::hamiltonian_simulation_flags[4] = 0;	// Dipolar

	// Setup Simulation Variables
	sim::total_time = 1000000;			// total simulation time (single run)
	sim::loop_time = 0;			// time in loop, eg hysteresis, Tc
	sim::partial_time=100;			// time between statistics collection
	sim::equilibration_time=100000;	// time for equilibration before main loop
	sim::temperature = 0.0;	// Constant system temperature

	// demag variables
	//demag::demag_resolution=2;
	//demag::update_rate=10000;
	
	//Integration parameters
	dt_SI = 1.0e-15;	// seconds
	dt = dt_SI*mp::gamma_SI; // Must be set before Hth
	half_dt = 0.5*dt;

	// MPI Mode (Assume decomposition)
	//vmpi::mpi_mode=1;
	//mpi_create_variables::mpi_interaction_range=2; // Unit cells
	//mpi_create_variables::mpi_comms_identify=true;
	
	//------------------------------------------------------------------------------
	// Material Definitions
	//------------------------------------------------------------------------------
	num_materials=1;
	material.resize(num_materials);

	//-------------------------------------------------------
	// Material 0
	//-------------------------------------------------------
	material[0].name="Co";
	material[0].alpha=0.1;
	material[0].Jij_matrix_SI[0]=-11.2e-21;
	material[0].mu_s_SI=1.5*9.27400915e-24;
	material[0].Ku1_SI=-4.644e-24;
	material[0].gamma_rel=1.0;
	material[0].element="Ag ";

	// Disable Error Checking
	err::check=false;
	
	// Initialise random number generator
	mtrandom::grnd.seed(1234);

	return EXIT_SUCCESS;
}

int single_spin_system(){

	// Reset system creation flags to zero
	for (int i=0;i<10;i++){
		cs::system_creation_flags[i] = 0;
	}
	
	// Set system dimensions !Angstroms
	cs::unit_cell_size[0] = 3.0;
	cs::unit_cell_size[1] = 3.0;
	cs::unit_cell_size[2] = 3.0;

	cs::system_dimensions[0] = 2.0;
	cs::system_dimensions[1] = 2.0;
	cs::system_dimensions[2] = 2.0;

	cs::particle_scale   = 50.0;
	cs::particle_spacing = 10.0;
	
	cs::particle_creation_parity=0;
	cs::crystal_structure = "sc";
	
	// Turn off multi-spin Flags
	sim::hamiltonian_simulation_flags[0] = 0;	// Exchange
	sim::hamiltonian_simulation_flags[4] = 0;	// Dipolar

	// MPI Mode (Homogeneous execution)
	//vmpi::mpi_mode=0;
	//mpi_create_variables::mpi_interaction_range=2; // Unit cells
	//mpi_create_variables::mpi_comms_identify=false;

	return EXIT_SUCCESS;
}

int set_derived_parameters(){
		
	// Set integration constants
	mp::dt = mp::dt_SI*mp::gamma_SI; // Must be set before Hth
	mp::half_dt = 0.5*mp::dt;

	// Ensure H vector is unit length
	double mod_H=1.0/sqrt(sim::H_vec[0]*sim::H_vec[0]+sim::H_vec[1]*sim::H_vec[1]+sim::H_vec[2]*sim::H_vec[2]);
	sim::H_vec[0]*=mod_H;
	sim::H_vec[1]*=mod_H;
	sim::H_vec[2]*=mod_H;

	// Calculate moment, magnetisation, and anisotropy constants
	/*for(int mat=0;mat<mp::num_materials;mat++){
		double V=cs::unit_cell_size[0]*cs::unit_cell_size[1]*cs::unit_cell_size[2];
		// Set magnetisation from mu_s and a
		if(material[mat].moment_flag==true){
			//material[mat].magnetisation=num_atoms_per_unit_cell*material[mat].mu_s_SI/V;
		}
		// Set mu_s from magnetisation and a
		else {
			//material[mat].mu_s_SI=material[mat].magnetisation*V/num_atoms_per_unit_cell;
		}
		// Set K as energy/atom
		if(material[mat].anis_flag==false){
			material[mat].Ku1_SI=material[mat].Ku1_SI*V/num_atoms_per_unit_cell;
			std::cout << "setting " << material[mat].Ku1_SI << std::endl;
		}
	}*/
	const string blank="";
	// Set derived material parameters
	for(int mat=0;mat<mp::num_materials;mat++){
		mp::material[mat].one_oneplusalpha_sq			=-mp::material[mat].gamma_rel/(1.0+mp::material[mat].alpha*mp::material[mat].alpha);
		mp::material[mat].alpha_oneplusalpha_sq			= mp::material[mat].alpha*mp::material[mat].one_oneplusalpha_sq;
		
		// set initial spins to unit length
		double sx = mp::material[mat].initial_spin[0];
		double sy = mp::material[mat].initial_spin[1];
		double sz = mp::material[mat].initial_spin[2];

		double modS = 1.0/sqrt(sx*sx+sy*sy+sz*sz);
		mp::material[mat].initial_spin[0]*=modS;
		mp::material[mat].initial_spin[1]*=modS;
		mp::material[mat].initial_spin[2]*=modS;
			
		for(int j=0;j<mp::num_materials;j++){
			material[mat].Jij_matrix[j]				= mp::material[mat].Jij_matrix_SI[j]/mp::material[mat].mu_s_SI;
		}
		mp::material[mat].Ku									= mp::material[mat].Ku1_SI/mp::material[mat].mu_s_SI;
		mp::material[mat].Kc									= mp::material[mat].Kc1_SI/mp::material[mat].mu_s_SI;
		mp::material[mat].Ks									= mp::material[mat].Ks_SI/mp::material[mat].mu_s_SI;
		mp::material[mat].H_th_sigma						= sqrt(2.0*mp::material[mat].alpha*1.3806503e-23/
																  (mp::material[mat].mu_s_SI*mp::material[mat].gamma_rel*dt));
	}
		// Check for which anisotropy function(s) are to be used		
		if(sim::TensorAnisotropy==true){
			sim::UniaxialScalarAnisotropy=false; // turn off scalar anisotropy calculation
			// loop over materials and convert all scalar anisotropy to tensor (along z)
			for(int mat=0;mat<mp::num_materials; mat++){
				
				const double one_o_mu=1.0/mp::material[mat].mu_s_SI;

				// If tensor is unset
				if(mp::material.at(mat).KuVec_SI.size()==0){
					const double ex = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(0);
					const double ey = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(1);
					const double ez = mp::material.at(mat).UniaxialAnisotropyUnitVector.at(2);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ex*ex);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ex*ey);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ex*ez);

					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ey*ex);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ey*ey);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ey*ez);

					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ez*ex);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ez*ey);
					mp::material.at(mat).KuVec.push_back(mp::material[mat].Ku*ez*ez);
				}
				else if(mp::material.at(mat).KuVec_SI.size()==9){
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(0)*one_o_mu);
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(1)*one_o_mu);
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(2)*one_o_mu);

					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(3)*one_o_mu);
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(4)*one_o_mu);
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(5)*one_o_mu);

					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(6)*one_o_mu);
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(7)*one_o_mu);
					mp::material.at(mat).KuVec.push_back(mp::material.at(mat).KuVec_SI.at(8)*one_o_mu);
				}
			}
		}
		
		// Unroll anisotropy values for speed
		if(sim::UniaxialScalarAnisotropy==true){
			zlog << zTs() << "Setting scalar uniaxial anisotropy." << std::endl;
			// Set global anisotropy type
			sim::AnisotropyType=0;
			MaterialScalarAnisotropyArray.resize(mp::num_materials);
			for(int mat=0;mat<mp::num_materials; mat++) MaterialScalarAnisotropyArray[mat].K=mp::material[mat].Ku;
		}
		else if(sim::TensorAnisotropy==true){
			zlog << zTs() << "Setting tensor uniaxial anisotropy." << std::endl;
			// Set global anisotropy type
			sim::AnisotropyType=1;
			MaterialTensorAnisotropyArray.resize(mp::num_materials);
			for(int mat=0;mat<mp::num_materials; mat++){
				MaterialTensorAnisotropyArray[mat].K[0][0]=mp::material.at(mat).KuVec.at(0);
				MaterialTensorAnisotropyArray[mat].K[0][1]=mp::material.at(mat).KuVec.at(1);
				MaterialTensorAnisotropyArray[mat].K[0][2]=mp::material.at(mat).KuVec.at(2);

				MaterialTensorAnisotropyArray[mat].K[1][0]=mp::material.at(mat).KuVec.at(3);
				MaterialTensorAnisotropyArray[mat].K[1][1]=mp::material.at(mat).KuVec.at(4);
				MaterialTensorAnisotropyArray[mat].K[1][2]=mp::material.at(mat).KuVec.at(5);

				MaterialTensorAnisotropyArray[mat].K[2][0]=mp::material.at(mat).KuVec.at(6);
				MaterialTensorAnisotropyArray[mat].K[2][1]=mp::material.at(mat).KuVec.at(7);
				MaterialTensorAnisotropyArray[mat].K[2][2]=mp::material.at(mat).KuVec.at(8);

			}
		}
		// Unroll cubic anisotropy values for speed
		if(sim::CubicScalarAnisotropy==true){
			zlog << zTs() << "Setting scalar cubic anisotropy." << std::endl;
			MaterialCubicAnisotropyArray.resize(mp::num_materials);
			for(int mat=0;mat<mp::num_materials; mat++) MaterialCubicAnisotropyArray.at(mat)=mp::material[mat].Kc;
		}

		for(int mat=0;mat<mp::num_materials;mat++){
		//std::cout << "checking range exclusivity" << std::endl;
		// Check for exclusivity of range
		if(material[mat].geometry!=0){
			const double lmin=material[mat].min;
			const double lmax=material[mat].max;
			for(int nmat=0;nmat<mp::num_materials;nmat++){
				if(nmat!=mat){
					double min=material[nmat].min;
					double max=material[nmat].max;
					std::cout << lmin << "\t" << min << "\t" << max << std::endl;
					std::cout << lmax << "\t" << min << "\t" << max << std::endl;
					if(((lmin>min) && (lmin<max)) || ((lmax>min) && (lmax<max))){
						std::cerr << "Warning - material " << mat << " overlaps material " << nmat << " - possibly use alloy keyword instead" << std::endl;
						std::cerr << " Material "<< mat << ":min = " << lmin << std::endl;
						std::cerr << " Material "<< mat << ":max = " << lmax << std::endl;
						std::cerr << " Material "<< nmat << ":min = " << min << std::endl;
						std::cerr << " Material "<< nmat << ":max = " << max << std::endl;
						//exit(EXIT_FAILURE);
					}
				}
			}
		}
	}
	
	return EXIT_SUCCESS;
}

} // end of namespace mp
