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
#include <iostream>
#include "public.hpp"
#include "demag.hpp"
#include "voronoi.hpp"
#include "material.hpp"
//#include "multilayer.hpp"
#include "sim.hpp"
#include "random.hpp"
#include "vio.hpp"
#include "vmpi.hpp"


//==========================================================
// Namespace material_parameters
//==========================================================
namespace mp{
	//----------------------------------
	// Material Container
	//----------------------------------

	//const int max_materials=100;

	int num_materials=1;


	valarray <materials_t> material(1);


	
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
	
	//----------------------------------
	//Input System Parameters
	//----------------------------------
	int particle_creation_parity;	// Offset of particle centre (odd/even)
	int int_system_dimensions[3];
	double system_dimensions[3];	// Size of system (A)
	double particle_scale;			// Diameter of particles/grains (A)
	double particle_spacing;		// Spacing Between particles (A)
	
	//----------------------------------
	//Derived System Parameters
	//----------------------------------
	double lattice_constant[3];
	double lattice_space_conversion[3];
	string crystal_structure;
	bool single_spin=false;
	string hamiltonian_type; 	// generic
								// LR_FePt
								// SR_FePt
								// LR_Co
								// LR_Fe
								// LR_Ni
								
	string atomic_element[4]; // Different atomic species
	
	int num_nearest_neighbours = 0;
	int hamiltonian_num_neighbours = 0;
	
	//----------------------------------
	// System creation flags
	//----------------------------------
	
	int system_creation_flags[10];
	
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
	if(error_checking::error_check==true){std::cout << "initialise_variables has been called" << std::endl;}

	// Setup default system settings
	mp::default_system();
	
	// Read values from input files
	int iostat = vin::read(infile);
	if(iostat==EXIT_FAILURE){
		std::cerr << "Error - input file \'" << infile << "\' not found, exiting" << std::endl;
		exit(1);
	}
	
	// Check for keyword parameter overide
	if(mp::single_spin==true){
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
		mp::system_creation_flags[i] = 0;
		sim::hamiltonian_simulation_flags[i] = 0; 
	}
	
	// Set system dimensions !Angstroms
	mp::lattice_constant[0] = 3.0;
	mp::lattice_constant[1] = 3.0;
	mp::lattice_constant[2] = 3.0;

	mp::system_dimensions[0] = 100.0;
	mp::system_dimensions[1] = 100.0;
	mp::system_dimensions[2] = 100.0;

	mp::particle_scale   = 50.0;
	mp::particle_spacing = 10.0;
	
	mp::particle_creation_parity=0;
	mp::crystal_structure = "sc";
	mp::hamiltonian_type = "generic";

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
	demag::demag_resolution=2;
	demag::update_rate=10000;
	
	//Integration parameters
	dt_SI = 1.0e-15;	// seconds
	dt = dt_SI*mp::gamma_SI; // Must be set before Hth
	half_dt = 0.5*dt;

	// MPI Mode (Assume decomposition)
	vmpi::mpi_mode=2;
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
	material[0].hamiltonian_type="generic";
	material[0].element="Ag ";

	// Disable Error Checking
	error_checking::error_check=false;
	
	// Initialise random number generator
	mtrandom::grnd.seed(1234);

	return EXIT_SUCCESS;
}

int single_spin_system(){

	// Reset system creation flags to zero
	for (int i=0;i<10;i++){
		mp::system_creation_flags[i] = 0;
	}
	
	// Set system dimensions !Angstroms
	mp::lattice_constant[0] = 3.0;
	mp::lattice_constant[1] = 3.0;
	mp::lattice_constant[2] = 3.0;

	mp::system_dimensions[0] = 2.0;
	mp::system_dimensions[1] = 2.0;
	mp::system_dimensions[2] = 2.0;

	mp::particle_scale   = 50.0;
	mp::particle_spacing = 10.0;
	
	mp::particle_creation_parity=0;
	mp::crystal_structure = "sc";
	mp::hamiltonian_type = "generic";
	
	// Turn off multi-spin Flags
	sim::hamiltonian_simulation_flags[0] = 0;	// Exchange
	sim::hamiltonian_simulation_flags[4] = 0;	// Dipolar

	// MPI Mode (Homogeneous execution)
	vmpi::mpi_mode=0;
	//mpi_create_variables::mpi_interaction_range=2; // Unit cells
	//mpi_create_variables::mpi_comms_identify=false;

	return EXIT_SUCCESS;
}

int set_derived_parameters(){

	//----------------------------------
	//Derived System Parameters
	//----------------------------------
	mp::lattice_space_conversion[0] = mp::lattice_constant[0]*0.5;
	mp::lattice_space_conversion[1] = mp::lattice_constant[1]*0.5*0.333333333333333;
	mp::lattice_space_conversion[2] = mp::lattice_constant[2]*0.5;

	mp::int_system_dimensions[0] = 2*round(mp::system_dimensions[0]/mp::lattice_constant[0]);
	mp::int_system_dimensions[1] = 6*round(mp::system_dimensions[1]/mp::lattice_constant[1]);
	mp::int_system_dimensions[2] = 2*round(mp::system_dimensions[2]/mp::lattice_constant[2]);
		
	double num_atoms_per_unit_cell=0; 
	
	if(mp::crystal_structure=="sc"){
		mp::num_nearest_neighbours = 6;
		num_atoms_per_unit_cell=1.0;
	}
	else if(mp::crystal_structure=="bcc"){
		mp::num_nearest_neighbours = 8;
		num_atoms_per_unit_cell=2.0;
	}
	else if(mp::crystal_structure=="fct"){
		mp::num_nearest_neighbours = 4;
		num_atoms_per_unit_cell=2.0;
	}
	else if(mp::crystal_structure=="fcc"){
		mp::num_nearest_neighbours = 12;
		num_atoms_per_unit_cell=4.0;
	}
	else{
		 std::cout << "Error in determining num_nearest_neighbours - unknown crystal type \'" << mp::crystal_structure << "\'" << std::endl;
		 exit(1);
	}
	
	if(mp::hamiltonian_type=="generic")	mp::hamiltonian_num_neighbours = mp::num_nearest_neighbours;
	if(mp::hamiltonian_num_neighbours==0){
		 std::cout << "Error in determining hamiltonian_num_neighbours - unknown Hamiltonian type \'" << mp::hamiltonian_type << "\'" << std::endl;
		 exit(1);
	}

	// Set integration constants
	mp::dt = mp::dt_SI*mp::gamma_SI; // Must be set before Hth
	mp::half_dt = 0.5*mp::dt;

	// Calculate moment, magnetisation, and anisotropy constants
	for(int mat=0;mat<mp::num_materials;mat++){
		double V=mp::lattice_constant[0]*mp::lattice_constant[1]*mp::lattice_constant[2];
		// Set magnetisation from mu_s and a
		if(material[mat].moment_flag==true){
			material[mat].magnetisation=num_atoms_per_unit_cell*material[mat].mu_s_SI/V;
		}
		// Set mu_s from magnetisation and a
		else {
			material[mat].mu_s_SI=material[mat].magnetisation*V/num_atoms_per_unit_cell;
		}
		// Set K as energy/atom
		if(material[mat].anis_flag==false){
			material[mat].Ku1_SI=material[mat].Ku1_SI*V/num_atoms_per_unit_cell;
			std::cout << "setting " << material[mat].Ku1_SI << std::endl;
		}
	}
	const string blank="";
	// Set derived material parameters
	for(int mat=0;mat<mp::num_materials;mat++){
		mp::material[mat].hamiltonian_type="generic";
		mp::material[mat].one_oneplusalpha_sq			=-mp::material[mat].gamma_rel/(1.0+mp::material[mat].alpha*mp::material[mat].alpha);
		mp::material[mat].alpha_oneplusalpha_sq			= mp::material[mat].alpha*mp::material[mat].one_oneplusalpha_sq;
		for(int j=0;j<mp::num_materials;j++){
			material[mat].Jij_matrix[j]				= mp::material[mat].Jij_matrix_SI[j]/mp::material[mat].mu_s_SI;
		}
		mp::material[mat].Ku									= mp::material[mat].Ku1_SI/mp::material[mat].mu_s_SI;
		mp::material[mat].H_th_sigma						= sqrt(2.0*mp::material[mat].alpha*1.3806503e-23/
																  (mp::material[mat].mu_s_SI*mp::material[mat].gamma_rel*dt));
		// If local crystal is unset, use global type
		if(mp::material[mat].crystal_structure==blank){
			mp::material[mat].crystal_structure=mp::crystal_structure;
		}
		// calculate number of neighbours for each material
		if(mp::material[mat].hamiltonian_type=="generic"){
			if(mp::material[mat].crystal_structure=="sc"){
				mp::material[mat].num_nearest_neighbours		= 6;
				mp::material[mat].hamiltonian_num_neighbours	= 6;
				mp::material[mat].cutoff = 1.01;
			}
			else if(mp::material[mat].crystal_structure=="bcc"){
				mp::material[mat].num_nearest_neighbours		= 8;
				mp::material[mat].hamiltonian_num_neighbours	= 8;
				mp::material[mat].cutoff = sqrt(3.0)*0.5*1.01;

			}
			else if(mp::material[mat].crystal_structure=="fcc"){
				mp::material[mat].num_nearest_neighbours		= 12;
				mp::material[mat].hamiltonian_num_neighbours	= 12;
				mp::material[mat].cutoff = sqrt(2.0)*0.5*1.01;

			}
			else{
				std::cerr << "Error in determining num_nearest_neighbours - unknown crystal type \'";
				std::cerr << mp::material[mat].crystal_structure << "\'" << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		else{
			std::cerr << "Error, only generic hamiltonians are implemented at present, exiting" << std::endl;
			exit(EXIT_FAILURE);
		}
		
		
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
						std::cerr << "Error - material " << mat << " overlaps material " << nmat << " - use alloy keyword instead" << std::endl;
						std::cerr << " Material "<< mat << ":min = " << lmin << std::endl;
						std::cerr << " Material "<< mat << ":max = " << lmax << std::endl;
						std::cerr << " Material "<< nmat << ":min = " << min << std::endl;
						std::cerr << " Material "<< nmat << ":max = " << max << std::endl;
						exit(EXIT_FAILURE);
					}
				}
			}
		}
	}
	
	return EXIT_SUCCESS;
}

} // end of namespace mp //
















/*
int initialise(std::string const infile){
	using namespace mp;
	
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(error_checking::error_check==true){std::cout << "initialise_variables has been called" << std::endl;}

	// Initialise system creation flags to zero
	for (int i=0;i<10;i++){
		system_creation_flags[i] = 0;
	}
	
	//-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
	//  Flag options
	//-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
	
	//-------------------------------------------------------------------
	// system_creation_flags[0] - Create system or read from file
	//-------------------------------------------------------------------
	//		0	Create system using code
	//	x	1	Read system generation information from source file (not yet implemented)
	//	x	2	Read system information from source file package (not yet implemented)
	//-------------------------------------------------------------------
	// system_creation_flags[1] - Set system particle shape
	//-------------------------------------------------------------------
	//		0	Full (use all atoms in lattice)
	//	x	1	Cube
	//		2	Cylinder
	//	x	3	Ellipsinder
	//		4	Sphere
	//		5	Truncated Octahedron
	//-------------------------------------------------------------------
	// system_creation_flags[2] - Set system type
	//-------------------------------------------------------------------
	//		0	Particle
	//		1	Particle Array
	//		2	Hexagonal Particle Array
	//		3	Voronoi Granular Film
	//		4	Grain Growth 2D Film
	//-------------------------------------------------------------------
	// system_creation_flags[3] - Set neighbourlist type
	//-------------------------------------------------------------------
	//		0	Explicit Jij(i,j,k)
	//	x	1	Range    Jij(r)
	//-------------------------------------------------------------------
	// system_creation_flags[4] - Set Multilayer Flag
	//-------------------------------------------------------------------
	//		0	Single Material
	//		1	Multilayer
	//		2	Random Intermixing
	system_creation_flags[0] = 0; // Create new system from scratch
	system_creation_flags[1] = 5; // Set system particle shape
	system_creation_flags[2] = 0; // Set system type
	system_creation_flags[3] = 0; // Set neighbourlist type
	system_creation_flags[4] = 0; // Set multilayer flag
	
	// output flags 
	// output flag increment mask
	
	//===================================================================
	//	Set hamiltonian Flags
	//===================================================================
	//
	//-------------------------------------------------------------------
	// hamiltonian_simulation_flags[0] - Calculate Exchange interaction
	//-------------------------------------------------------------------
	//		0	Disable Exchange Interaction
	//		1	Enable Exchange Interaction
	//-------------------------------------------------------------------
	// hamiltonian_simulation_flags[1] - Calculate Anisotropy
	//-------------------------------------------------------------------
	//		0	Disable Global Anisotropy
	//		1	Enable Global Uniaxial Anisotropy
	//		2	Enable Global Cubic Anisotropy
	//		3	Enable Local Anisotropy
	//-------------------------------------------------------------------
	// hamiltonian_simulation_flags[2] - Calculate Applied Fields
	//-------------------------------------------------------------------
	//		0	Disable Global Applied Field
	//		1	Enable Global Applied Field
	//		2	Enable Local Applied Field
	//-------------------------------------------------------------------
	// hamiltonian_simulation_flags[3] - Calculate Thermal Fields
	//-------------------------------------------------------------------
	//		0	Disable Thermal Fields
	//		1	Enable Thermal Fields
	//		2	Enable Local Thermal Fields
	//-------------------------------------------------------------------
	// hamiltonian_simulation_flags[4] - Calculate Dipolar Fields
	//-------------------------------------------------------------------
	//		0	Disable Dipolar Fields
	//		1	Enable Dipolar Fields
	//-------------------------------------------------------------------
	// hamiltonian_simulation_flags[5-9] - Calculate Additional Fields
	//-------------------------------------------------------------------
	//		0	Disable Additional Fields
	//		1	Enable Additional Field 1
	//		2	Enable Additional Field 2 ...

	sim::hamiltonian_simulation_flags[0] = 1;	// Exchange
	sim::hamiltonian_simulation_flags[1] = 1;	// Anisotropy
	sim::hamiltonian_simulation_flags[2] = 1;	// Applied
	sim::hamiltonian_simulation_flags[3] = 1;	// Thermal
	sim::hamiltonian_simulation_flags[4] = 0;	// Dipolar
	sim::hamiltonian_simulation_flags[5] = 0;	// Extra Term 2?
	sim::hamiltonian_simulation_flags[6] = 0;	// Extra Term 3?
	sim::hamiltonian_simulation_flags[7] = 0;	// Extra Term 4?
	sim::hamiltonian_simulation_flags[8] = 0;	// Extra Term 5?
	sim::hamiltonian_simulation_flags[9] = 0;	// Extra Term 6?
	
	sim::total_time = 1;			// total simulation time (single run)
	sim::loop_time = 0;			// time in loop, eg hysteresis, Tc
	sim::partial_time=1;			// time between statistics collection
	sim::equilibration_time=0;	// time for equilibration before main loop
	sim::temperature = 300.0;

	demag::demag_resolution=2;
	demag::update_rate=10000;
	//Integration parameters
	//alpha = 0.1;
	dt_SI = 1.0e-15;	// seconds
	dt = dt_SI*material_parameters::gamma_SI; // Must be set before Hth
	half_dt = 0.5*dt;


	//System Parameters
	//lattice_constant[0] = 2.5;
	//lattice_constant[1] = 2.5;
	//lattice_constant[2] = 2.5;

	lattice_constant[0] = 3.54;
	lattice_constant[1] = 3.54;
	lattice_constant[2] = 3.54;

	particle_creation_parity=1;

	
	crystal_structure = "fcc";
	hamiltonian_type = "generic";
	atomic_element[0]="Co ";
	atomic_element[1]="Gd ";
	
	system_dimensions[0] = 110.0;
	system_dimensions[1] = 110.0;
	system_dimensions[2] = 110.0;
	
	particle_scale   = 80.1;
	particle_spacing = 10.0;

	//----------------------------------------------------------
	// Multilayer Variables
	//----------------------------------------------------------
	multilayers::num_layers=1;
	multilayers::layer_dir=2; // x,y,z
	multilayers::interface_roughness=0.0;
	multilayers::spin_glass_density=0.0;
	multilayers::layer_limits_array[0]=0.0;
	multilayers::layer_limits_array[1]=0.2;
	multilayers::layer_limits_array[2]=0.6;
	multilayers::layer_limits_array[3]=0.8;
	//multilayers::layer_is_continuous[1]=true;
	//----------------------------------------------------------
	// Voronoi Variables
	//----------------------------------------------------------
	create_voronoi::voronoi_sd=0.0;
	create_voronoi::parity=0;
	//------------------------------------------------------------------------------
	// MPI Mode
	//------------------------------------------------------------------------------
	// 0 - Serial (Identical code on all nodes)
	// 1 - Parallel Stats (Same initialisation, different RNG seed, MPI statistics)
	// 2 - Parallel Decomposition (Different initialisation, code, full comms,etc)
	mpi_generic::mpi_mode=2;
	mpi_create_variables::mpi_interaction_range=4; // Unit cells
	mpi_create_variables::mpi_comms_identify=false;

	//------------------------------------------------------------------------------
	// Material Definitions
	//------------------------------------------------------------------------------
	num_materials=1;
	material.resize(num_materials);

	//-------------------------------------------------------
	// Material 0
	//-------------------------------------------------------
	material[0].name="CoPt";
	material[0].alpha=1.0;
	material[0].Jij_matrix_SI[0]=-3.0e-21;
	material[0].Jij_matrix_SI[1]=-1.4e-21;
	material[0].mu_s_SI=1.5*9.27400915e-24;
	material[0].Ku_SI=-3.58838e-23;
	material[0].gamma_rel=1.0;
	material[0].hamiltonian_type="generic";
	material[0].element="Ag ";
	//-------------------------------------------------------
	// Material 1
	//-------------------------------------------------------
	if(num_materials>1){
	material[1].name="Exchange Layer";
	material[1].alpha=1.0;
	material[1].Jij_matrix_SI[0]=-1.4e-21; // RE-TM
	material[1].Jij_matrix_SI[1]=-3.0e-21;  // RE-RE
	//material[1].Jij_matrix_SI[2]=1.6e-21;
	material[1].mu_s_SI=1.5*9.27400915e-24;
	material[1].Ku_SI=-0.807246e-25;
	material[1].gamma_rel=1.0;
	material[1].hamiltonian_type="generic";
	material[1].element="Li ";
	}
	//-------------------------------------------------------
	// Material 2
	//-------------------------------------------------------
	if(num_materials>2){
	material[2].name="Spin Glass";
	material[2].alpha=0.1;
	material[2].Jij_matrix_SI[0]=1.6e-21;
	material[2].Jij_matrix_SI[1]=1.6e-21;
	material[2].Jij_matrix_SI[2]=-1.6e-21;
	material[2].mu_s_SI=1.407e-23;
	material[2].Ku_SI=4.644e-23;
	material[2].gamma_rel=1.0;
	material[2].hamiltonian_type="generic";
	material[2].element="Fe ";
	}
	//print_mat();
	
	// Open main input file
	int iostat = vin::read(infile);
	if(iostat==EXIT_FAILURE){
		std::cerr << "Error - input file \'" << infile << "\' not found, exiting" << std::endl;
		exit(1);
	}
	
	
	lattice_space_conversion[0] = lattice_constant[0]*0.5;
	lattice_space_conversion[1] = lattice_constant[1]*0.5*0.333333333333333;
	lattice_space_conversion[2] = lattice_constant[2]*0.5;

	// Open material files
	
	
	
	
	//----------------------------------
	//Derived System Parameters
	//----------------------------------

	if(crystal_structure=="sc") num_nearest_neighbours = 6;
	if(crystal_structure=="fcc") num_nearest_neighbours = 12;
	if(crystal_structure=="rs") num_nearest_neighbours = 18;
	if(num_nearest_neighbours==0){
		 std::cout << "Error in determining num_nearest_neighbours for unknown crystal type\'" << crystal_structure << "\'" << std::endl;
		 exit(1);
	}
	
	if(hamiltonian_type=="generic")	hamiltonian_num_neighbours = num_nearest_neighbours;
	if(hamiltonian_num_neighbours==0){
		 std::cout << "Error in determining hamiltonian_num_neighbours - unknown Hamiltonian type!" << std::endl;
		 exit(1);
	}

	// Set derived material parameters
	for(int mat=0;mat<material_parameters::num_materials;mat++){
		material[mat].num_nearest_neighbours		= num_nearest_neighbours;
		material[mat].hamiltonian_num_neighbours	= hamiltonian_num_neighbours;
		material[mat].one_oneplusalpha_sq			=-material[mat].gamma_rel/(1.0+material[mat].alpha*material[mat].alpha);
		material[mat].alpha_oneplusalpha_sq			= material[mat].alpha*material[mat].one_oneplusalpha_sq;
		for(int j=0;j<num_materials;j++){
			material[mat].Jij_matrix[j]				= material[mat].Jij_matrix_SI[j]/material[mat].mu_s_SI;
		}
		material[mat].Ku									= material[mat].Ku_SI/material[mat].mu_s_SI;;
		material[mat].H_th_sigma						= sqrt(2.0*material[mat].alpha*1.3806503e-23/
																  (material[mat].mu_s_SI*material[mat].gamma_rel*dt));
	}

	//exit(0);

	int_system_dimensions[0] = 2*round(system_dimensions[0]/lattice_constant[0]);
	int_system_dimensions[1] = 6*round(system_dimensions[1]/lattice_constant[1]);
	int_system_dimensions[2] = 2*round(system_dimensions[2]/lattice_constant[2]);

	//----------------------------------
	// Enable/Disable Error Checking
	//----------------------------------

	error_checking::error_check=false;
	//error_checking::error_check=true;
	
	// Initialise random number generator
	
	mtrandom::grnd.seed(1234);
	//MTRand(1234);
	//MTRand grnd;
	//for(int i=0;i<10;i++){
	//	std::cout << drand() << std::endl;
	//}
	//exit (0);
return 0;
}
*/
