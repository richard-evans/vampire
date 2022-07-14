//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans and Sarah Jenkins 2021. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "montecarlo.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"


namespace program{

//--------------------------------------------------------------------------------
// constants but can be moved to input parameters if need be
//--------------------------------------------------------------------------------
const double exchange_stiffness_max_constraint_angle   = 180.01; // degrees
const double exchange_stiffness_delta_constraint_angle =  5; // 22.5 degrees
const double pi180 = M_PI/180.0;

//--------------------------------------------------------------------------------
// forward function declarations
//--------------------------------------------------------------------------------
void set_constraint_mask(const std::vector<double>& coordinates, // atomic coordinates along a principal direction x,y, or z
								 const std::vector<int>& material,       // material ID of each atom
								 std::vector<double>& fractional,        // fractional coordinate value (for initialisation)
								 std::vector<int>& mask,                 // mask showing which atoms are to be constrained
								 const int constrained);                 // material ID of atoms type to be constrained

void calculate_torque(const std::vector<int>& mask,
							 const std::vector<int>& material,
							 std::vector<double>& total_torques,
							 std::vector<double>& total_magnetizations);

//--------------------------------------------------------------------------------
// Program to calculate exchange stiffness for ferro, ferri and antiferromagnets
//--------------------------------------------------------------------------------
//
//      System is set up as follows:
//
//    |               |               :
//    |               |               :
//    |     free      |     free      :
//    |               |               :
//    |               |               :
//    ^               ^
//    constrained     constrained
//
//--------------------------------------------------------------------------------
//
//   The code assumes that no anisotropy is included, and performs a coherent
//   rotation of the central plane of atoms, while constraining the direction of
//   the first plane of the constrained material. The other materials are
//   integrated freely (standard Monte Carlo) allowing the determination of
//   exchange stiffness in ferri and collinear and non-collinear
//   antiferromagnets. The direction of the first plane of atoms is chosen
//   along z. The middle plane is the rotated in the y-z plane for small angles
//   (up to 30 degrees) to determine the torque curve <T>(theta). Linear
//   regression is then used to determine the gradient d<T>/dtheta and by
//   quadrature extract the change in free energy DF(theta).
//
//   Assuming the relationship
//
//         Eex = A <m1><m2> cos theta
//
//   we extract the effective exchange coupling for a given system and
//   temperature. Note, this program only works in serial due to the use of the
//   constrained monte carlo method.
//
//--------------------------------------------------------------------------------
void exchange_stiffness(){

 	// check calling of routine if error checking is activated
 	if(err::check==true) std::cout << "program::exchange_stiffness has been called" << std::endl;

	// Check integrator is hybrid CMC, if not then exit disgracefully
	if( sim::integrator != sim::hybrid_cmc ){
		err::zexit("Program exchange-anisotropy requires Hybrid Constrained Monte Carlo as the integrator. Exiting.");
	}

	//---------------------------------------------------------------------------
	// determine constrained material
	//---------------------------------------------------------------------------
	bool found = false;
	int constrained_material_id = 0;
	for(int m = 0; m < mp::num_materials; m++){
		if(mp::material[m].constrained){
			// check for previous constrained material
			if(found){
				std::cerr     << "Error: more than one constrained material defined in exchange stiffness program, where only one constrained material is allowed. Please define only a single constrained material." << std::endl;
				zlog << zTs() << "Error: more than one constrained material defined in exchange stiffness program, where only one constrained material is allowed. Please define only a single constrained material." << std::endl;
				err::vexit();
			}
			// ok this is the first constrained material
			found = true;
			// set constrained material
			constrained_material_id = m;
		}
	}

	// now check that at least one material is constrained
	if(!found){
		std::cerr     << "Error: no constrained material defined for exchange stiffness program, where one constrained material is required. Please define a single constrained material." << std::endl;
		zlog << zTs() << "Error: no constrained material defined for exchange stiffness program, where one constrained material is required. Please define a single constrained material." << std::endl;
		err::vexit();
	}

	//---------------------------------------------------------------------------
	// set up constrained atoms mask and vectors
	//---------------------------------------------------------------------------
	// mask for constraed atoms (0 = plane 1, 1 = plane 2, 2 = unconstrained)
	std::vector<int> constraint_mask(atoms::num_atoms);

	// fractional coordates of atom along constraint direction 0 ->1 ->0
	std::vector<double> fractional_coordinates(atoms::num_atoms);

	// x-direction
 	if( sim::domain_wall_axis == 0 ) set_constraint_mask(atoms::x_coord_array, atoms::type_array, fractional_coordinates, constraint_mask, constrained_material_id);
	if( sim::domain_wall_axis == 1 ) set_constraint_mask(atoms::y_coord_array, atoms::type_array, fractional_coordinates, constraint_mask, constrained_material_id);
	if( sim::domain_wall_axis == 2 ) set_constraint_mask(atoms::z_coord_array, atoms::type_array, fractional_coordinates, constraint_mask, constrained_material_id);

	// count number of atoms in each plane
	int n_atm_p1 = 0;
	int n_atm_p2 = 0;
	for( int atom = 0 ; atom < atoms::num_atoms; atom++){
		const int mask_id = constraint_mask[atom];
		// if atom is plane 1 or 2 then add them to count
		if(      mask_id == 0 ) n_atm_p1++;
		else if( mask_id == 1 ) n_atm_p2++;
	}
	// normalise and inverse for later use
	const double inv_n_atm_p1 = 1.0 / double( n_atm_p1 );
	const double inv_n_atm_p2 = 1.0 / double( n_atm_p2 );

	// Data structure to store torque curves data[torques][temperatures]
	std::vector< std::vector<double> > torque_data;
	std::vector< std::vector<double> > m1_data; // mean magnetization (layer 1)
	std::vector< std::vector<double> > m2_data; // mean magnetization (layer 2)
	std::vector< double > angles;
	std::vector< double > temperatures;
	//int count = 0;
	//for(double constraint_theta = 0.0; constraint_theta < mt; constraint_theta += dt) count++;
	//const int num_angles = count;
	//torque_data.resize(num_angle)

	// Open output file for torque data
	std::ofstream ofile;
	ofile.open("exchange-stiffness-torques.txt");

	//---------------------------------------------------------------------------
	// Main exchange calculation program
	//---------------------------------------------------------------------------
	const double mt = exchange_stiffness_max_constraint_angle;
	const double dt = exchange_stiffness_delta_constraint_angle;

	// set constraint phi component
	const double constraint_phi = 90.0;
	const double cosphi = cos(constraint_phi*pi180);
	const double sinphi = sin(constraint_phi*pi180);

	// loop over constraint angles ct (constraint_theta)
	for(double constraint_theta = 0.0; constraint_theta < mt; constraint_theta += dt){

		// push back data to store calculated torques and angles
		angles.push_back(constraint_theta*pi180); // radians
		m1_data.push_back( std::vector<double>() );
		m2_data.push_back( std::vector<double>() );
		torque_data.push_back( std::vector<double>() );

		// initialise new spin positions (half rotation from 1st plane to mid plane, then second half rotation to max plane)
		for( int atom = 0 ; atom < atoms::num_atoms; atom++ ){

			// get material ID
			const int mat = atoms::type_array[atom];

			// if constrained material initialise normal profile (constant angle)
			if( mat == constrained_material_id){
				const double theta = double(constraint_theta); // angle from z per fractional coordinate
				atoms::x_spin_array[atom] = cosphi*sin( theta * pi180 * fractional_coordinates[atom] );
				atoms::y_spin_array[atom] = sinphi*sin( theta * pi180 * fractional_coordinates[atom] );
				atoms::z_spin_array[atom] = cos( theta * pi180 * fractional_coordinates[atom] );
			}
			// otherwise assume initial spin direction for sublattice (should work OK for most ferro, ferri and antiferromagnets)
			else{
				atoms::x_spin_array[atom] = mp::material[mat].initial_spin[0];
				atoms::y_spin_array[atom] = mp::material[mat].initial_spin[1];
				atoms::z_spin_array[atom] = mp::material[mat].initial_spin[2];
			}

		}

		// reset hybrid CMC constraints
		std::vector<double> phi_theta_constraints(6, 0.0); // sets of (theta,phi)
		phi_theta_constraints[0] = 0.0; // set first plane along z
		phi_theta_constraints[1] = 0.0; // set first plane along z
		phi_theta_constraints[2] = constraint_theta; // set middle plane angle from z to theta
		phi_theta_constraints[3] = constraint_phi;   // set middle plane angle from z to theta
		std::vector<bool> constrained(3,false);
		constrained[0] = true; // constrain first plane
		constrained[1] = true; // constrain middle plane
		montecarlo::initialise_masked_cmc_mc(constraint_mask.size(), constraint_mask, constrained, phi_theta_constraints);

		// initialise temperature
		sim::temperature=sim::Tmin;

		// equilibrate ground state structure at zero kelvin
		for(int ztime = 0; ztime < 500; ztime++) sim::integrate(1);

		// reset temperature array
		temperatures.resize(0);

		// Perform Temperature Loop
		while( sim::temperature <= sim::Tmax){

			// reset torque and magnetization averages and counter
			std::vector<double> torques(6, 0.0);
			std::vector<double> magnetizations(2, 0.0);
			double counter = 0.0;

			// Equilibrate system
			sim::integrate(sim::equilibration_time);

			// Reset mean magnetisation counters
			stats::reset();

			// Reset start time
			int start_time=sim::time;

			// Simulate system
			while(sim::time<sim::loop_time+start_time){

				// Integrate system
				sim::integrate(sim::partial_time);

				// Calculate statistics
				stats::update();
				calculate_torque(constraint_mask, atoms::type_array, torques, magnetizations);
				counter += 1.0;

			}

			// Output data
			vout::data();

			// calculate and store mean torques
			ofile << sim::temperature << "\t" << double(constraint_theta) << "\t" << magnetizations[0]/counter << "\t" <<  magnetizations[1]/counter << "\t" <<
						torques[0]*inv_n_atm_p1/counter << "\t" << torques[1]*inv_n_atm_p1/counter << "\t" << torques[2]*inv_n_atm_p1/counter << "\t" <<
						torques[3]*inv_n_atm_p2/counter << "\t" << torques[4]*inv_n_atm_p2/counter << "\t" << torques[5]*inv_n_atm_p2/counter << std::endl;

			// store computed net torque in data array
			temperatures.push_back(sim::temperature);
			const int end = torque_data.size()-1;
			const double net_torque = torques[3]*inv_n_atm_p2/counter - torques[0]*inv_n_atm_p1/counter;
			m1_data[end].push_back(magnetizations[0]/counter);
			m2_data[end].push_back(magnetizations[1]/counter);
			torque_data[end].push_back(net_torque);

			// Increment temperature
			sim::temperature += sim::delta_temperature;

		} // end of temperature loop

	} // end of angle loop

	// close output file
	ofile.close();

	//--------------------------------------------------------------------------------------------------
	// perform linear regression to extract T(theta)
	//--------------------------------------------------------------------------------------------------

	// output data for checking
	/*for(int i=0; i<torque_data.size(); i++){
		std::cout << angles[i] << "\t";
		for(int j=0; j<torque_data[i].size(); j++){
			std::cout << torque_data[i][j] << "\t" << m1_data[i][j] << "\t" << m2_data[i][j] << "\t";
		}
		std::cout << std::endl;
	}*/

	// data structure to store torques in linear memory
	std::vector<double> torque1D( angles.size() );

	std::cout << "---------------------------------------------------------------" << std::endl;
	std::cout << " Final exchange fitting" << std::endl;
	std::cout << "---------------------------------------------------------------" << std::endl;

	// loop over all temperatures
	for(size_t j = 0; j < temperatures.size(); j++){

		// populate 1D data for fitting
		for(size_t i = 0; i < torque_data.size(); i++) torque1D[i] = torque_data[i][j];

		double m = 0.0;
		double c = 0.0;

		// compute gradient and intercept
		vmath::regression(angles, torque1D, m, c);

		// output resulting gradient to screen
		std::cout << temperatures[j] << "\t" << m << "\t" << c << std::endl;

	}



	return;

}

//-----------------------------------------------------------------------------------
// Function to determine planes of atoms to constrain
//-----------------------------------------------------------------------------------
void set_constraint_mask(const std::vector<double>& coordinates, // atomic coordinates along a principal direction x,y, or z
								 const std::vector<int>& material,       // material ID of each atom
								 std::vector<double>& fractional,        // fractional coordinate value (for initialisation)
								 std::vector<int>& mask,                 // mask showing which atoms are to be constrained
								 const int constrained){                 // material ID of atoms type to be constrained

	// set tolerance parameter for detemining unique plane of atoms
	const double tolerance = 0.1; // Angstroms

	// set constant number of atoms
	const int num_atoms = coordinates.size();

	//-----------------------------------------
	// determine first unique plane of atoms
	//-----------------------------------------

	// determine minimum and maximum system coordinates of constrained atoms
	double min_coord = 1.0e10; // minimum coordinates of constrained atom type
	double max_coord = -1.0e10; // maximum coordinates of constrained atom type
	for( int atom = 0 ; atom < num_atoms; atom++){
		// only consider positions of constrained atoms
		if(material[atom] == constrained){
			if(coordinates[atom] < min_coord) min_coord = coordinates[atom];
			if(coordinates[atom] > max_coord) max_coord = coordinates[atom];
		}
	}

	// determine coordinate of constrained atoms nearest middle
	double middle = min_coord + (max_coord - min_coord) * 0.5;
	double mid_coord = max_coord;
	double difference = 1.0e10;
	for( int atom = 0 ; atom < num_atoms; atom++){
		// only consider positions of constrained atoms
		if(material[atom] == constrained){
			// if atom is closer than current closest
			if(fabs(coordinates[atom] - middle) < difference){
				// set atom as new middle coordinate
				mid_coord = coordinates[atom];
				// set up new difference, so only atoms closer than this one are considered
				difference = fabs(coordinates[atom] - middle);
			}
		}
	}

	// make sure that there is a reasonable separation between planes (at least 1 Angstrom)
	const double separation_A = mid_coord - min_coord;
	const double separation_B = max_coord - mid_coord;

	if( fabs(separation_A) < 1.0 || fabs(separation_B) < 1.0 ){
		std::cerr     << "Error: insufficent separation between atoms along domain wall direction to enable sensible calculation of exchange stiffness. Please increase system size along domain wall direction." << std::endl;
		zlog << zTs() << "Error: insufficent separation between atoms along domain wall direction to enable sensible calculation of exchange stiffness. Please increase system size along domain wall direction." << std::endl;
		err::vexit();
	}

	//---------------------------------------------------------------------------
	// Determine constrained mask and number of atoms in each plane
	//---------------------------------------------------------------------------
	int num_atoms_in_plane_1 = 0;
	int num_atoms_in_plane_2 = 0;

	// calculate fractional extent of atoms
	// (this is not absolute, so there may be some atoms with fractional coordinates < 0 and >1)
	const double inv_delta = 2.0 / (max_coord - min_coord);

	for( int atom = 0 ; atom < num_atoms; atom++){

		const double c = coordinates[atom];

		if(material[atom] == constrained){
			// check for first plane
			if( c > min_coord - tolerance && c < min_coord + tolerance ){
				mask[atom] = 0; // set mask
				num_atoms_in_plane_1++; // add counter
			}
			// check for second plane
			else if( c > mid_coord - tolerance && c < mid_coord + tolerance ){
				mask[atom] = 1; // set mask
				num_atoms_in_plane_2++; // add counter
			}
			else{
				mask[atom] = 2; // set mask (unconstrained)
			}
		}
		else{
			// always set unconstrained materials to unconstrained
			mask[atom] = 2; // set mask
		}
		// calculate fractional coordinates along domain wall direction (reversing for c > 0.5)
		fractional[atom] = c * inv_delta <= 1.0 ? c * inv_delta : 2.0 - c * inv_delta;
		//std::cout << atom << "\t" << c << "\t" << c*inv_delta << "\t" << fractional[atom] << std::endl;
	}

	//---------------------------------------------------------------------------
	// print information to log file for debugging
	//---------------------------------------------------------------------------
	zlog << zTs() << "Statistics for exchange-stiffnes initialisation:" << std::endl;
	zlog << zTs() << "   number of atoms in plane 1: " << num_atoms_in_plane_1 << std::endl;
	zlog << zTs() << "   number of atoms in plane 2: " << num_atoms_in_plane_2 << std::endl;
	zlog << zTs() << "   coordinates of plane 1:     " << min_coord << std::endl;
	zlog << zTs() << "   coordinates of plane 2:     " << mid_coord << std::endl;
	zlog << zTs() << "   coordinates of plane 3:     " << max_coord << std::endl;
	zlog << zTs() << "   separation plane 1-2:       " << separation_A << std::endl;
	zlog << zTs() << "   separation plane 2-3:       " << separation_B << std::endl;

	return;

}
//------------------------------------------------------------------------------
// Function to calculate equal and opposite torques on each constrained plane
//------------------------------------------------------------------------------
void calculate_torque(const std::vector<int>& mask,
							 const std::vector<int>& material,
							 std::vector<double>& total_torques,
							 std::vector<double>& total_magnetizations){

	// calculate net fields on all spins
	sim::calculate_spin_fields(0,atoms::num_atoms);
	sim::calculate_external_fields(0,atoms::num_atoms);

	double mm[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double tt[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double counter[2] = { 0.0, 0.0 };

	for(size_t atom = 0; atom < mask.size(); atom++){

		const int mask_id = mask[atom];

		// if atom is plane 1 or 2 then compute total torque
		if( mask_id == 0 || mask_id == 1 ){

			// get atomic moment
			const int mat = material[atom];
			const double mu = mp::material[mat].mu_s_SI;

			// Store local spin in Sand local field in H
			const double S[3] = { atoms::x_spin_array[atom],
										 atoms::y_spin_array[atom],
										 atoms::z_spin_array[atom]};

			const double H[3] = { atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
										 atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
										 atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] };

			// compute torques
			const double ttx = mu*(S[1]*H[2] - S[2]*H[1]);
			const double tty = mu*(S[2]*H[0] - S[0]*H[2]);
			const double ttz = mu*(S[0]*H[1] - S[1]*H[0]);

			// accumulate torques for local sum
			tt[3*mask_id + 0] += ttx;
			tt[3*mask_id + 1] += tty;
			tt[3*mask_id + 2] += ttz;

			// accumulate spin moments for local sum
			mm[3*mask_id+0] += S[0];
			mm[3*mask_id+1] += S[1];
			mm[3*mask_id+2] += S[2];

			// update counter for each mask
			counter[mask_id] += 1.0;

		}

	}

	// copy local torques to main total torque array
	for(int i=0; i<6; i++) total_torques[i] += tt[i];

	// calculate magnetizations for planes 1 and 2
	total_magnetizations[0] += sqrt(mm[0]*mm[0] + mm[1]*mm[1] + mm[2]*mm[2]) / counter[0];
	total_magnetizations[1] += sqrt(mm[3]*mm[3] + mm[4]*mm[4] + mm[5]*mm[5]) / counter[1];

	return;

}

} //end of namespace program
