/// Program to convert vampire cfg files to povray format
///
/// ./cfg2rasmol 

// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct material_t{
	
	double mx,my,mz,mu_s,magm;

};

void rgb( double ireal, double &red, double &green, double &blue){

			if(ireal>0.8){
				red = 0.0;
				green = 0.0;
				blue = 1.0;
			}
			else if(ireal>=0.0){
				red = 1.0-ireal*1.2;
				green = 1.0-ireal*1.2;
				blue = 1.0;
			}
			else if(ireal>=-0.8){
				red = 1.0;
				green = 1.0+ireal*1.2;
				blue = 1.0+ireal*1.2;
			}
			else if(ireal<-0.8){
				red = 1.0;
				green = 0.0;
				blue = 0.0;
			}
			else{
				red = 1.0;
				green = 1.0;
				blue = 1.0;
			}

			if(blue<0.0) blue=0.0;
			if(red<0.0) red=0.0;
			if(green<0.0) green=0.0;
			
}

int main(){
	
	std::vector <int> mat(0);
	std::vector <int> cat(0); 
	std::vector <double> cx(0);
	std::vector <double> cy(0);
	std::vector <double> cz(0);
	std::vector <std::string> type(0);
	std::vector <std::string> filenames(0);
	
	// open coordinate file
	std::ifstream coord_file;
	coord_file.open("atoms-coords.cfg");

	// read in file header
	std::string dummy;
	
	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	// get number of atoms
	unsigned int n_atoms;
	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;	
	dummy.erase (dummy.begin(), dummy.begin()+17);
	n_atoms=atoi(dummy.c_str());

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	// get number of subsidiary files
	unsigned int n_files;
	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;	
	dummy.erase (dummy.begin(), dummy.begin()+22);
	n_files=atoi(dummy.c_str());

	for(int file=0; file<n_files; file++){
		getline(coord_file,dummy);
		filenames.push_back(dummy);
		//std::cout << filenames[file] << std::endl;
	}

	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;

	unsigned int n_local_atoms;
	getline(coord_file,dummy);
	//std::cout << dummy << std::endl;	
	n_local_atoms=atoi(dummy.c_str());
	
	// resize arrays
	mat.resize(n_atoms);
	cat.resize(n_atoms); 
	cx.resize(n_atoms);
	cy.resize(n_atoms);
	cz.resize(n_atoms);
	type.resize(n_atoms);
	
	unsigned int counter=0;
	
	// finish reading master file coordinates
	for(int i=0;i<n_local_atoms;i++){
		coord_file >> mat[counter] >> cat[counter] >> cx[counter] >> cy[counter] >> cz[counter] >> type[counter];
		counter++;
	}
	
	// close master file
	coord_file.close();
	
	// now read subsidiary files
	for(int file=0; file<n_files; file++){
		std::ifstream infile;
		infile.open(filenames[file].c_str());
		
		// read number of atoms in this file
		getline(infile,dummy);
		n_local_atoms=atoi(dummy.c_str());
		for(int i=0;i<n_local_atoms;i++){
			infile >> mat[counter] >> cat[counter] >> cx[counter] >> cy[counter] >> cz[counter] >> type[counter];
			counter++;
		}
		// close subsidiary file
		infile.close();
	}
	
	// check for correct read in of coordinates
	if(counter!=n_atoms) std::cerr << "Error in reading in coordinates" << std::endl;

	// Loop over all possible spin config files
	int ios = 0;
	int spinfile_counter=0;
	
	while(ios==0){
		
		std::stringstream file_sstr;
		file_sstr << "atoms-";
		file_sstr << std::setfill('0') << std::setw(8) << spinfile_counter;
		file_sstr << ".cfg";
		std::string cfg_file = file_sstr.str();
		//const char* cfg_filec = cfg_file.c_str();

		std::ifstream spinfile;
		spinfile.open(cfg_file.c_str());

		if(spinfile.is_open()){
		
			std::cout << "Processing file: " << cfg_file << std::endl;
			
		// Read in file header
		getline(spinfile,dummy);
			
	//std::cout << dummy << std::endl;

	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;

	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;

	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;

	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;

	// get number of atoms
	unsigned int n_spins;
	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;	
	dummy.erase (dummy.begin(), dummy.begin()+17);
	n_spins=atoi(dummy.c_str());
	if(n_spins!=n_atoms) std::cerr << "Error! - mismatch between number of atoms in coordinate and spin files" << std::endl;
	
	getline(spinfile,dummy); // sys dimensions
	//std::cout << dummy << std::endl;
	char const field_delim = '\t';
	dummy.erase (dummy.begin(), dummy.begin()+18);
	std::istringstream ss(dummy);
	std::vector<double> val(0);
	for (std::string num; getline(ss, num, field_delim); ) {
		val.push_back(atof(num.c_str()));
		//std::cout << num << "\t" << val << std::endl;
	}
	double dim[3] = {val[0],val[1],val[2]};
	//std::cout << dim[0] << "\t" << dim[1] << "\t" << dim[2] << std::endl;

	getline(spinfile,dummy); // coord file

	getline(spinfile,dummy); // time
	//std::cout << dummy << std::endl;	
	dummy.erase (dummy.begin(), dummy.begin()+5);
	double time = atof(dummy.c_str());

	getline(spinfile,dummy); // field
	//std::cout << dummy << std::endl;	
	dummy.erase (dummy.begin(), dummy.begin()+6);
	double field = atof(dummy.c_str());

	getline(spinfile,dummy); // temp
	//std::cout << dummy << std::endl;	
	dummy.erase (dummy.begin(), dummy.begin()+12);
	double temperature = atof(dummy.c_str());
	
	// magnetisation
	getline(spinfile,dummy);
	dummy.erase (dummy.begin(), dummy.begin()+14);
	std::istringstream ss2(dummy);
	val.resize(0);
	for (std::string num; getline(ss2, num, field_delim); ) {
		val.push_back(atof(num.c_str()));
		//std::cout << num << "\t" << val << std::endl;
	}
	double mx = val[0];
	double my = val[1];
	double mz = val[2];
	
	// get number of materials
	unsigned int n_mat;
	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;	
	dummy.erase (dummy.begin(), dummy.begin()+20);
	n_mat=atoi(dummy.c_str());

	std::vector <material_t> material(n_mat); 
	for(int imat=0;imat<n_mat;imat++){
		getline(spinfile,dummy);
		//dummy.erase (dummy.begin(), dummy.begin()+22);
		std::istringstream ss(dummy);
		std::vector<double> val(0);
		for (std::string num; getline(ss, num, field_delim); ) {
			val.push_back(atof(num.c_str()));
		}
		material[imat].mu_s = val[0];
		material[imat].mx = val[1];
		material[imat].my = val[2];
		material[imat].mz = val[3];
		material[imat].magm = val[4];
		//std::cout << material[imat].mu_s << "\t" << material[imat].mx << "\t" << material[imat].my << "\t"  << material[imat].mz << "\t"  << material[imat].magm << std::endl;
	}
	
	// line
	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;
	
	// get number of subsidiary files
	unsigned int n_files;
	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;	
	dummy.erase (dummy.begin(), dummy.begin()+22);
	n_files=atoi(dummy.c_str());
	filenames.resize(0);
	
	for(int file=0; file<n_files; file++){
		getline(spinfile,dummy);
		filenames.push_back(dummy);
		//std::cout << filenames[file] << std::endl;
	}

	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;

	unsigned int n_local_atoms;
	getline(spinfile,dummy);
	//std::cout << dummy << std::endl;	
	n_local_atoms=atoi(dummy.c_str());

	// Open Povray Include File
	std::stringstream incpov_file_sstr;
	incpov_file_sstr << "atoms-";
	incpov_file_sstr << std::setfill('0') << std::setw(8) << spinfile_counter;
	incpov_file_sstr << ".inc";
	std::string incpov_file = incpov_file_sstr.str();
	
	// Open Povray Output file
	std::stringstream pov_file_sstr;
	pov_file_sstr << "atoms-";
	pov_file_sstr << std::setfill('0') << std::setw(8) << spinfile_counter;
	pov_file_sstr << ".pov";
	std::string pov_file = pov_file_sstr.str();

	std::ofstream pfile;
	pfile.open(pov_file.c_str());
	
	// Ouput povray file header
			double size, mag_vec;
			double vec[3];

			size = sqrt(dim[0]*dim[0] + dim[1]*dim[1] + dim[2]*dim[2]);
			vec[0] = (1.0/dim[0]);
			vec[1] = (1.0/dim[1]);
			vec[2] = (1.0/dim[2]);
			mag_vec = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
			vec[0]/=mag_vec;
			vec[1]/=mag_vec;
			vec[2]/=mag_vec;

			pfile << "#include \"colors.inc\"" << std::endl;
			pfile << "#include \"metals.inc\""	<< std::endl;
			pfile << "#include \"screen.inc\""	<< std::endl;
			pfile << "#declare LX=" << dim[0]*0.5 << ";" << std::endl;
			pfile << "#declare LY=" << dim[1]*0.5 << ";" << std::endl;
			pfile << "#declare LZ=" << dim[2]*0.5 << ";" << std::endl;
			pfile << "#declare CX=" << size*vec[0]*6.0 << ";" << std::endl;
			pfile << "#declare CY=" << size*vec[1]*6.0 << ";" << std::endl;
			pfile << "#declare CZ=" << size*vec[2]*6.0 << ";" << std::endl;
	 		pfile << "#declare ref=0.05;" << std::endl;
			pfile << "global_settings { assumed_gamma 2.0 }" << std::endl;
			pfile << "background { color Gray30 }" << std::endl;

			pfile << "Set_Camera(<CX,CY,CZ>, <LX,LY,LZ>, 15)" << std::endl;
			pfile << "Set_Camera_Aspect(4,3)" << std::endl;
			pfile << "Set_Camera_Sky(<0,0,1>)" << std::endl;
			pfile << "light_source { <2*CX, 2*CY, 2*CZ> color White}" << std::endl;

			for(int imat=0;imat<n_mat;imat++){
				pfile << "#declare sscale"<< imat << "=2.0;" << std::endl;
				pfile << "#declare rscale"<< imat << "=1.2;" << std::endl;
				pfile << "#declare cscale"<< imat << "=3.54;" << std::endl;
				pfile << "#declare cones"<< imat << "=0;" << std::endl;
				pfile << "#declare arrows"<< imat << "=1;" << std::endl;
				pfile << "#declare spheres"<< imat << "=1;" << std::endl;
				pfile << "#declare cubes" << imat << "=0;" << std::endl;
				pfile << "#declare spincolors"<< imat << "=1;" << std::endl;
				pfile << "#declare spincolor"<< imat << "=pigment {color rgb < 0.1 0.1 0.1 >};" << std::endl;
				pfile << "#macro spinm"<< imat << "(cx,cy,cz,sx,sy,sz, cr,cg,cb)" << std::endl;
				pfile << "union{" << std::endl;
				pfile << "#if(spheres" << imat << ") sphere {<cx,cy,cz>,0.5*rscale"<< imat << "} #end" << std::endl;
				pfile << "#if(cubes" << imat << ") box {<cx-cscale"<< imat << "*0.5,cy-cscale" << imat << "*0.5,cz-cscale"<< imat << "*0.5>,<cx+cscale"<< imat << "*0.5,cy+cscale" << imat << "*0.5,cz+cscale"<< imat << "*0.5>} #end" << std::endl;
				pfile << "#if(cones"<< imat << ") cone {<cx+0.5*sx*sscale0,cy+0.5*sy*sscale"<< imat << ",cz+0.5*sz*sscale"<< imat << ">,0.0 <cx-0.5*sx*sscale"<< imat << ",cy-0.5*sy*sscale"<< imat << ",cz-0.5*sz*sscale"<< imat << ">,sscale0*0.5} #end" << std::endl;
				pfile << "#if(arrows" << imat << ") cylinder {<cx+sx*0.5*sscale"<< imat <<",cy+sy*0.5*sscale"<< imat <<",cz+sz*0.5*sscale"<< imat <<
							">,<cx-sx*0.5*sscale"<< imat <<",cy-sy*0.5*sscale"<< imat <<",cz-sz*0.5*sscale"<< imat <<">,sscale"<< imat <<"*0.12}";
				pfile << "cone {<cx+sx*0.5*1.6*sscale"<< imat <<",cy+sy*0.5*1.6*sscale"<< imat <<",cz+sz*0.5*1.6*sscale"<< imat <<">,sscale"<< imat <<"*0.0 <cx+sx*0.5*sscale"<< imat <<
							",cy+sy*0.5*sscale"<< imat <<",cz+sz*0.5*sscale"<< imat <<">,sscale"<< imat <<"*0.2} #end" << std::endl;
				pfile << "#if(spincolors"<< imat << ") texture { pigment {color rgb <cr cg cb>}finish {reflection {ref} diffuse 1 ambient 0}}" << std::endl;
				pfile << "#else texture { spincolor"<< imat << " finish {reflection {ref} diffuse 1 ambient 0}} #end" << std::endl;
				pfile << "}" << std::endl;
				pfile << "#end" << std::endl;
				
				
			}
			pfile << "#include \"" << incpov_file_sstr.str() << "\"" << std::endl;
	
		pfile.close();


	std::ofstream incpfile;
	incpfile.open(incpov_file.c_str());
	
	double sx,sy,sz,red,green,blue,ireal;
	unsigned int si=0;
	
	// Read in spin coordinates and output to povray file
	for(int i=0; i<n_local_atoms;i++){
		spinfile >> sx >> sy >> sz;
		rgb(sz,red,green,blue);
		incpfile << "spinm"<< mat[si] << "(" << cx[si] << "," << cy[si] << "," << cz[si] << "," 
		<< sx << "," << sy << "," << sz << "," << red << "," << green << "," << blue << ")" << std::endl;
		si++;
		
	}
	
		// close master file
	spinfile.close();
	
	// now read subsidiary files
	for(int file=0; file<n_files; file++){
		std::ifstream infile;
		infile.open(filenames[file].c_str());
		
		// read number of atoms in this file
		getline(infile,dummy);
		n_local_atoms=atoi(dummy.c_str());
		for(int i=0;i<n_local_atoms;i++){
			infile >> sx >> sy >> sz;
			rgb(sz,red,green,blue);
			incpfile << "spinm"<< mat[si] << "(" << cx[si] << "," << cy[si] << "," << cz[si] << "," 
			<< sx << "," << sy << "," << sz << "," << red << "," << green << "," << blue << ")" << std::endl;
			si++;
			
		}
		// close subsidiary file
		infile.close();
	}
	
	// close povray inc file
	incpfile.close();
	
	
		spinfile_counter++;
		}
		else ios=1;
	}


	
	// finished
	return 0;

}
