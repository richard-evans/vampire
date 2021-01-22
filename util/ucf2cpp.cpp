//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

namespace unitcell{

   //---------------------------------------------------------------------------
   // Unit cell atom class definition
   //---------------------------------------------------------------------------
   class atom_t {
	public:
      double x; /// atom x-coordinate
      double y; /// atom y-coordinate
      double z; /// atom z-coordinate
      unsigned int mat; /// material
      unsigned int lc; /// lattice category
      unsigned int hc; /// height category
      unsigned int ni; /// number of interactions
	};

   //---------------------------------------------------------------------------
   // Unit cell interaction class definition
   //---------------------------------------------------------------------------
	class interaction_t {
	public:
      unsigned int i; /// atom unit cell id
      unsigned int j; /// neighbour atom unit cell id
      int dx; /// delta x in unit cells
      int dy; /// delta y in unit cells
      int dz; /// delta z in unit cells
      double rij; // interaction range (unit cells)
      double Jij[3][3]; /// Exchange tensor
	};

   //---------------------------------------------------------------------------
   // Unit cell class definition
   //---------------------------------------------------------------------------
	class unit_cell_t {
	public:

		double dimensions[3];
		double shape[3][3];
      double cutoff_radius; // nearest neighbours
      unsigned int interaction_range; /// maximum range in unit cells

		unsigned int lcsize; /// number of local categories
		unsigned int hcsize; /// number of height categories
		unsigned int surface_threshold; /// threshold for surface atoms

		// list of atoms in each unit cell
		std::vector <unitcell::atom_t> atom;

      //unitcell::exchange_template_t bilinear;
      //unitcell::exchange_template_t biquadratic;
      //exchange_template_t fourspin_interaction; // tbc

	} unit_cell; // evil global variable
}

//-----------------------------------------------------------------------------------------
// function to translate atomic position but wrap-around to remain in the range 0-0.999999
//-----------------------------------------------------------------------------------------
double translate(double value, double shift){

   value = value+shift;
   if(value < 0.0) value += 1.0;
   else if(value > 0.999999) value = value - 1.0;

   // check values are sensible
   if(value < 0.0 || value > 1.0){
      std::cerr << "Error! shift is illegal - initial value must be in the range 0-1" << std::endl;
      exit(1);
   }

   return value;

}

int translate(int value, int shift, int max){

   value = value+shift;
   if(value < 0) value = max;
   else if(value > max) value = 1;

   // check values are sensible
   if(value < 0 || value > max){
      std::cerr << "Error! int shift is illegal - initial value must be in the range 0-1" << std::endl;
      exit(1);
   }

   return value;

}

namespace vin {

   //---------------------------------------------------------------------------
   // Function to open file on master process and return string on all
   // processes containing file contents
   //---------------------------------------------------------------------------
   std::string get_string(std::string const filename, std::string source_file_name, int line){

      // boolean variable specifying root process
      bool root = true;

      // number of characters in file (needed by all processors)
      uint64_t len = 0;

      // message buffer to store processed string as characters suitable for MPI Broadcast
      std::vector<char> message(0);

      // Read in file on root
      if (root){

         // ifstream declaration
         std::ifstream inputfile;

         // Open file
         inputfile.open(filename.c_str());

         // Check for correct opening
         if(!inputfile.is_open()){
            if(line >= 0) std::cerr << "Error opening input file \"" << filename << "\" specified on line " << line << " of " << source_file_name << " : File does not exist or cannot be opened! Exiting!" << std::endl;
            else          std::cerr << "Error opening input file \"" << filename << "\" : File does not exist or cannot be opened! Exiting!" << std::endl;
            exit(1);
         }

         // load file directly into std::string
         std::string contents( (std::istreambuf_iterator<char>(inputfile)),
                                std::istreambuf_iterator<char>());

         // get total number of characters in file
         len = contents.length();

         // reserve correct amount of storage for message
         message.reserve(len);

         // copy contents to message buffer for broadcast
         std::copy(contents.begin(), contents.end(), std::back_inserter(message));

      }

      // return message array cast to a std::string on all processors
      std::string result_str(message.begin(),message.end());
      return result_str;

   }

} // end of namespace vin

//------------------------------------------------------------------------------
// Simple code to convert vampire ucf file to code to generate nearest
// neighbour model in VAMPIRE
//------------------------------------------------------------------------------
int main(int argc, char* argv[]){

   using namespace unitcell;

   std::string filename;

   if(argc > 1){
      filename = std::string(argv[1]);
   }
   else{
      std::cerr << "Error - no file name specified" << std::endl;
      exit(EXIT_FAILURE);
   }

   // stringstream stream declaration
	std::stringstream inputfile;

   // fill input file stream with contents of file opened on master process
   inputfile.str( vin::get_string(filename.c_str(), "input", -1) );

   std::cout << "done!\nProcessing unit cell data..." << std::flush;

	// keep record of current line
	unsigned int line_counter=0;
	unsigned int line_id=0;

   std::string exchange_type_string; // string defining exchange type

   // defaults for interaction list
   int interaction_range = 1; // assume +-1 unit cell as default

	// Loop over all lines
	while (! inputfile.eof() ){
		line_counter++;
		// read in whole line
		std::string line;
		getline(inputfile,line);
		//std::cout << line.c_str() << std::endl;

		// ignore blank lines
		std::string empty="";
		if(line==empty) continue;

		// set character triggers
		const char* hash="#";	// Comment identifier

		bool has_hash=false;
		// Determine if line is a comment line
		for(unsigned int i=0;i<line.length();i++){
			char c=line.at(i);

			if(c== *hash){
					has_hash=true;
					break;
			}
		}
		// if hash character found then read next line
		if(has_hash==true) continue;

		// convert line to string stream
		std::istringstream iss(line,std::istringstream::in);

		// non-comment line found - check for line number
		switch(line_id){
			case 0:
				iss >> unit_cell.dimensions[0] >> unit_cell.dimensions[1] >> unit_cell.dimensions[2];
				break;
			case 1:
				iss >> unit_cell.shape[0][0] >> unit_cell.shape[0][1] >> unit_cell.shape[0][2];
				break;
			case 2:
				iss >> unit_cell.shape[1][0] >> unit_cell.shape[1][1] >> unit_cell.shape[1][2];
				break;
			case 3:
				iss >> unit_cell.shape[2][0] >> unit_cell.shape[2][1] >> unit_cell.shape[2][2];
				break;
			case 4:
				int num_uc_atoms;
				iss >> num_uc_atoms;
				//std::cout << "Reading in " << num_uc_atoms << " atoms" << std::endl;
				// resize unit_cell.atom array if within allowable bounds
				if( (num_uc_atoms >0) && (num_uc_atoms <= 1000000)) unit_cell.atom.resize(num_uc_atoms);
				else {
					std::cerr << "Error! Requested number of atoms " << num_uc_atoms << " on line " << line_counter
					<< " of unit cell input file " << filename.c_str() << " is outside of valid range 1-1,000,000. Exiting" << std::endl;
				}

            std::cout << "\nProcessing data for " << unit_cell.atom.size() << " atoms..." << std::flush;


            // loop over all atoms and read into class
            for(unsigned int i = 0; i < unit_cell.atom.size(); i++){

					line_counter++;

					// declare safe temporaries for atom input
					int id=i;
					double cx=2.0, cy=2.0,cz=2.0; // coordinates - default will give an error
					int mat_id=0, lcat_id=0, hcat_id=0; // sensible defaults if omitted
					// get line
					std::string atom_line;
					getline(inputfile,atom_line);
					std::istringstream atom_iss(atom_line,std::istringstream::in);
					atom_iss >> id >> cx >> cy >> cz >> mat_id >> lcat_id >> hcat_id;
					//std::cout << id << "\t" << cx << "\t" << cy << "\t" << cz<< "\t"  << mat_id << "\t" << lcat_id << "\t" << hcat_id << std::endl;
					//inputfile >> id >> cx >> cy >> cz >> mat_id >> lcat_id >> hcat_id;
					// now check for mostly sane input
					if(cx>=0.0 && cx <=1.0) unit_cell.atom[i].x=cx;
					else{
						std::cerr << "Error! atom x-coordinate for atom " << id << " on line " << line_counter
									 << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
					}
					if(cy>=0.0 && cy <=1.0) unit_cell.atom[i].y=cy;
					else{
						std::cerr << "Error! atom y-coordinate for atom " << id << " on line " << line_counter
									 << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
					}
					if(cz>=0.0 && cz <=1.0) unit_cell.atom[i].z=cz;
					else{
						std::cerr << "Error! atom z-coordinate for atom " << id << " on line " << line_counter
						<< " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
					}
					if(mat_id >=0) unit_cell.atom[i].mat=mat_id;
					else{
						std::cerr << "Error! Requested material id " << mat_id << " for atom number " << id <<  " on line " << line_counter
									 << " of unit cell input file " << filename.c_str() << " is greater than the number of materials ( " << 100 << " ) specified in the material file. Exiting" << std::endl;
               }
					unit_cell.atom[i].lc=lcat_id;
					unit_cell.atom[i].hc=hcat_id;
					//std::cout << i << "\t" << id << "\t" << cx << "\t" << cy << "\t" << cz << "\t" << mat_id << "\t" << lcat_id << "\t" << hcat_id << std::endl;
				}
				break;
			/*case 5:{

            // read (bilinear) exchange interactions
            unit_cell.bilinear.read_interactions(num_uc_atoms, inputfile, iss, filename, line_counter, interaction_range);
				break;

         }
         case 6:{

            // read biquadratic exchange interactions
            unit_cell.biquadratic.read_interactions(num_uc_atoms, inputfile, iss, filename, line_counter, interaction_range);
            break;

         }*/

			default:
            // ignore lines
				//terminaltextcolor(RED);
				//std::cerr << "Error! Unknown line type on line " << line_counter
				//	<< " of unit cell input file " << filename.c_str() << ". Exiting" << std::endl; err::vexit();
				//terminaltextcolor(WHITE);
            break;
		}
		line_id++;
	} // end of while loop

   //---------------------------------------------------------------------------
   // output data to screen as code
   //---------------------------------------------------------------------------
   std::cout << "\t// Set basic unit cell properties\n" << std::endl;
   // always use normalised cell sizes"
   std::cout << "\tunit_cell.dimensions[0] = 1.0;" << std::endl;
   std::cout << "\tunit_cell.dimensions[1] = 1.0;" << std::endl;
   std::cout << "\tunit_cell.dimensions[2] = 1.0;\n" << std::endl;

   std::cout << "\tunit_cell.shape[0][0] = " << unit_cell.shape[0][0] << ";" << std::endl;
   std::cout << "\tunit_cell.shape[0][1] = " << unit_cell.shape[0][1] << ";" << std::endl;
   std::cout << "\tunit_cell.shape[0][2] = " << unit_cell.shape[0][2] << ";\n" << std::endl;

   std::cout << "\tunit_cell.shape[1][0] = " << unit_cell.shape[1][0] << ";" << std::endl;
   std::cout << "\tunit_cell.shape[1][1] = " << unit_cell.shape[1][1] << ";" << std::endl;
   std::cout << "\tunit_cell.shape[1][2] = " << unit_cell.shape[1][2] << ";\n" << std::endl;

   std::cout << "\tunit_cell.shape[2][0] = " << unit_cell.shape[2][0] << ";" << std::endl;
   std::cout << "\tunit_cell.shape[2][1] = " << unit_cell.shape[2][1] << ";" << std::endl;
   std::cout << "\tunit_cell.shape[2][2] = " << unit_cell.shape[2][2] << ";\n" << std::endl;

   // determine max lc and hc
   int max_lc = 0;
   int max_hc = 0;
   int max_mat = 0;
   for(unsigned int i = 0; i < unit_cell.atom.size(); i++){
      if(unit_cell.atom[i].mat > max_mat) max_mat = unit_cell.atom[0].mat;
      if(unit_cell.atom[i].lc > max_lc) max_lc = unit_cell.atom[i].lc;
      if(unit_cell.atom[i].hc > max_hc) max_hc = unit_cell.atom[i].hc;
   }

   std::cout << "\tunit_cell.lcsize = " << max_lc+1 << ";" << std::endl;
   std::cout << "\tunit_cell.hcsize = " << max_hc+1 << ";" << std::endl;

   // always assume 1 unit cell interactions range
   std::cout << "\tunit_cell.interaction_range = 1;" << std::endl;
   std::cout << "\tunit_cell.atom.resize(" << unit_cell.atom.size() << ");" << std::endl;

   // surface threshold is not walways simple...
   std::cout << "\tunit_cell.surface_threshold = 8;\n" << std::endl;

   for(unsigned int i = 0; i < unit_cell.atom.size(); i++){
      std::cout << "\t//-----------------------------" << std::endl;
      std::cout << "\tunit_cell.atom[" << i << "].x   = " << translate(unit_cell.atom[i].x,0.0) << ";" << std::endl;
      std::cout << "\tunit_cell.atom[" << i << "].y   = " << translate(unit_cell.atom[i].y,0.0) << ";" << std::endl;
      std::cout << "\tunit_cell.atom[" << i << "].z   = " << translate(unit_cell.atom[i].z,0.0) << ";" << std::endl;
      std::cout << "\tunit_cell.atom[" << i << "].mat = " << unit_cell.atom[i].mat << ";" << std::endl;
      std::cout << "\tunit_cell.atom[" << i << "].lc  = " << unit_cell.atom[i].lc << ";" << std::endl;
      std::cout << "\tunit_cell.atom[" << i << "].hc  = " << translate(unit_cell.atom[i].hc, 0, max_hc) << ";" << std::endl;
      std::cout << "\tunit_cell.atom[" << i << "].ni  = " << 12 << ";" << std::endl;
      std::cout << "\tunit_cell.atom[" << i << "].nm  = " << false << ";" << std::endl;
   }

   return 0;

}
