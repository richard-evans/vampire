//-------------------------------------------------------
//
//  Unit Cell creator
//
//  Takes a primitive unit cell and replicates it, 
//  creating the neighbourlist and populating atomic 
//  properties for input into vampire
//
//  (C) R.F.L.Evans 22/04/2015
//
//
//-------------------------------------------------------
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

class uc_atom_t{

public:

  double cx;
  double cy;
  double cz;

  int material;
  int hc;
  int lc;

};

class material_t{

public:

  double mu_s;
  double alpha;
  double Ku;
  double Sx;
  double Sy;
  double Sz;
  std::string name;
  std::string element;


};

class nn_t{

public:
  int i;
  int j;
  int dx;
  int dy;
  int dz;
  double Jij;
};

int main(){

  // system constants
  //unit cell sizes
  double unit_cell_size[3]={3.54,3.54,3.54};

  const int num_materials=2;
  // exchange constants
  std::vector<std::vector<double> > exchange_constants;
  exchange_constants.resize(num_materials);
  for(int m=0;m<num_materials;m++) exchange_constants.at(m).resize(num_materials);

  // material parameters
  std::vector<material_t> materials(num_materials);
  
  materials.at(0).mu_s=7.63; // mu_B's
  materials.at(0).alpha=0.1;
  materials.at(0).Ku=-8.07246e-24; // J/atom
  materials.at(0).name="RE";
  materials.at(0).element="Ag";
  materials.at(0).Sx=0.0;
  materials.at(0).Sy=0.0;
  materials.at(0).Sz=-1.0;

  materials.at(1).mu_s=1.92;
  materials.at(1).alpha=0.1;
  materials.at(1).Ku=-8.07246e-24;
  materials.at(1).name="TM";
  materials.at(1).element="Fe";
  materials.at(1).Sx=0.0;
  materials.at(1).Sy=0.0;
  materials.at(1).Sz=1.0;

  exchange_constants.at(0).at(0)=1.26e-21; // Gd-Gd
  exchange_constants.at(0).at(1)=-1.09e-21; // Gd-Fe
  exchange_constants.at(1).at(0)=-1.09e-21; // Fe-Gd
  exchange_constants.at(1).at(1)=2.835e-21; // Fe-Fe
  
  // create atoms in unit cell
  std::vector<uc_atom_t> unit_cell(0);

  unit_cell.resize(4);

  unit_cell.at(0).cx=0.0;
  unit_cell.at(0).cy=0.0;
  unit_cell.at(0).cz=0.0;
  unit_cell.at(0).material=0;
  unit_cell.at(0).hc=0;
  unit_cell.at(0).lc=1;
  
  unit_cell.at(1).cx=0.5;
  unit_cell.at(1).cy=0.5;
  unit_cell.at(1).cz=0.0;
  unit_cell.at(1).material=1;
  unit_cell.at(1).hc=0;
  unit_cell.at(1).lc=0;
  
  unit_cell.at(2).cx=0.5;
  unit_cell.at(2).cy=0.0;
  unit_cell.at(2).cz=0.5;
  unit_cell.at(2).material=1;
  unit_cell.at(2).hc=1;
  unit_cell.at(2).lc=0;
  
  unit_cell.at(3).cx=0.0;
  unit_cell.at(3).cy=0.5;
  unit_cell.at(3).cz=0.5;
  unit_cell.at(3).material=1;
  unit_cell.at(3).hc=1;
  unit_cell.at(3).lc=0;
  
  // store vector of unit cells
  std::vector< std::vector < std::vector < std::vector<uc_atom_t > > > >crystal;
  
  crystal.resize(3);
  for(int i=0;i<3;i++){
    crystal.at(i).resize(3);
    for(int j=0;j<3;j++){
      crystal.at(i).at(j).resize(3);
      for(int k=0;k<3;k++){
	crystal.at(i).at(j).at(k).resize(unit_cell.size());
      }
    }
  }
    
  // replicate unit cell
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	for(int a=0;a<unit_cell.size();a++){
	  crystal.at(i).at(j).at(k).at(a).cx=unit_cell.at(a).cx+double(i);
          crystal.at(i).at(j).at(k).at(a).cy=unit_cell.at(a).cy+double(j);
	  crystal.at(i).at(j).at(k).at(a).cz=unit_cell.at(a).cz+double(k);
	  crystal.at(i).at(j).at(k).at(a).material=unit_cell.at(a).material;
	  crystal.at(i).at(j).at(k).at(a).hc=unit_cell.at(a).hc+2*k-2;
          crystal.at(i).at(j).at(k).at(a).lc=unit_cell.at(a).lc;	  
	}
      }
    }
  }

  // create neighbour list
  double nn_range=0.5*0.5+0.5*0.5;
  std::vector<nn_t> nn_list;

  // loop over all atoms in unit cell
  for(int ai=0;ai<unit_cell.size();ai++){
    double icx=crystal[1][1][1].at(ai).cx;
    double icy=crystal[1][1][1].at(ai).cy;
    double icz=crystal[1][1][1].at(ai).cz;
    int imat=crystal[1][1][1].at(ai).material;
    int ihc=crystal[1][1][1].at(ai).hc;
    int ilc=crystal[1][1][1].at(ai).lc;

    // loop over all other atoms
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	for(int k=0;k<3;k++){
	  for(int aj=0;aj<unit_cell.size();aj++){
	    double jcx=crystal[i][j][k].at(aj).cx;
	    double jcy=crystal[i][j][k].at(aj).cy;
	    double jcz=crystal[i][j][k].at(aj).cz;
	    int jmat=crystal[i][j][k].at(aj).material;
	    int jhc=crystal[i][j][k].at(aj).hc;
	    int jlc=crystal[i][j][k].at(aj).lc;
	    double range_sq=(jcx-icx)*(jcx-icx)+(jcy-icy)*(jcy-icy)+(jcz-icz)*(jcz-icz);
	    bool same_atom=(ai==aj && i==1 && j==1 && k==1);
	    if(range_sq<=nn_range && same_atom==false){
	      nn_t temp;
	      temp.i=ai;
	      temp.j=aj;
	      temp.dx=i-1;
	      temp.dy=j-1;
	      temp.dz=k-1;
	      temp.Jij=exchange_constants.at(imat).at(jmat);
	      nn_list.push_back(temp);
	      //std::cout << ai << "\t" << aj << "\t" << i-1 << "\t" << j-1 << "\t" << k-1 << "\t" << sqrt(range_sq) << std::endl;
	    } 
	  }
	}
      }
    }

  }


  // output to files
  // declare outfile file stream
  std::ofstream ucf_file;
  // open it (file_name)
  ucf_file.open ("UC.ucf");

  ucf_file << "# Unit cell size:" << std::endl;
  ucf_file << unit_cell_size[0] << "\t" << unit_cell_size[1] << "\t" << unit_cell_size[2] << std::endl;
  ucf_file << "# Unit cell vectors: " << std::endl;
  ucf_file << "1.0 0.0 0.0 " << std::endl;
  ucf_file << "0.0 1.0 0.0 " << std::endl;
  ucf_file << "0.0 0.0 1.0 " << std::endl;
  ucf_file << "# Atoms num, id cx cy cz mat lc hc " << std::endl;
  ucf_file << unit_cell.size() << std::endl;
  // loop over all atoms
  for(int atom=0; atom<unit_cell.size(); atom++){
    ucf_file << atom << "\t";
    ucf_file << unit_cell.at(atom).cx << "\t";
    ucf_file << unit_cell.at(atom).cy << "\t";
    ucf_file << unit_cell.at(atom).cz << "\t";
    ucf_file << unit_cell.at(atom).material << "\t";
    ucf_file << unit_cell.at(atom).lc << "\t";
    ucf_file << unit_cell.at(atom).hc << std::endl;
  }
  ucf_file << "#Interactions n exctype, id i j dx dy   dz        Jij"<< std::endl;
  ucf_file << nn_list.size() << "\t" << 0 << std::endl;
  // loop over all interactions
  for(unsigned int nn=0; nn < nn_list.size(); nn++){
    ucf_file << nn << "\t";
    ucf_file << nn_list[nn].i << "\t";
    ucf_file << nn_list[nn].j << "\t"; 
    ucf_file << nn_list[nn].dx << "\t";
    ucf_file << nn_list[nn].dy << "\t"; 
    ucf_file << nn_list[nn].dz << "\t"; 
    ucf_file << nn_list[nn].Jij << std::endl;
  }
  
  // material file
  std::ofstream mat_file;
  // open it (file_name)
  mat_file.open ("MAT.mat");
  mat_file << "#================================================" << std::endl;
  mat_file << "# Generated material file for input into vampire" << std::endl;
  mat_file << "#================================================" << std::endl;
  mat_file << "#" << std::endl;
  mat_file << "# File timestamp: " << std::endl;
  mat_file << "#" << std::endl;
  mat_file << "#------------------------------------------------" << std::endl;
  mat_file << "# Number of Materials" << std::endl;
  mat_file << "#------------------------------------------------" << std::endl;
  mat_file << "material:num-materials=" << num_materials << std::endl;
  mat_file << "#------------------------------------------------" << std::endl;
  
  // Loop over all materials
  for(int m=0;m<materials.size();m++){
    mat_file << "# Material " << m+1 << " (" << materials.at(m).name << ")" << std::endl;
    mat_file << "#------------------------------------------------" << std::endl;
    mat_file << "material[" << m+1 << "]:material-name=" << materials.at(m).name << std::endl;
    mat_file << "material[" << m+1 << "]:damping-constant=" << materials.at(m).alpha << std::endl;
    mat_file << "material[" << m+1 << "]:atomic-spin-moment="<< materials.at(m).mu_s << " !muB" << std::endl;
    mat_file << "material[" << m+1 << "]:uniaxial-anisotropy-constant=" << materials.at(m).Ku << std::endl;
    //mat_file << "material[" << m << "]:uniaxial-anisotropy-direction=" << 0 << "," << 1 << ","<< 0 << std::endl;
    mat_file << "material[" << m+1 << "]:material-element="<< materials.at(m).element << ""<< std::endl;
    mat_file << "material[" << m+1 << "]:initial-spin-direction="<< materials.at(m).Sx << "," << materials.at(m).Sy << "," << materials.at(m).Sz << std::endl;
    mat_file << "#------------------------------------------------" << std::endl;
  }
  return 0;
}
