#include <iostream>
#include <vector>
#include <cmath>

class atom_t{

public:
  double x;
  double y;
  double z;
  int id;
  int no;
};

class nn_t{

public:
  int i;
  int j;
  int dx;
  int dy;
  int dz;

};
    
int main(){

  // size in unit cells
  int nx=5,ny=5,nz=5;

  // atoms
  std::vector<atom_t> cellatoms(2);

  // bcc
  cellatoms[0].x=0.0;
  cellatoms[0].y=0.0;
  cellatoms[0].z=0.0;
  cellatoms[0].id=0;

  cellatoms[1].x=0.5;
  cellatoms[1].y=0.5;
  cellatoms[1].z=0.5;
  cellatoms[1].id=1;

  // cutoff radius (fraction of a)
  double cutoff=sqrt(3.0)/2.0+0.05; 

  int counter=0;

  std::vector<std::vector <std::vector <std::vector <atom_t> > > > atoms(0);
  std::vector<nn_t> neighbours(0);

  // generate lattice
  atoms.resize(nx);
  for(int i=0;i<nx;i++){
    atoms[i].resize(ny);
    for(int j=0;j<ny;j++){
      atoms[i][j].resize(nz);
      for(int k=0;k<nz;k++){
	for(int atom=0;atom<cellatoms.size();atom++){
	  atom_t tmpatom;
	  tmpatom.x=cellatoms[atom].x+double(i);
          tmpatom.y=cellatoms[atom].y+double(j);
          tmpatom.z=cellatoms[atom].z+double(k);
	  tmpatom.id=atom;
	  tmpatom.no=counter;
	  counter++;
	  atoms[i][j][k].push_back(tmpatom);
	  //std::cout << atoms[i][j][k][atom].x << " " << atoms[i][j][k][atom].y << " " << atoms[i][j][k][atom].z << " " << atoms[i][j][k][atom].id << " " << atoms[i][j][k][atom].no << " " <<std::endl;
	}
      }
    }
  }


  //std::cout << counter << " atoms generated" << std::endl;

  // calculate neighbourlist for central cell and output to screen
  int cx=nx/2+1;
  int cy=ny/2+1;
  int cz=nz/2+1;

  // loop over all atoms in local unit cell
  for(int atom=0;atom<cellatoms.size();atom++){
    double ix=atoms[cx][cy][cz][atom].x;
    double iy=atoms[cx][cy][cz][atom].y;
    double iz=atoms[cx][cy][cz][atom].z;
    int iid=atoms[cx][cy][cz][atom].id;
    int in=atoms[cx][cy][cz][atom].no;
    //std::cout << ix << " " << iy << " " << iz << " " << iid << " " << in << std::endl;
    // loop over all other unit cells
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
	for(int k=0;k<nz;k++){
	  for(int jatom=0;jatom<cellatoms.size();jatom++){
	    double jx=atoms[i][j][k][jatom].x;
	    double jy=atoms[i][j][k][jatom].y;
	    double jz=atoms[i][j][k][jatom].z;
	    int jid=atoms[i][j][k][jatom].id;
	    int jn=atoms[i][j][k][jatom].no;
	    //std::cout << jx << " " << jy << " " << jz << " " << jid << " " << jn << std::endl;
	    
	    // exclude self-interaction
	    if(in!=jn){
	      double range=sqrt((jx-ix)*(jx-ix)+(jy-iy)*(jy-iy)+(jz-iz)*(jz-iz));
	      if(range<cutoff){
		//std::cout << iid << "\t" << jid << "\t" << i-cx << "\t" << j-cy << "\t" << k-cz << std::endl;
		nn_t tmpnn;
		tmpnn.i=iid;
		tmpnn.j=jid;
		tmpnn.dx=i-cx;
		tmpnn.dy=j-cy;
                tmpnn.dz=k-cz;
		neighbours.push_back(tmpnn);
	      }
	    }
	  }
	}
      }
    }
  }

  // output data to screen
  std::cout << "\t\t\tunit_cell.lcsize=" << cellatoms.size() << ";" << std::endl;
  std::cout << "\t\t\tunit_cell.hcsize=" << 2 << ";" << std::endl;
  std::cout << "\t\t\tunit_cell.interaction_range=" << 1 << ";" << std::endl;
  std::cout << "\t\t\tunit_cell.atom.resize(" << cellatoms.size() << ");" << std::endl;
  for(int atom=0;atom<cellatoms.size();atom++){
    std::cout << "\t\t\t//-----------------------------" << std::endl;
    std::cout << "\t\t\tunit_cell.atom[" << atom << "].x=" << cellatoms[atom].x << ";"<< std::endl;
    std::cout << "\t\t\tunit_cell.atom[" << atom << "].y=" << cellatoms[atom].y << ";"<< std::endl;
    std::cout << "\t\t\tunit_cell.atom[" << atom << "].z=" << cellatoms[atom].z << ";"<< std::endl;
    std::cout << "\t\t\tunit_cell.atom[" << atom << "].lc=" << atom << ";"<< std::endl;
    std::cout << "\t\t\tunit_cell.atom[" << atom << "].hc=" << atom << ";"<< std::endl;
  }
  std::cout << "\t\t\t//-----------------------------" << std::endl;
  std::cout << "\t\t\tunit_cell.interaction.resize(" << neighbours.size() << ");" << std::endl;
  
  for(int nn=0;nn<neighbours.size();nn++){
    std::cout << "\t\t\t//-----------------------------" << std::endl;
  std::cout << "\t\t\tunit_cell.interaction["<< nn <<"].i="<< neighbours[nn].i<<";" << std::endl;
  std::cout << "\t\t\tunit_cell.interaction["<< nn <<"].j="<< neighbours[nn].j<<";" << std::endl;
  std::cout << "\t\t\tunit_cell.interaction["<< nn <<"].dx="<< neighbours[nn].dx<<";" << std::endl;
  std::cout << "\t\t\tunit_cell.interaction["<< nn <<"].dy="<< neighbours[nn].dy<<";" << std::endl;
  std::cout << "\t\t\tunit_cell.interaction["<< nn <<"].dz="<< neighbours[nn].dz<<";" << std::endl;
}
  return 0; 

}
