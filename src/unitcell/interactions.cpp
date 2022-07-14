//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016, Jack B Collings 2021. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>

// Vampire headers
#include "errors.hpp"
#include "unitcell.hpp"
#include "vio.hpp"

// unitcell module headers
#include "internal.hpp"

namespace local{

   class atom_t{

   public:
      int mat; // material
      double x; // positions
      double y;
      double z;
      int id; // atom number in unit cell
      int idx; // unit cell number
      int idy;
      int idz;
      bool nm;
   };

}

namespace unitcell{
namespace internal{

//------------------------------------------------------------------------------
//  Function to calculate neighbour interactions assuming fractional unit cell
//------------------------------------------------------------------------------
void calculate_interactions(unit_cell_t& unit_cell){

   // Resize material-exchange-nn-cutoff tensor
   unsigned int num_uc_materials = 0;
   for (size_t i = 0; i < unit_cell.atom.size(); ++i){
         if (unit_cell.atom[i].mat > num_uc_materials) num_uc_materials = unit_cell.atom[i].mat;
      }
   ++num_uc_materials; // since unit cell category has the -1 shift
   nn_cutoff_range.resize(num_uc_materials, std::vector<double>(num_uc_materials));
   interaction_cutoff_range.resize(num_uc_materials, std::vector<double>(num_uc_materials));
   if (exchange_function == exponential || exchange_function == material_exponential || exchange_function == RKKY){
      material_exchange_parameters.resize(num_uc_materials, std::vector<exchange_parameters_t>(num_uc_materials));
   }
   else material_exchange_parameters.clear();

   // Set nn and interaction ranges using cutoff factors
   double max_rcut = 0; // Max rcut used to find number of unit cells to get all interactions/nearest neighbours
   for (unsigned int i = 0; i < num_uc_materials; ++i){
      for (unsigned int j = 0; j < num_uc_materials; ++j){
         nn_cutoff_range[i][j] *= 1.000001*unit_cell.cutoff_radius;
         interaction_cutoff_range[i][j] *= nn_cutoff_range[i][j]*exchange_interaction_range;
         double tmp_r = std::max(nn_cutoff_range[i][j], interaction_cutoff_range[i][j]);
         if (tmp_r > max_rcut) max_rcut = tmp_r;
      }
   }
   if (exchange_function == internal::shell || exchange_function == internal::nearest_neighbour) nn_cutoff_range.clear();

   // temporary for unit cell size
   const double ucsx = unit_cell.dimensions[0];
   const double ucsy = unit_cell.dimensions[1];
   const double ucsz = unit_cell.dimensions[2];

   // save number of atoms in unit cell
   unit_cell.bilinear.num_unit_cell_atoms = unit_cell.atom.size();
   unit_cell.biquadratic.num_unit_cell_atoms = unit_cell.atom.size();

   // set number of interactions per atom used for exchange normalisation
   for (size_t i = 0; i < unit_cell.atom.size(); ++i){
      unit_cell.bilinear.ni.push_back(unit_cell.atom[i].ni);
      unit_cell.biquadratic.ni.push_back(unit_cell.atom[i].ni);
   }

   // determine number of unit cells in x,y and z
   const int nx = 1 + 2*ceil(max_rcut); // number of replicated cells in x,y,z
   const int ny = 1 + 2*ceil(max_rcut);
   const int nz = 1 + 2*ceil(max_rcut);

   zlog << zTs() << "Generating neighbour interactions for a lattice of " << nx << " x " << ny << " x " << nz << " unit cells for neighbour calculation" << std::endl;

   // temporary array for replicated atoms
   std::vector<local::atom_t> ratoms;

   // replicate in x,y,z
   for(int x = 0; x < nx; ++x){
      for(int y = 0; y < ny; ++y){
         for(int z = 0; z < nz; ++z){
            for(size_t a=0; a < unit_cell.atom.size(); ++a){
               local::atom_t tmp;
               tmp.mat = unit_cell.atom[a].mat;
               tmp.x = (unit_cell.atom[a].x + double(x))*ucsx;
               tmp.y = (unit_cell.atom[a].y + double(y))*ucsy;
               tmp.z = (unit_cell.atom[a].z + double(z))*ucsz;
               tmp.id = a;
               tmp.idx = x; // unit cell id
               tmp.idy = y;
               tmp.idz = z;
               tmp.nm = unit_cell.atom[a].nm;
               ratoms.push_back(tmp);
            }
         }
      }
   }

   // calculate neighbours and exchange for central cell
   int mid_cell_x = (nx-1)/2;
   int mid_cell_y = (ny-1)/2;
   int mid_cell_z = (nz-1)/2;

   //const double nnrcut_sq = unit_cell.cutoff_radius*unit_cell.cutoff_radius*1.001*1.001; // nearest neighbour cutoff radius

   // loop over all i atoms
   for(size_t i=0; i < ratoms.size(); ++i){

      // check for i atoms only in central cell
      if( (ratoms[i].idx == mid_cell_x) && (ratoms[i].idy == mid_cell_y) && (ratoms[i].idz == mid_cell_z) ){

         // loop over all j atoms
         for(size_t j=0; j < ratoms.size(); ++j){

            // calculate interatomic radius_sq
            const double rx = ratoms[j].x - ratoms[i].x;
            const double ry = ratoms[j].y - ratoms[i].y;
            const double rz = ratoms[j].z - ratoms[i].z;
            double range_sq = rx*rx + ry*ry + rz*rz;
            double range = sqrt(range_sq);
            bool magnetic = !(ratoms[i].nm || ratoms[j].nm);
            // check for rij < rcut and i!= j
            if( range < interaction_cutoff_range[ratoms[i].mat][ratoms[j].mat] && i != j && magnetic ){
               // Neighbour found
               uc::interaction_t tmp;

               // Determine unit cell id for i and j atoms
               tmp.i = ratoms[i].id;
               tmp.j = ratoms[j].id;

               // Determine unit cell material for i and j atoms
               tmp.mat_i = ratoms[i].mat;
               tmp.mat_j = ratoms[j].mat;

               // Determine unit cell offsets
               tmp.dx = ratoms[j].idx - ratoms[i].idx;
               tmp.dy = ratoms[j].idy - ratoms[i].idy;
               tmp.dz = ratoms[j].idz - ratoms[i].idz;

               // save interaction range
               tmp.rij = range;

               // save initial shell
               tmp.shell = 0;

               // Determine normalised exchange constants
               tmp.Jij[0][0] = uc::internal::exchange(range, interaction_cutoff_range[ratoms[i].mat][ratoms[j].mat], ratoms[i].mat, ratoms[j].mat); // xx
               tmp.Jij[0][1] = 0.0; // xy
               tmp.Jij[0][2] = 0.0; // xz

               tmp.Jij[1][0] = 0.0; // yx
               tmp.Jij[1][1] = uc::internal::exchange(range, interaction_cutoff_range[ratoms[i].mat][ratoms[j].mat], ratoms[i].mat, ratoms[j].mat); // yy
               tmp.Jij[1][2] = 0.0; // yz

               tmp.Jij[2][0] = 0.0; // zx
               tmp.Jij[2][1] = 0.0; // zy
               tmp.Jij[2][2] = uc::internal::exchange(range, interaction_cutoff_range[ratoms[i].mat][ratoms[j].mat], ratoms[i].mat, ratoms[j].mat); // zz

               unit_cell.bilinear.interaction.push_back(tmp);

               // save same interactions for biquadratic exchange
               unit_cell.biquadratic.interaction.push_back(tmp);

            }
         }
      }
   }

   // Set calculated interactions range
   int interaction_range=0;
   for(size_t i=0; i<unit_cell.bilinear.interaction.size(); i++){
      if(abs(unit_cell.bilinear.interaction[i].dx)>interaction_range) interaction_range=abs(unit_cell.bilinear.interaction[i].dx);
      if(abs(unit_cell.bilinear.interaction[i].dy)>interaction_range) interaction_range=abs(unit_cell.bilinear.interaction[i].dy);
      if(abs(unit_cell.bilinear.interaction[i].dz)>interaction_range) interaction_range=abs(unit_cell.bilinear.interaction[i].dz);
   }
   unit_cell.interaction_range = interaction_range;

   // Normalise exchange interactions

   unit_cell.bilinear.normalise_exchange(nn_cutoff_range);
   unit_cell.biquadratic.normalise_exchange(nn_cutoff_range);

   // Find shells for neighbours
   unit_cell.bilinear.find_shells();

   // Output interactions to screen
   /*for(int i=0; i<unit_cell.bilinear.interaction.size(); i++){
      std::cerr << i << "\t" << unit_cell.bilinear.interaction[i].i << "\t"
                << unit_cell.bilinear.interaction[i].j << "\t"
                << unit_cell.bilinear.interaction[i].dx << "\t"
                << unit_cell.bilinear.interaction[i].dy << "\t"
                << unit_cell.bilinear.interaction[i].dz << "\t"
                << unit_cell.bilinear.interaction[i].rij << "\t"
                << unit_cell.bilinear.interaction[i].Jij[0][0] << std::endl;
   }*/

   // Output interactions including shell ID and lattice vectors to screen
   /*for(int i=0; i<unit_cell.bilinear.interaction.size(); i++){
      // Determine unit cell id for i and j atoms
      int idi = unit_cell.bilinear.interaction[i].i;
      int idj = unit_cell.bilinear.interaction[i].j;

      double jx = unit_cell.atom[idj].x + double(unit_cell.bilinear.interaction[i].dx);
      double jy = unit_cell.atom[idj].y + double(unit_cell.bilinear.interaction[i].dy);
      double jz = unit_cell.atom[idj].z + double(unit_cell.bilinear.interaction[i].dz);

      // Determine unit cell offsets
      double dx = jx - unit_cell.atom[idi].x;
      double dy = jy - unit_cell.atom[idi].y;
      double dz = jz - unit_cell.atom[idi].z;

      std::cerr << i << "\t" << unit_cell.bilinear.interaction[i].i << "\t"
                << unit_cell.bilinear.interaction[i].j << "\t"
                << unit_cell.bilinear.interaction[i].shell << "\t"
                << dx << "\t"
                << dy << "\t"
                << dz << "\t"
                << unit_cell.bilinear.interaction[i].rij << "\t"
                << unit_cell.bilinear.interaction[i].Jij[0][0] << std::endl;
   }*/

   // Check for interactions
   if(unit_cell.bilinear.interaction.size()==0){
      terminaltextcolor(RED);
      std::cerr << "Error! No interactions generated for " << uc::internal::crystal_structure << " crystal structure. Try increasing the interaction range. Aborting." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error! No interactions generated for " << uc::internal::crystal_structure << " crystal structure. Try increasing the interaction range. Aborting." << std::endl;
      err::vexit();
   }

   return;

}


   } // end of internal namespace

} // end of unitcell namespace
