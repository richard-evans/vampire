//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <array>
#include <cstdint>
#include <map>
#include <tuple>

// module headers
#include "internal.hpp"

namespace{

// dipole fields (T) written by one run of vampire
struct dipole_fields_t{
   std::vector<std::array<double,3>> atom_field; // atomistic fields, sorted by atomic position
   std::vector<std::array<double,3>> cell_field; // cell fields of cells containing atoms, in cell order
   double energy = 0.0;                          // magnetostatic energy (J) at time 0
};

// read the atomistic and cell dipole fields written by a run. The atoms are sorted by position,
// as the order of atoms in the output depends on the number of processes. Non-magnetic atoms
// (zero moment) are skipped.
bool read_dipole_fields(dipole_fields_t& fields){

   std::ifstream positions("atomistic_dipole_positions.txt");
   std::ifstream atom_fields("atomistic_dipole_field.txt");
   std::ifstream cell_fields("dipole-field.txt");
   std::ifstream output("output");
   if(!positions.is_open() || !atom_fields.is_open() || !cell_fields.is_open() || !output.is_open()) return false;

   // atomic positions rounded to 0.001 A as the sort key
   std::map<std::tuple<int64_t,int64_t,int64_t>, std::array<double,3>> sorted;
   double x, y, z, m;
   std::array<double,3> b;
   while(positions >> x >> y >> z >> m && atom_fields >> b[0] >> b[1] >> b[2]){
      if(m == 0.0) continue;
      sorted[std::make_tuple(std::llround(x*1000.0), std::llround(y*1000.0), std::llround(z*1000.0))] = b;
   }
   for(const auto& atom : sorted) fields.atom_field.push_back(atom.second);

   // cell id, number of atoms, position, moment and field
   int id, num_atoms;
   while(cell_fields >> id >> num_atoms >> x >> y >> z >> m >> b[0] >> b[1] >> b[2]){
      if(num_atoms > 0) fields.cell_field.push_back(b);
   }

   // first data line of output: time step, mx, my, mz, |m|, magnetostatic energy
   std::string line;
   while(std::getline(output, line) && (line.empty() || line[0] == '#'));
   std::stringstream data(line);
   int64_t step;
   if(!(data >> step >> x >> y >> z >> m >> fields.energy)) return false;

   return !fields.atom_field.empty() && !fields.cell_field.empty();

}

// relative rms difference sqrt( sum |b - b_ref|^2 / sum |b_ref|^2 ) between two field lists
double relative_rms_difference(const std::vector<std::array<double,3>>& b, const std::vector<std::array<double,3>>& b_ref){
   double sum_d2 = 0.0;
   double sum_ref2 = 0.0;
   for(size_t i = 0; i < b.size(); i++){
      for(int k = 0; k < 3; k++){
         sum_d2 += (b[i][k] - b_ref[i][k])*(b[i][k] - b_ref[i][k]);
         sum_ref2 += b_ref[i][k]*b_ref[i][k];
      }
   }
   return std::sqrt(sum_d2/sum_ref2);
}

// remove files written by a run
void cleanup(){
   vt::remove({"output", "log", "atomistic_dipole_positions.txt", "atomistic_dipole_field.txt", "dipole-field.txt"});
}

}

//------------------------------------------------------------------------------
// Test that the tensor and hierarchical dipole solvers agree with the exact
// atomistic solver for a two atom thick cylinder (sc, a = 2.87 A, mu_s = 2.2 muB,
// 1 nm cells) magnetised along (1,0,1), e.g. dipole/cylinder (10 nm diameter,
// 1898 atoms), dipole/cylinder-fill (the same in a kept non-magnetic fill) or
// dipole/disc (50 nm diameter, 47698 atoms).
//
//    input-atomistic, input-tensor, input-hierarchical
//
// The relative rms differences from the atomistic fields must be below
// atom_tol for the atomistic fields and cell_tol for the cell fields (the
// atomistic solver's cell field is the mean field of the atoms in the cell),
// and the hierarchical cell fields must agree with the tensor cell fields
// within hierarchical_tol. The mean atomistic field must give demagnetising
// factors with Nxx + Nyy + Nzz = 1 within 0.01, where Nxx = 1 - <Bx>/(mu0 Mx).
// The magnetostatic energy must equal (1/2) mu0 Ms^2 V m.N.m with these
// factors, within 1e-6 for the atomistic solver (where it is an identity) and
// within energy_tol for the tensor and hierarchical solvers.
//------------------------------------------------------------------------------
bool dipole_solver_test(const std::string dir, double atom_tol, double cell_tol, double hierarchical_tol, double energy_tol, const std::string executable){

   // get root directory
   std::string path = std::filesystem::current_path();

   // fixed-width output for prettiness
   std::stringstream test_name;
   test_name << "Testing dipole solvers for " << dir;
   std::cout << std::setw(60) << std::left << test_name.str() << " : " << std::flush;

   // change directory
   if( !vt::chdir(path+"/data/"+dir) ) return false;
   cleanup();

   // run each solver and read the fields
   const std::vector<std::string> solvers = {"atomistic", "tensor", "hierarchical"};
   std::vector<dipole_fields_t> fields(solvers.size());
   std::string error;
   for(size_t s = 0; s < solvers.size() && error.empty(); s++){
      if(vt::system(executable + " --input-file input-" + solvers[s]) != 0) error = "error running vampire (" + solvers[s] + ")";
      else if(!read_dipole_fields(fields[s])) error = "missing dipole field output (" + solvers[s] + ")";
      cleanup();
   }

   // return to parent directory
   if( !vt::chdir(path) ) return false;

   if(error.empty()){
      for(size_t s = 1; s < solvers.size(); s++){
         if(fields[s].atom_field.size() != fields[0].atom_field.size() || fields[s].cell_field.size() != fields[0].cell_field.size()){
            error = "different number of atoms or cells (" + solvers[s] + ")";
         }
      }
   }
   if(!error.empty()){
      std::cout << "FAIL | " << error << std::endl;
      return false;
   }

   // demagnetising factors from the mean atomistic field, mu0 Ms = mu0 mu_s / a^3 and m = (1,0,1)/sqrt(2)
   const double mu0_ms = 4.0e-7*M_PI * 2.2*9.2740100783e-24 / (2.87e-10*2.87e-10*2.87e-10);
   const double m = 1.0/std::sqrt(2.0);
   std::array<double,3> mean_b = {0.0, 0.0, 0.0};
   for(const auto& b : fields[0].atom_field) for(int k = 0; k < 3; k++) mean_b[k] += b[k]/double(fields[0].atom_field.size());
   const double nxx = 1.0 - mean_b[0]/(mu0_ms*m);
   const double nzz = 1.0 - mean_b[2]/(mu0_ms*m);
   const double trace = 2.0*nxx + nzz; // Nyy = Nxx by symmetry

   // magnetostatic energy (1/2) mu0 Ms^2 V (Nxx mx^2 + Nzz mz^2), with V = number of atoms * a^3
   const double volume = double(fields[0].atom_field.size()) * 2.87e-10*2.87e-10*2.87e-10;
   const double energy = 0.5 * mu0_ms*mu0_ms/(4.0e-7*M_PI) * volume * (nxx*m*m + nzz*m*m);
   std::vector<double> energy_error(solvers.size());
   for(size_t s = 0; s < solvers.size(); s++) energy_error[s] = std::abs(fields[s].energy/energy - 1.0);

   const double tensor_atom = relative_rms_difference(fields[1].atom_field, fields[0].atom_field);
   const double hierarchical_atom = relative_rms_difference(fields[2].atom_field, fields[0].atom_field);
   const double tensor_cell = relative_rms_difference(fields[1].cell_field, fields[0].cell_field);
   const double hierarchical_cell = relative_rms_difference(fields[2].cell_field, fields[0].cell_field);
   const double hierarchical_tensor = relative_rms_difference(fields[2].cell_field, fields[1].cell_field);

   if(tensor_atom < atom_tol && hierarchical_atom < atom_tol && tensor_cell < cell_tol && hierarchical_cell < cell_tol &&
      hierarchical_tensor < hierarchical_tol && std::abs(trace - 1.0) < 0.01 &&
      energy_error[0] < 1.0e-6 && energy_error[1] < energy_tol && energy_error[2] < energy_tol){
      std::cout << "OK" << std::endl;
      return true;
   }
   else{
      std::cout << "FAIL | rms difference from atomistic: atoms (tensor, hierarchical) " << tensor_atom << ", " << hierarchical_atom
                << "; cells " << tensor_cell << ", " << hierarchical_cell << "; hierarchical-tensor cells " << hierarchical_tensor
                << "; Nxx + Nyy + Nzz = " << trace << "; relative energy error (atomistic, tensor, hierarchical) "
                << energy_error[0] << ", " << energy_error[1] << ", " << energy_error[2] << std::endl;
      return false;
   }

}
