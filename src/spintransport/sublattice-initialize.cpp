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

// Vampire headers
#include "create.hpp"
#include "spintransport.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{

   namespace internal{

      //--------------------------------------------------------------------------------
      // Function to generate cell data in the case of sublattice resolved torques
      //--------------------------------------------------------------------------------
      void initialise_sublattice_cell_data(
         const uint64_t num_atoms,   // number of local atoms
         const int num_materials, // number of materials
         const std::vector<int>& atoms_type_array, // material types of atoms
         const std::vector<double>& atoms_m_spin_array,  // moments of atoms (muB)
         const std::vector<double>& material_damping_array, // array of material level damping constants
         const std::vector<bool>& is_magnetic_material, // array of size num_mat to state whether material is magnetic (true) or not (false)
         const std::vector<cs::nm_atom_t>& non_magnetic_atoms_array, // list of non-magnetic atoms
         internal::v3cell3D_t& cells3D, // Cell data
         const int stack_x, const int stack_y, const int stack_z // direction of stacks relative to current direction

      ){

         // copy cell sizes to local variable
         const double cell_size[3] = { st::internal::cell_size_x,
                                       st::internal::cell_size_y,
                                       st::internal::cell_size_z };

         //---------------------------------------------------------------------------------
         // Determine number of unique sublattices and atom associations
         //---------------------------------------------------------------------------------
         std::vector<int> sublattice_list; // list of unique sublattices to calculate ST for
         // get list of unique sublattices to calculate ST for
         for(int mat = 0; mat < num_materials; mat++){
            int stsl = st::internal::mp[mat].stsl;
            if( std::find(sublattice_list.begin(), sublattice_list.end(), stsl) == sublattice_list.end() ){
               sublattice_list.push_back(stsl);
            }
         }

         //for(int i=0; i< sublattice_list.size(); i++){
         //   std::cout << "Sublattice " << sublattice_list[i] << " mapped to index " << i << std::endl;
         //}

         // relabel list in linear order for memory efficiency
         for(int mat = 0; mat < num_materials; mat++){
            for(int stsl = 0; stsl < sublattice_list.size(); stsl++){
               if(st::internal::mp[mat].stsl == sublattice_list[stsl]){
                  st::internal::mp[mat].stsl = stsl;
               }
            }
         }

         //for(int mat = 0; mat < num_materials; mat++){
         //   std::cout << "Material " << mat << ": stsl = " << st::internal::mp[mat].stsl << std::endl;
         //}

         // now associate each atom with the correct sublattice for spin transport calculation
         st::internal::atom_sublattice.resize(num_atoms);
         for(uint64_t atom = 0; atom < num_atoms; atom++){
            int mat = atoms_type_array[atom];
            st::internal::atom_sublattice[atom] = st::internal::mp[mat].stsl;
         }

         // set num sublattices to number of unique sublattices in system
         st::internal::num_sublattices = sublattice_list.size();

         // define local constant for number of sublattices to avoid repeated access to global variable
         const unsigned int num_sublattices = st::internal::num_sublattices;

         //---------------------------------------------------------------------------------
         // Initialise suplattice specific calculation of spin-transport
         //---------------------------------------------------------------------------------

         // get total number of cells
         const int num_cells = st::internal::total_num_cells;

         // resize boolean array to determine if cell is magnetic or not
         st::internal::sl_magnetic.resize(num_cells);

         // resize arrays to store average resistance in each cell
         st::internal::cell_sl_resistance.resize(num_cells); // also stores accumulated resistance before averaging
         st::internal::cell_sl_spin_resistance.resize(num_cells); // also stores accumulated spin resistance before averaging

         // resize array to store average damping constant alpha in each cell
         st::internal::cell_sl_alpha.resize(num_cells);

         // resize array to store inverse total moment m_s (T = 0)
         st::internal::cell_sl_isaturation.resize(num_cells); // also stores accumulated moment before inversion

         // resize arrays to store slonczewski prefactors
         st::internal::cell_sl_relaxation_torque_rj.resize(num_cells); // also stores accumulated parameters before averaging
         st::internal::cell_sl_precession_torque_pj.resize(num_cells); // also stores accumulated parameters before averaging

         // loop over all sublattices to resize into num_cells
         for( int cell = 0; cell < num_cells; cell++ ){
            st::internal::sl_magnetic[cell].resize(                 num_sublattices,true); // initially assume all cells magnetic
            st::internal::cell_sl_resistance[cell].resize(          num_sublattices,0.0);
            st::internal::cell_sl_spin_resistance[cell].resize(     num_sublattices,0.0);
            st::internal::cell_sl_alpha[cell].resize(               num_sublattices,0.0); // also stores accumulated damping constant before averaging
            st::internal::cell_sl_isaturation[cell].resize(         num_sublattices,1.0); // also stores accumulated moment before inversion
            st::internal::cell_sl_relaxation_torque_rj[cell].resize(num_sublattices,0.0); // also stores accumulated parameters before averaging
            st::internal::cell_sl_precession_torque_pj[cell].resize(num_sublattices,0.0); // also stores accumulated parameters before averaging
         }

         // Temporary arrays to accumulate intermediate values
         std::vector<uint64_t> total_num_atoms_in_cell( num_cells, 0 ); // array to store total number of magnetic and non-magnetic atoms in cell
         //std::vector<double>   total_num_magnetic_atoms(   st::internal::total_num_cells, 0  ); // array to store total number of actually magnetic atoms (ignoring nm = keep atoms) [double format for easy normalisation]
         //std::vector<double>   total_resistivity_sq(       st::internal::total_num_cells, 0.0); // array to store total resistivity^2 calculated from constituent atoms
         //std::vector<double>   total_spin_resistivity_sq(  st::internal::total_num_cells, 0.0); // array to store total spin resistivity^2 calculated from constituent atoms

         // Temporary arrays to accumulate intermediate values for material-resolved calculation
         std::vector < std::vector<uint64_t> > cell_sl_num_atoms_in_cell(  num_cells);   // array to store total number of magnetic and non-magnetic atoms in cell
         std::vector < std::vector< double > > cell_sl_num_magnetic_atoms( num_cells);  // array to store total number of actually magnetic atoms (ignoring nm = keep atoms) [double format for easy normalisation]
         //std::vector < std::vector< double > > cell_sl_resistivity(        num_cells);         // array to store total resistivity calculated from constituent atoms
         //std::vector < std::vector< double > > cell_sl_resistivity_sq(     num_cells);      // array to store total resistivity^2 calculated from constituent atoms
         //std::vector < std::vector< double > > cell_sl_spin_resistivity(   num_cells);    // array to store total spin resistivity calculated from constituent atoms
         //std::vector < std::vector< double > > cell_sl_spin_resistivity_sq(num_cells); // array to store total spin resistivity^2 calculated from constituent atoms*/

         // allocate memory for sublattice properties for each cell
         for( int cell = 0; cell < num_cells; cell++ ){
            cell_sl_num_atoms_in_cell[cell].resize(  num_sublattices,  0  );
            cell_sl_num_magnetic_atoms[cell].resize( num_sublattices, 0.0 );
            //cell_sl_resistivity[cell].resize(        num_sublattices, 0.0 );
            //cell_sl_resistivity_sq[cell].resize(     num_sublattices, 0.0 );
            //cell_sl_spin_resistivity[cell].resize(   num_sublattices, 0.0 );
            //cell_sl_spin_resistivity_sq[cell].resize(num_sublattices, 0.0 );
         }

         std::vector<uint64_t> sl_num_atoms_in_cell(num_sublattices);   // array to store total number of magnetic and non-magnetic atoms in cell
         std::vector< double > sl_num_magnetic_atoms(num_sublattices);  // array to store total number of actually magnetic atoms (ignoring nm = keep atoms)
         std::vector< double > sl_resistivity(num_sublattices); // array to store total resistivity calculated from constituent atoms
         //std::vector< double > sl_resistivity_sq(num_sublattices);      // array to store total resistivity^2 calculated from constituent atoms
         std::vector< double > sl_spin_resistivity(num_sublattices); // array to store total spin resistivity calculated from constituent atoms
         //std::vector< double > sl_spin_resistivity_sq(num_sublattices); // array to store total spin resistivity^2 calculated from constituent atoms

         // sublattice arrays to store total moments, alpha constants, and STT prefactors
         std::vector< double > sl_total_moment(num_sublattices, 0.0);          // array to store total moment in each cell for each sublattice
         std::vector< double > sl_total_alpha(num_sublattices, 0.0);          // array to store total moment in each cell for each sublattice
         std::vector< double > sl_total_rj(num_sublattices, 0.0);          // array to store total moment in each cell for each sublattice
         std::vector< double > sl_total_pj(num_sublattices, 0.0);          // array to store total moment in each cell for each sublattice

         // loop over all xy-cells (stacks)
         for(unsigned int i = 0; i < cells3D.size(); i++){
            for(unsigned int j = 0; j < cells3D[i].size(); j++){

               // loop over all cells in stack
               for(unsigned int k = 0; k < cells3D[i][j].size(); k++){

                  // determine cell ID
                  const uint64_t cell = cells3D[i][j][k].id;

                  // determine total number of local atoms
                  const uint64_t num_atoms_in_cell = cells3D[i][j][k].atom.size() + cells3D[i][j][k].nm_atom.size();
                  total_num_atoms_in_cell[cell] = num_atoms_in_cell;

                  // variables to accumulate resistances for cell
                  for(int sl = 0 ; sl < num_sublattices; sl++) sl_num_atoms_in_cell  [sl] =  0 ; // array to store total number of magnetic and non-magnetic atoms in cell
                  for(int sl = 0 ; sl < num_sublattices; sl++) sl_num_magnetic_atoms [sl] = 0.0; // array to store total number of actually magnetic atoms (ignoring nm = keep atoms)
                  for(int sl = 0 ; sl < num_sublattices; sl++) sl_resistivity        [sl] = 0.0; // array to store total resistivity calculated from constituent atoms
                  //for(int sl = 0 ; sl < num_sublattices; sl++) sl_resistivity_sq     [sl] = 0.0; // array to store total resistivity^2 calculated from constituent atoms
                  for(int sl = 0 ; sl < num_sublattices; sl++) sl_spin_resistivity   [sl] = 0.0; // array to store total spin resistivity calculated from constituent atoms
                  //for(int sl = 0 ; sl < num_sublattices; sl++) sl_spin_resistivity_sq[sl] = 0.0; // array to store total spin resistivity^2 calculated from constituent atoms

                  for(int sl = 0 ; sl < num_sublattices; sl++) sl_total_moment[sl] = 0.0;
                  for(int sl = 0 ; sl < num_sublattices; sl++) sl_total_alpha[sl]  = 0.0; // damping constant
                  for(int sl = 0 ; sl < num_sublattices; sl++) sl_total_rj[sl]     = 0.0; // spin torque prefactor parameter
                  for(int sl = 0 ; sl < num_sublattices; sl++) sl_total_pj[sl]     = 0.0; // spin torque prefactor parameter

                  // check for cells with only non-magnetic atoms = keep
                  //uint64_t num_magnetic_atoms = 0; // counter for number of actually magnetic atoms

                  // magnetic atoms
                  for(unsigned int atom = 0; atom < cells3D[i][j][k].atom.size(); atom++ ){
                     // get atom number in total list
                     const int atomID = cells3D[i][j][k].atom[atom];
                     // get atom material of atom
                     const int mat = atoms_type_array[atomID]; // get material type
                     // get sublattice of atom
                     const int sl = st::internal::atom_sublattice[atom];

                     // get resistivities
                     const double r  = st::internal::mp[mat].resistivity.get();
                     const double sr = st::internal::mp[mat].spin_resistivity.get();

                     // add resistivities
                     sl_resistivity[sl]         += r;     // add resistivity to total
                     //sl_resistivity_sq[sl]      += r*r;   // add resistivity^2 to total
                     sl_spin_resistivity[sl]    += sr;    // add spin resistivity to total
                     //sl_spin_resistivity_sq[sl] += sr*sr; // add spin resistivity^2 to total

                     // add number of total atoms (magnetic and non-magnetic) for each sublattice for each cell
                     sl_num_atoms_in_cell[sl] += 1;

                     // check that atom is magnetic
                     if(is_magnetic_material[mat] == true){

                        // increment number of atoms in sublattice counter
                        sl_num_magnetic_atoms[sl] += 1.0;

                        // calculate total moment
                        sl_total_moment[sl] += atoms_m_spin_array[atomID];
                        sl_total_alpha [sl] += material_damping_array[mat];

                        // calculate total prefactor parameters
                        sl_total_rj[sl] += st::internal::mp[mat].stt_rj.get();
                        sl_total_pj[sl] += st::internal::mp[mat].stt_pj.get();

                     }

                  }

                  // non-magnetic (remove) atoms
                  for(unsigned int atom = 0; atom < cells3D[i][j][k].nm_atom.size(); atom++ ){

                     const int mat = non_magnetic_atoms_array[atom].mat; // get material type
                     const int sl = st::internal::atom_sublattice[atom]; // get sublattice of atom

                     // add number of total atoms (magnetic and non-magnetic) for each sublattice for each cell
                     sl_num_atoms_in_cell[sl] += 1;

                     // get resistivities
                     const double r  = st::internal::mp[mat].resistivity.get();
                     const double sr = st::internal::mp[mat].spin_resistivity.get();

                     sl_resistivity[sl]         += r;     // add resistivity to total
                     //sl_resistivity_sq[sl]      += r*r;   // add resistivity^2 to total
                     sl_spin_resistivity[sl]    += sr;    // add spin resistivity to total
                     //sl_spin_resistivity_sq[sl] += sr*sr; // add spin resistivity^2 to total

                  }

                  // save accumulated values
                  for( int sl = 0; sl < num_sublattices; sl++ ){
                     st::internal::cell_sl_resistance          [cell][sl] = sl_resistivity     [sl]; // temporarily store rho_eff for each sublattice here
                     st::internal::cell_sl_spin_resistance     [cell][sl] = sl_spin_resistivity[sl]; // temporarily store spin-dependent rho_eff for each sublattice here
                     st::internal::cell_sl_relaxation_torque_rj[cell][sl] = sl_total_rj        [sl];
                     st::internal::cell_sl_precession_torque_pj[cell][sl] = sl_total_pj        [sl];
                     st::internal::cell_sl_isaturation         [cell][sl] = sl_total_moment    [sl]; // (mu_B)
                     st::internal::cell_sl_alpha               [cell][sl] = sl_total_alpha     [sl]; // total damping

                     // temporary values
                     cell_sl_num_atoms_in_cell[cell][sl]   = sl_num_atoms_in_cell[sl];   // store total number of magnetic and non-magnetic atoms in cell
                     cell_sl_num_magnetic_atoms[cell][sl]  = sl_num_magnetic_atoms[sl];  // store total number of magnetic atoms in cell in double format for normalisation
                     //cell_sl_resistivity[sl][cell]         = sl_resistivity[sl];
                     //cell_sl_resistivity_sq[cell][sl]      = sl_resistivity_sq[sl];      // store total resistivity^2 calculated from constituent atoms
                     //cell_sl_spin_resistivity[sl][cell]    = sl_spin_resistivity[sl];
                     //cell_sl_spin_resistivity_sq[cell][sl] = sl_spin_resistivity_sq[sl]; // store total spin resistivity^2 calculated from constituent atoms
                  }

               } // end of cell loop
            } // end of stack y loop
         } // end of stack x loop

         //-----------------------------------------------------------------------------------------------
         // reduce on all processors to enable unique determination of empty cells and average parameters
         //-----------------------------------------------------------------------------------------------
         #ifdef MPICF

            // reduce total num atoms in cell
            MPI_Allreduce(MPI_IN_PLACE, &total_num_atoms_in_cell[0], int(num_cells), MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

            // cast to int for MPI
            int bufsize = num_sublattices;
            // reduce all cell totals onto all processors
            for( int cell = 0; cell < num_cells; cell++ ){
               MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_resistance          [cell][0], bufsize, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
               MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_spin_resistance     [cell][0], bufsize, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
               MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_relaxation_torque_rj[cell][0], bufsize, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
               MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_precession_torque_pj[cell][0], bufsize, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
               MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_isaturation         [cell][0], bufsize, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
               MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_alpha               [cell][0], bufsize, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
               MPI_Allreduce(MPI_IN_PLACE, &cell_sl_num_atoms_in_cell                 [cell][0], bufsize, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
               MPI_Allreduce(MPI_IN_PLACE, &cell_sl_num_magnetic_atoms                [cell][0], bufsize, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
               //MPI_Allreduce(MPI_IN_PLACE, &cell_sl_resistivity_sq                    [cell][0], bufsize, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
               //MPI_Allreduce(MPI_IN_PLACE, &cell_sl_spin_resistivity_sq               [cell][0], bufsize, MPI_DOUBLE,   MPI_SUM, MPI_COMM_WORLD);
            }
         #endif

         //-----------------------------------------------------------------------------------------------
         // Calculate average cell parameters
         //-----------------------------------------------------------------------------------------------
         // The idea is that each sublattice takes up a fractional cross-sectional area of each cell. The
         // total resistance follows the total resistance formula 1/R_T = 1/R_1 + 1/R_2 + ...
         //
         // The fractional cross sectional area of each sublattice is:
         //
         //     frac_xsa = num_atoms_in_sl_in_cell / total_num_atoms_in_cell
         //
         // The resistance of each sublattice in the cell is then given by
         //
         //             R_i = rho_eff L / (frac_xsa * A)
         //
         // where rho_eff is the mean resistivity of the constiuent materials, L is the length of the cell,
         // A is the cross-sectional area of the cell, and i = 1,2,... is the sublattice index.
         //
         // The total resistance then follows the parallel resistor formula:
         //
         //                  1/R_tot = 1/R_1 + 1/R_2 + ... 1/R_n
         //
         // calculating the total resistance for sublattices in each cell, then serial resistor formula
         // for all cells in a stack
         //


         // cell size parameters for resitivity to resistance calculation
         const double iA = 1.0 / (cell_size[stack_x] * cell_size[stack_y]); // Angstroms^-2
         const double L = cell_size[stack_z];

         // loop over all xy-cells (stacks)
         for(unsigned int i = 0; i < cells3D.size(); i++){
            for(unsigned int j = 0; j < cells3D[i].size(); j++){

               // loop over all cells in stack
               for(unsigned int k = 0; k < cells3D[i][j].size(); k++){

                  // determine cell ID
                  const uint64_t cell = cells3D[i][j][k].id;

                  // determine total number of atoms on all processors
                  const uint64_t num_atoms_in_cell = total_num_atoms_in_cell[cell];

                  // loop over all sublattices
                  for( int sl = 0; sl < num_sublattices; sl++ ){

                     // determine total number of sublattice atoms on all processors
                     const double num_sl_magnetic_atoms_in_cell = cell_sl_num_magnetic_atoms[cell][sl];
                     const uint64_t num_sl_atoms_in_cell = cell_sl_num_atoms_in_cell[cell][sl];

                     //------------------------------------------------------------
                     // if sublattice cell is empty space assume uniform padding resistance
                     //------------------------------------------------------------
                     if(num_sl_atoms_in_cell == 0){
                        // special case for sublattices with no atoms - assume zero additional resistance
                        // if the cell resistance is zero this will cause a "short" if in the entire stack
                        // so this is checked later when calculating the total resistance
                        st::internal::cell_sl_resistance[cell][sl]      = 0.0;
                        st::internal::cell_sl_spin_resistance[cell][sl] = 0.0;
                        st::internal::cell_sl_isaturation[cell][sl] = 0.0; // set inverse saturation to 0.0
                        st::internal::sl_magnetic[cell][sl] = false; // set sublattice in cell as non-magnetic
                     }
                     //-------------------------------------------------------------------------
                     // otherwise add contributions from atoms to calculate average resistivity
                     //-------------------------------------------------------------------------
                     else{
                        //------------------------------------------------------------
                        // check for empty cells or cells with only non-magnetic atoms (keep) and if so treat as non-magnetic
                        //------------------------------------------------------------
                        if(num_sl_magnetic_atoms_in_cell < 0.1){
                           // normalise resistance to number of atoms in sublattice
                           st::internal::cell_sl_resistance[cell][sl] = st::internal::cell_sl_resistance[cell][sl] / double(num_sl_atoms_in_cell);
                           st::internal::cell_sl_spin_resistance[cell][sl] = st::internal::cell_sl_spin_resistance[cell][sl] / double(num_sl_atoms_in_cell);
                           st::internal::sl_magnetic[cell][sl] = false; // no magnetic atoms -> non-magnetic cell
                           st::internal::cell_sl_isaturation[cell][sl] = 1.0; // assume 1 so inverse is still 1 (any value is fine but needs to be > 0)
                        }
                        //-------------------------------------------------------------------------
                        // finally consider magnetic atoms in each sublattice in each cell
                        //-------------------------------------------------------------------------
                        else{
                           // calculate mean spin torque prefactor parameters for magnetic cells
                           const double count = num_sl_atoms_in_cell; // number of magnetic atoms in sublattice in cell
                           st::internal::cell_sl_resistance[cell][sl] = st::internal::cell_sl_resistance[cell][sl] / count;
                           st::internal::cell_sl_spin_resistance[cell][sl] = st::internal::cell_sl_spin_resistance[cell][sl] / count;
                           st::internal::cell_sl_relaxation_torque_rj[cell][sl] = st::internal::cell_sl_relaxation_torque_rj[cell][sl] / count;
                           st::internal::cell_sl_precession_torque_pj[cell][sl] = st::internal::cell_sl_precession_torque_pj[cell][sl] / count;
                           st::internal::cell_sl_alpha[cell][sl]                = st::internal::cell_sl_alpha[cell][sl] / count;
                        }

                        // set inverse total moment
                        st::internal::cell_sl_isaturation[cell][sl] = 1.0/st::internal::cell_sl_isaturation[cell][sl]; // (1/mu_B)

                        // set inverse saturation of non-magnetic cells to zero
                        if(num_sl_magnetic_atoms_in_cell < 0.1){
                            st::internal::cell_sl_isaturation[cell][sl] = 0.0;
                        }

                        //------------------------------------------------------------
                        // calculate resistance values for each sublattice
                        //------------------------------------------------------------

                        // work out fractional area for each sublattice
                        const double fractional_area = double(num_sl_atoms_in_cell) / double(num_atoms_in_cell);

                        // get average resistivity for each sublattice in this cell
                        const double sublattice_cell_resistivity = st::internal::cell_sl_resistance[cell][sl];
                        const double sublattice_cell_spin_resistivity = st::internal::cell_sl_spin_resistance[cell][sl];

                        // set cell resistance R = rho * L / (f * A)
                        st::internal::cell_sl_resistance[cell][sl]      = sublattice_cell_resistivity * L * iA / fractional_area;
                        st::internal::cell_sl_spin_resistance[cell][sl] = sublattice_cell_spin_resistivity * L * iA / fractional_area;

                     } // end of if statement for cells containing atoms

                  } // end of sublattice loop
               } // end of cell loop
            } // end of stack y loop
         } // end of stack x loop

         //------------------------------------------------------------------------
         // Calculate cell data
         //------------------------------------------------------------------------
         st::internal::cell_position.resize(3*st::internal::total_num_cells);

         // convert to system coordinate system
         // translation from system coordinates to cell coordinates
         const int stack[3] = { stack_x, stack_y, stack_z };
         // example 1:             2        0        1
         // example 2:             0        1        2
         // example 3:             1        2        0

         // translation from cell coordinates to system coordinates
         // example 1: istack[3]   1        2        0
         // example 2: istack[3]   0        1        2
         // example 3: istack[3]   2        0        1

         // calculate inverse stack relations
         int istack[3];
         for(int i=0; i<3; i++){
            if(stack[i] == 0) istack[0] = i;
            if(stack[i] == 1) istack[1] = i;
            if(stack[i] == 2) istack[2] = i;
         }

         //std::cout << " stack: " << stack[0] << "\t" << stack[1] << "\t" << stack[2] << std::endl;
         //std::cout << "istack: " << istack[0] << "\t" << istack[1] << "\t" << istack[2] << std::endl;

         // loop over all xy-cells (stacks)
         for(unsigned int i = 0; i < cells3D.size(); i++){
            for(unsigned int j = 0; j < cells3D[i].size(); j++){

               // loop over all cells in stack
               for(unsigned int k = 0; k < cells3D[i][j].size(); k++){

                  // determine cell ID
                  const uint64_t cell = cells3D[i][j][k].id;

                  // determine cell counts in each direction
                  unsigned int xyz[3] = { i, j, k};

                  // calculate cell coordinates
                  const double x = double(xyz[istack[0]]) * cell_size[0];
                  const double y = double(xyz[istack[1]]) * cell_size[1];
                  const double z = double(xyz[istack[2]]) * cell_size[2];

                  st::internal::cell_position[3*cell+0] = x;
                  st::internal::cell_position[3*cell+1] = y;
                  st::internal::cell_position[3*cell+2] = z;

               }
            }
         }

         //------------------------------------------------------------------------
         // resize cell vector data arrays and set to zero
         //------------------------------------------------------------------------
         st::internal::cell_resistance.resize(num_cells,0.0); // also stores accumulated resistance before averaging

         st::internal::cell_sl_magnetization_x.resize(num_cells);
         st::internal::cell_sl_magnetization_y.resize(num_cells);
         st::internal::cell_sl_magnetization_z.resize(num_cells);
         st::internal::cell_sl_spin_torque_fields_x.resize(num_cells);
         st::internal::cell_sl_spin_torque_fields_y.resize(num_cells);
         st::internal::cell_sl_spin_torque_fields_z.resize(num_cells);
         for( int cell = 0; cell < num_cells; cell++ ){
            st::internal::cell_sl_magnetization_x[cell].resize(num_sublattices, 0.0);
            st::internal::cell_sl_magnetization_y[cell].resize(num_sublattices, 0.0);
            st::internal::cell_sl_magnetization_z[cell].resize(num_sublattices, 0.0);
            st::internal::cell_sl_spin_torque_fields_x[cell].resize(num_sublattices, 0.0);
            st::internal::cell_sl_spin_torque_fields_y[cell].resize(num_sublattices, 0.0);
            st::internal::cell_sl_spin_torque_fields_z[cell].resize(num_sublattices, 0.0);
         }

         return;

      }

   } // end of internal namespace

} // end of spin_transport namespace
