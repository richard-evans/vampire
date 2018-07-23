//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef CREATE_ATOMS_CLASS_H_
#define CREATE_ATOMS_CLASS_H_

//namespace create{
namespace cs{

//------------------------------------------------------------------
// Simple class to store atom properties in object oriented way
//------------------------------------------------------------------
class catom_t {

public:

   // Coordinates
   double x; // x-position of atom
   double y; // y-position of atom
   double z; // z-position of atom

   // Global atomic coordinates
   uint64_t uc_id;        // atom number of host unit cell
   int64_t scx;                   // supercell x coordinate of atom |
   int64_t scy;                   // supercell y coordinate of atom |
   int64_t scz;                   // supercell z coordinate of atom /

   // Flags
   bool include; // boolean to incude atom in structure (or not)
   bool boundary; // boolean to determine if atom interacts with MPI halo
   bool non_interacting_halo; // boolean to determine if atom is non-interacting halo

   // Integers
   int material;              // atom material belongs to
   int uc_category;           // atom category within unit cell
   int lh_category;           // atom height category within unit cell
   int grain;                 // grain id of atom
   //int supercell;       // supercell id of atom
   int mpi_type;              // mpi category of atom (core, boundary or halo)
   int mpi_cpuid;             // CPU id atom is located on
   int mpi_atom_number;       //
   int mpi_old_atom_number;   //
   unsigned int nn;           // number of neighbours
   unsigned int nbqn;         // number of biquadratic neighbours

   //----------------------------------
   // Class constructor
   //----------------------------------
   catom_t():
      x(0.0),
      y(0.0),
      z(0.0),
      uc_id(0),
      scx(0),
      scy(0),
      scz(0),
      include(false),
      boundary(false),
      non_interacting_halo(true),
      material(0),
      uc_category(0),
      lh_category(0),
      grain(0),
      //supercell(0),
      mpi_type(0),
      mpi_cpuid(0),
      mpi_atom_number(0),
      mpi_old_atom_number(0),

      nn(0),
      nbqn(0)
   {
      // Do nothing
      return;
   };

}; // end of class catom

} // end of namespace create

#endif /*CREATE_ATOMS_CLASS_H_*/
