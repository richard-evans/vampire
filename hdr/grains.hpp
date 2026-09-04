//-----------------------------------------------------------------------------
//
// This file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) / Free BSD licence (see licence file for details).
//
// (c) R F L Evans 2009-2026. All rights reserved.
//
// Email: richard.evans@york.ac.uk
//
//-----------------------------------------------------------------------------

#ifndef GRAINS_H_
#define GRAINS_H_

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers
#include "create.hpp" // cs::catom_t

//--------------------------------------------------------------------------------
// Namespace for variables and functions for the grains module: generation of
// Voronoi/Laguerre-tessellated granular structures (granular thin films, and
// nucleated grain substructures within a particle), every downstream step
// that operates on the resulting grain polygons (weight fitting, faceting,
// rounding, boundary spacing, atom assignment), and per-grain magnetic
// property bookkeeping once system creation (of any kind - a granular film,
// a particle array, a nanoparticle agglomerate, ...) has assigned every atom
// to a grain.
//--------------------------------------------------------------------------------
namespace grains{

   //-----------------------------------------------------------------------------
   // Per-grain magnetic-property bookkeeping, read/written directly by
   // system-creation code outside this module (particle arrays, the
   // nanoparticle agglomerate, random per-grain anisotropy): num_grains is
   // set by whichever system-creation path produced the grains, then
   // set_properties() (called once creation is complete) derives everything
   // else from atoms::grain_array.
   //-----------------------------------------------------------------------------
   extern int num_grains;
   extern bool random_anisotropy; // flag to control randomly oriented uniaxial anisotropy

   // Computes per-grain atom count, mean coordinates and summed saturation
   // magnetic moment from atoms::grain_array, after first compacting grain
   // numbers so no empty grain id remains. Writes grain-coordinates.txt on
   // the root process. Called once, after atom creation is complete.
   int set_properties();

   //-----------------------------------------------------------------------------
   // Function to initialize grains module
   //-----------------------------------------------------------------------------
   void initialize();

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for grains module
   //-----------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //-----------------------------------------------------------------------------
   // Function to process material parameters for grains module
   //-----------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

   //-----------------------------------------------------------------------------
   // Builds a granular thin film (create:granular-film / create:voronoi-film
   // system type): samples a grain-size distribution, places seeds by relaxed
   // disc packing, tessellates a Laguerre (or plain Voronoi) diagram, applies
   // optional faceting/rounding/boundary spacing, and assigns atoms to grains.
   //-----------------------------------------------------------------------------
   int voronoi_film(std::vector<cs::catom_t> &);

   //-----------------------------------------------------------------------------
   // Nucleates a grain substructure within an already-created particle
   // (create:voronoi-substructure).
   //-----------------------------------------------------------------------------
   void voronoi_substructure(std::vector<cs::catom_t> &);

} // end of namespace grains

#endif //GRAINS_H_
