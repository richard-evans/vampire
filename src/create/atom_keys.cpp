//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) S Jenkins 2026. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

// Vampire headers
#include "create.hpp"
#include "errors.hpp"
#include "vmpi.hpp"

// create module headers
#include "internal.hpp"

namespace create{

   namespace internal{

      // the key is a 64 bit number made of 4 16 bit chunks: the unit cell's x, y, and z coordinates followed by the atom's index in that unit cell. This gives each atom a unique key.
      const uint64_t key_limit = uint64_t(1) << 16;

      //--------------------------------------------------------------------------------------------------
      // Function to build an atom's key value from its number in the unit cell
      // and the unit cell's x,y, and z coordinates, accounting for periodic boundary conditions.
      // For the key to be valid, none of these values may exceed the 16 bit integer limit (see key_limits_ok())
      //--------------------------------------------------------------------------------------------------
      uint64_t atom_key(const cs::catom_t& atom){
         const int64_t sc[3] = {atom.scx, atom.scy, atom.scz};
         uint64_t key = atom.uc_id;
         for(int i = 0; i < 3; i++){
            const int64_t n = int64_t(cs::total_num_unit_cells[i]);
            key |= uint64_t(((sc[i] % n) + n) % n) << (16 * (i + 1));
         }
         return key;
      }

      //------------------------------------------------------------------------
      // Function to check that every atom of the system fits into a key
      //------------------------------------------------------------------------
      bool key_limits_ok(const std::vector<cs::catom_t>& catom_array){
         for(int i = 0; i < 3; i++) if(cs::total_num_unit_cells[i] >= key_limit) return false;
         for(size_t atom = 0; atom < catom_array.size(); atom++) if(catom_array[atom].uc_id >= key_limit) return false;
         return true;
      }

      //------------------------------------------------------------------------
      // Function to return a uniform random number in [0,1) determined only by
      // a seed, the key of an atom and the number of the draw for that atom
      // (splitmix64 finaliser), so that random numbers are decompositionally independent
      //------------------------------------------------------------------------
      double keyed_uniform(const uint64_t seed, const uint64_t key, const uint64_t draw){
         uint64_t z = seed * 0x9e3779b97f4a7c15ULL + key;
         z = (z ^ (draw + 0x632be59bd9b4e5fULL)) * 0xbf58476d1ce4e5b9ULL;
         for(int i = 0; i < 2; i++){
            z += 0x9e3779b97f4a7c15ULL;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            z ^= z >> 31;
         }
         return double(z >> 11) * (1.0 / 9007199254740992.0); // 53 random bits
      }

      //------------------------------------------------------------------------
      // Function to return a uniform random number for a random choice of
      // creation: the next number of the create module's generator, or with
      // create:decomposition-independent-random draw number <draw> of the atom
      //------------------------------------------------------------------------
      double draw_uniform(const int seed, const cs::catom_t& atom, const uint64_t draw){
         if(decomposition_independent_random) return keyed_uniform(uint64_t(seed), atom_key(atom), draw);
         return grnd();
      }

      //------------------------------------------------------------------------
      // Function to check that the keyed random numbers can be used
      // (create:decomposition-independent-random)
      //------------------------------------------------------------------------
      void check_keyed_random(const std::vector<cs::catom_t>& catom_array){
         if(!key_limits_ok(catom_array)){
            err::fatal("Error: create:decomposition-independent-random supports at most " + std::to_string(key_limit) +
                       " unit cells along each axis and atoms per unit cell");
         }
      }

   } // end of internal namespace

   //---------------------------------------------------------------------------
   // Function to store the key of every atom (core, boundary and halo). The keys
   // are valid if every atom fits into a key and the local atoms of this
   // process have distinct keys (true for structures built from a lattice)
   //---------------------------------------------------------------------------
   void set_atom_keys(const std::vector<cs::catom_t>& catom_array){

      atom_keys.clear();
      atom_keys_valid = internal::key_limits_ok(catom_array);
      if(!atom_keys_valid) return;

      atom_keys.resize(catom_array.size());
      for(size_t atom = 0; atom < catom_array.size(); atom++) atom_keys[atom] = internal::atom_key(catom_array[atom]);

      #ifdef MPICF
         const size_t num_local_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
      #else
         const size_t num_local_atoms = catom_array.size();
      #endif
      std::vector<uint64_t> local(atom_keys.begin(), atom_keys.begin() + num_local_atoms);
      std::sort(local.begin(), local.end());
      if(std::adjacent_find(local.begin(), local.end()) != local.end()){
         atom_keys.clear();
         atom_keys_valid = false;
      }

   }

} // end of create namespace
