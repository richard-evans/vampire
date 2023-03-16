//-------------------------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the Free BSD licence (see licence
//   file for details).
//
//   (c) Jack B Collings, Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//-------------------------------------------------------------------------------------------------
//

/**************************************************************************************************
* @file calculate_A.cpp
* 
* @brief Contains function to calulate the exchange stiffness tensor ( A_{ i j } ) for a microcell.
*
* @details The isotropic exchange stiffness is
*
*           A_{ i j } = (1 / 4) * sum_{ i != j } [ J_{ i j } * ( x_j - x_i )^2 ] / V
* 
* where V is the micromagnetic cell, microcell, volume. Index i denotes summation over all atoms in
* the microcell, and j denotes summation over all nearest neighbours to atoms in the microcell. The
* value x_i - x_j is the distance between an atom and its neighbour in the x-direction (for
* isotropic exchange summing over x distances is equivalent to any other direction).
*
* We take some components of the exchange field into A, thus A is technically not what is returned;
* a field coefficient is retunred instead. The field coefficients, C, are given by
*
*                 C = [ 2 * A / (M * delta^2)] = [ 2 * A * V / ( m * delta^2 ) ]
*                          C = [ A * V * 4 / ( 2 * m * delta^2 ) ]
*
* where delta is the inter-microcell distance equal to microcell length in the axial direction
* considered. M is the magnetization, which is m / V, where m is the magnetic moment of the
* microcell.
**************************************************************************************************/

// Vampire headers
#include "micromagnetic.hpp"
#include "atoms.hpp"
#include "internal.hpp"
#include "material.hpp"
#include "create.hpp"

// Micromagnetic module headers
#include "errors.hpp"
#include "vio.hpp"

// C Libraries
#include <stdlib.h>

// C++ Libraries
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>

// Define micromagnetic namespace
namespace micromagnetic
{

// Define internal namespace
namespace internal
{
   
/** @brief The function that calculates the field coefficients.   */
std::vector < double > calculate_a (
   int num_atoms,                                  /** @param num_atoms total number of atoms.  */
   int num_cells,                                  /** @param num_cells number of microcells.   */
   int num_local_cells,                            /** @param num_local_cells number of microcells local to processor   */
   std::vector < int > cell_array,                 // Microcells atoms are in
   std::vector < int > neighbour_list_array,       // Nearest neighbours for each atom
   std::vector < int > neighbour_list_start_index, // Start index for atom's neighbours
   std::vector < int > neighbour_list_end_index,   // End index for atom's neigbours
   const std::vector < int > type_array,           // Material for each atom
   std::vector < mp::materials_t > material,       // Atom material parameters
   std::vector < double > volume_array,            // Volume of each microcell
   std::vector < double > x_coord_array,           // x-positions of atoms in the system
   std::vector < double > y_coord_array,           // y-positions of atoms in the system
   std::vector < double > z_coord_array,           // z-positions of atoms in the system
   std::vector < double > pos_array,               // Positions of microcells in the system
   double macro_x,                                 // x-length of microcells
   double macro_y,                                 // y-length of microcells
   double macro_z,                                 // z-length of microcells
   double num_atoms_in_unit_cell,                  // Number of atoms per unit cell -- remove
   std::vector < int > local_cell_array            /** @param local_cell_array cell array local to each processor */
)
{

   // Stores the exchange constants in a 2D array (they will be multiplied by microcell
   // volume and 4 hence the V and 4 in variable name)
   std::vector< std::vector< double > > a2dV4;

   // A_{ i j } initially set to zero, and the array length should be the microcell number
   a2dV4.resize( num_cells, std::vector < double > ( 3, 0.0 ) );

   bool check_microcell_coupling = true;

   // For MPI version, only add local atoms
   #ifdef MPICF
      unsigned int num_local_atoms =  vmpi::num_core_atoms + vmpi::num_bdry_atoms;
   #else
      unsigned int num_local_atoms = atoms::num_atoms;
   #endif

   // This vector will hold all the field coefficients
   std::vector < double > a ( 0, 0.0 );

   // Determines type of exchange stiffness to be calculated
   int exchange_type = 0;
   std::cout << "neighbour_list_array_size:\t" << neighbour_list_array.size() << std::endl;
   //----------------------------------------------------------------------------------------
   // Undertakes the calculation of the exchange stiffness based on type selected -----------
   //----------------------------------------------------------------------------------------
   switch( exchange_type )
   {

      // Atomistically generated, diagonal, inhomogeneous exchange stiffness tensor
      case 0:
      {

         //----------------------------------------------------------------------------------
         // Identify exchange coupled microcells to validate microcell neighbour list -------
         //----------------------------------------------------------------------------------

         // A map to connect each microcell to a vector of their neighbour microcells
         std::map < int, std::vector < int > > exchange_coupled_microcell_map;

         // Loop over all atoms in the system
         for ( unsigned int atom = 0; atom < num_local_atoms; ++atom )
         {

            // Store cell of atom
            const unsigned int acell = cell_array[ atom ];

            // Get start and end index of atom's neighbours
            const unsigned int start   = atoms::neighbour_list_start_index[ atom ];
            const unsigned int end     = atoms::neighbour_list_end_index[ atom ];

            std::cout << "neighbour_list_start_index_array_size:\t" << atoms::neighbour_list_start_index.size() << std::endl;
            std::cout << "neighbour_list_end_index_array_size:\t" << atoms::neighbour_list_end_index.size() << std::endl;

            std::cout << "Start:\t" << start << "\tEnd:\t" << end << std::endl;
            // Loop through all nearest neighbours
            for ( unsigned int natom = start; natom <= end; ++natom ) /** @todo is this correct? */
            {

               std::cout << "natom:\t" << natom << std::endl;
               // Get microcell of neighbour atom
               const unsigned int ncell   = cell_array[ atoms::neighbour_list_array[ natom ] ];

               // If neighbour atom is in same microcell, ignore
               if ( ncell == acell ) continue;

               // Reject neighbour atom in microcell already included as a neighbour
               else if ( std::find( exchange_coupled_microcell_map[ acell ].begin() , exchange_coupled_microcell_map[ acell ].end(), ncell ) != exchange_coupled_microcell_map[ acell ].end() ) continue;
               
               // Not a duplicate, so add to the vector of neighbours
               else exchange_coupled_microcell_map[ acell ].push_back( ncell );

            }

         }

         // Print the exchange coupled microcell map to check for errors
         /*
         for ( unsigned int cell = 0; cell < num_cells; ++cell )
         {

            // Print the microcell
            std::cout << "Cell:\t" << cell << "\t Neighbours:\t";

            // Loop through this cell's vector of neighbours
            for ( unsigned int ncell = 0; ncell < exchange_coupled_microcell_map[ cell ].size(); ++ncell )
            {

               // Print the coupled neighbour cell index
               std::cout << exchange_coupled_microcell_map[ cell ][ ncell ] << "\t";

            }

            std::cout << std::endl;

         }
         */

         //----------------------------------------------------------------------------------
         // The diagonal exchange stiffness tensor is calculated for each cell --------------
         //----------------------------------------------------------------------------------

         // Loop over all atoms in the system
         for ( unsigned int atom = 0; atom < num_local_atoms; ++atom )
         {

            // Save the material type of the atom
            const unsigned int mat  = type_array[ atom ];

            // Save the cell atom belongs to
            const unsigned int cell = cell_array[ atom ];
            
            // Find the start and end index for the atom's neighbours
            const unsigned int start   = atoms::neighbour_list_start_index[ atom ];
            const unsigned int end     = atoms::neighbour_list_end_index[ atom ];


            // Loop over all nearest neighbours
            for( unsigned int nn = start; nn < end; ++nn )
            {

               // Calcualte the atom id of the nn atom
               const unsigned int natom   = atoms::neighbour_list_array[ nn ];

               // Calculate component square distances between the two atoms
               const double dx   = x_coord_array[ atom ] - x_coord_array[ natom ];
               const double dx2  = dx * dx;
               const double dy   = y_coord_array[ atom ] - y_coord_array[ natom ];
               const double dy2  = dy * dy;
               const double dz   = z_coord_array[ atom ] - z_coord_array[ natom ];
               const double dz2  = dz * dz;

               // J_{ i j } is stored as J_{ i j } / mu_s so have to multiply by mu_s
               const double Jij  = atoms::i_exchange_list[ atoms::neighbour_interaction_type_array[ nn ] ].Jij * mp::material[ mat ].mu_s_SI;
               
               // ( Exchange stiffness ) * ( Volume ) * 4 is
               // ( distance squared along axis ) * ( J_{ i j } ) summed for all interactions
               // in the microcell
               a2dV4[ cell ][ 0 ] += Jij * dx2;
               a2dV4[ cell ][ 1 ] += Jij * dy2;
               a2dV4[ cell ][ 2 ] += Jij * dz2;
               
            }

         }

         // Check for homogeneous isotropic override to simulate standard micromagnetics
         if ( homogeneous_isotropic_exchange )
         {
            
            // Loop over number of cells to override value of exchange stiffness tensor
            for ( unsigned int cell = 0; cell < a2dV4.size(); ++cell )
            {
               
               // Set the values to the override homogneneous exchange values given
               a2dV4[ cell ][ 0 ] = homogeneous_isotropic_exchange_value * 4.0 * volume_array[cell];
               a2dV4[ cell ][ 1 ] = homogeneous_isotropic_exchange_value * 4.0 * volume_array[cell];
               a2dV4[ cell ][ 2 ] = homogeneous_isotropic_exchange_value * 4.0 * volume_array[cell];
            
            }
         
         }

         // Print out exchange stiffness values for each cell to check for errors
         /*
         std::cout << "Exchange Stiffness Values Given for Each Cell" << std::endl;
         for ( unsigned int cell = 0; cell < num_cells; ++cell )
         {
            
            const double cell_volume = volume_array[ cell ];
            
            std::cout << "Cell:\t" << cell << std::endl;
            std::cout << "A_xx\t" << a2dV4[ cell ][ 0 ] / ( cell_volume * 4.0 ) << std::endl;
            std::cout << "A_yy\t" << a2dV4[ cell ][ 1 ] / ( cell_volume * 4.0 ) << std::endl;
            std::cout << "A_zz\t" << a2dV4[ cell ][ 2 ] / ( cell_volume * 4.0 ) << std::endl;

         }
         std::cout << "----------------------------------------------" << std::endl;
         */

         //----------------------------------------------------------------------------------
         // Create microcell neighbour list for second-order approximation to second order mu
         // derivative with respect to space, and attribute correct exchange stiffness ------
         // tensor to each microcell --------------------------------------------------------
         //----------------------------------------------------------------------------------
         
         // Index for setting start and end indices of microcell neighbour list
         unsigned int array_index = 0;

         // Useful constants to determine if neighbour microcell is a nearest neighbour
         const double err_factor = 0.00001;
         const double errx_offset = err_factor * macro_x;
         const double erry_offset = err_factor * macro_y;
         const double errz_offset = err_factor * macro_z;

         // Useful constants to calculate the field coefficients
         const double halfomacrox2 = 0.5 / ( macro_x * macro_x );
         const double halfomacroy2 = 0.5 / ( macro_y * macro_y );
         const double halfomacroz2 = 0.5 / ( macro_z * macro_z );

         // Loop through all microcells in the system
         for ( unsigned int cell = 0; cell < num_cells; ++cell )
         {

            // Set neighbour list start index for this cell
            macro_neighbour_list_start_index[ cell ] = array_index;

            // Define useful constants for parsing cell position array
            const unsigned int cell0 = cell * 3;
            const unsigned int cell1 = cell0 + 1;
            const unsigned int cell2 = cell0 + 2;
            
            // Loop through all microcells and test to find the neighbouring microcells
            for ( unsigned int neighbour = 0; neighbour < num_cells; ++neighbour )
            {

               // Check microcell and neighbour microcell are not the same
               if ( cell != neighbour )
               {

                  // If neighbour is not exchange coupled, ignore
                  if ( std::find( exchange_coupled_microcell_map[ cell ].begin(), exchange_coupled_microcell_map[ cell ].end(), neighbour ) == exchange_coupled_microcell_map[ cell ].end() ) continue;

                  // Define useful constants for parsing pos_array
                  const unsigned int neighbour0 = neighbour * 3;
                  const unsigned int neighbour1 = neighbour0 + 1;
                  const unsigned int neighbour2 = neighbour0 + 2;

                  // Obtain inter-microcell distances
                  const double abs_deltax = abs( pos_array[ cell0 ] - pos_array[ neighbour0 ] );
                  const double abs_deltay = abs( pos_array[ cell1 ] - pos_array[ neighbour1 ] );
                  const double abs_deltaz = abs( pos_array[ cell2 ] - pos_array[ neighbour2 ] );
                  
                  // Calculate the field coefficients factor to AV4
                  const double m_s = ms[ cell ];
                  const double factorx = halfomacrox2 / m_s;
                  const double factory = halfomacroy2 / m_s;
                  const double factorz = halfomacroz2 / m_s;

                  // Only consider neighbour microcells in nearest neighbour range
                  if ( ( abs_deltax - errx_offset < macro_x ) && ( abs_deltay - erry_offset < macro_y ) && ( abs_deltaz - errz_offset < macro_z ) )
                  {
                     
                     // Find linear x, y neighbours, i.e., must be aligned in z direction
                     if ( abs_deltaz - errz_offset < 0 )
                     {

                        // The [ neighbour != cell ] condition used earlier ensures x, y, z
                        // positions of neighbour microcell and microcell cannot be identical
                        
                        // It is an x linear neighbour if aligned in y direction
                        if ( abs_deltay - erry_offset < 0 )
                        {
                           
                           // Add to neighbour list array
                           macro_neighbour_list_array.push_back( neighbour );
                        
                           // Push field coefficient into array
                           a.push_back( a2dV4[ cell ][ 0 ] * factorx );
                        
                           // Increment array index for next neighbour position
                           ++array_index;
                        
                           // Neighbour identified, now move to the next candidate
                           continue;

                        }

                        // It is a y linear neighbour if aligned in the x direction
                        else if ( abs_deltax - errx_offset < 0 )
                        {

                           // Add to neighbour list array
                           macro_neighbour_list_array.push_back( neighbour );

                           // Push field coefficient into array
                           a.push_back( a2dV4[ cell ][ 1 ] * factory );

                           // Increment array index for next neighbour position
                           ++array_index;

                           // Neighbour identified, now move on to the next candidate
                           continue;

                        }

                     }

                     // Only other case of interest is if the neighbour microcell is a
                     // linear z neighbour
                     else if ( ( abs_deltax - errx_offset < 0 ) && ( abs_deltay - erry_offset < 0 ) )
                     {

                        // Add to neighbour list array
                        macro_neighbour_list_array.push_back( neighbour );

                        // Push field coefficient into array
                        a.push_back( a2dV4[ cell ][ 2 ] * factorz );

                        // Increment array index for next neighbour position
                        ++array_index;

                        // Neighbour identified, now move on to the next candidate
                        continue;

                     }

                     // Neighbour cannot be a linear neighbour so move on to next
                     // candidate
                     else continue;

                  }

                  // Ignore neighbour microcells outside nearest neighbour range
                  else continue;

               }

            }

            // Update neighbour list end index
            macro_neighbour_list_end_index[ cell ] = array_index - 1;

         }

         break;

      }

      case 1: // Vector
      {

         // Set terminal colour to red for the error statement
         terminaltextcolor( RED );

         // Print out the error statement
         std::cerr << "Error! Vectoral exchange calculation not yet implemented in micromagnetic mode" << std::endl;
         
         // Set terminal colour back to white
         terminaltextcolor( WHITE );
         
         // Add error statement to the log
         zlog << zTs() << "Error! Vectoral exchange calculation not yet implemented in micromagnetic mode" << std::endl;
         
         // Exit as error
         err::vexit();
         
         break;

      }

      case 2: // Tensor
      {

         // Set terminal colour to red for the error statement
         terminaltextcolor( RED );

         // Print out the error statement
         std::cerr << "Error! Tensor exchange calculation not yet implemented in micromagnetic mode" << std::endl;
         
         // Set terminal colour back to white
         terminaltextcolor( WHITE );

         // Add error statement to the log
         zlog << zTs() << "Error! Tensor exchange calculation not yet implemented in micromagnetic mode" << std::endl;
         
         // Exit as error
         err::vexit();
      
      break;
      
      }
   
   }

   // Print out field coefficients to check for errors
   // Loop through all cells
   std::cout << "a.size():\t" << a.size() << std::endl;
   // std::cout << "a[0]:\t" << a[0] << std::endl;
   /*
   for( int i = 0; i < num_cells; ++i )
   {
      
      for( int j = macro_neighbour_list_start_index[ i ]; j < macro_neighbour_list_end_index[ i ]; ++j )
      {
         
         std::cout << "Cell:\t" << i << "\tNeighbour cell:\t" <<  macro_neighbour_list_array[ j ] << std::endl;
         std::cout << "Field coefficient value:\t" << a[ macro_neighbour_list_array[ j ] ] << std::endl;
      
      }
   
   }
   */
   if ( a.size() == 0 ) a.push_back(0.0);
   // Return the 1D vector of the cellular exchange field coefficients
   return a;

}

} // End internal namespace

} // End micromagnetic namespace
