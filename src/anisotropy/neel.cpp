//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: sw766@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //---------------------------------------------------------------------------------
      // Function to add second order uniaxial anisotropy along vector e
      //
      // Example 1:
      //                               o -- ø -- o
      //                           y        |
      //                           ^        o
      //                           |__ > x
      //
      //
      //       Energy (scalar) = -ks/2 sum_j ( S_i . e_ij )^2            (1)
      //                       = -ks/2 ( S_i sum_j e_ij )^2
      //                       = -ks/2 ( S_i S_i \sum_j e_ij e_ij )      (2)
      //
      //       Vectors e_ij = [-1,0,0],[1,0,0],[0,-1,0]
      //
      //       Energy from Eq. 1 = -ks/2 [ ( -Sx . - Sx )^2 + ( Sx . Sx )^2 + ( -Sy . - Sy )^2 ]
      //                         = -ks/2 [ 2 Sx^2 + Sy^2 ]
      //
      //       Energy from Eq. 2 = -ks/2 [ Sx^2 ( -1.-1 + 1.1 ) + Sy^2 (-1.-1) ]
      //                         = -ks/2 [ 2 Sx^2 + Sy^2 ]
      //
      //       Tensor T[ij] = sum_j e_ij[i] e_ij[j] = [  2   0   0  ]
      //                                              [  0   1   0  ]
      //                                              [  0   0   0  ]
      //
      //       Energy (tensor) = S . T . S
      //
      //                       = [ Sx Sy Sz ] [  2   0   0  ] [ Sx ]
      //                                      [  0   1   0  ] [ Sy ]
      //                                      [  0   0   0  ] [ Sz ]
      //
      //                       = Sx Txx Sx + Sx Txy Sy + Sx Txz Sz +
      //                         Sy Tyx Sx + Sy Tyy Sy + Sy Tyz Sz +
      //                         Sz Tzx Sx + Sz Tzy Sy + Sz Tzz Sz
      //
      //                       = Sx 2 Sx + Sy Sy
      //
      //---------------------------------------------------------------------------------
      void neel_anisotropy(const unsigned int num_atoms){

         //----------------------------------------------------------------------------------
         // Loop over all atoms and add neel tensor
         //----------------------------------------------------------------------------------
         for (int atom=0; atom < num_atoms; ++atom){

            //std::cout << atom << std::endl;

            for (int i = 0; i < 3; ++i){
               for (int j = 0; j < 3; ++j){
                  internal::second_order_tensor[ index(atom, i, j) ] += internal::neel_tensor[ index(atom, i, j) ];
               }

               // print tensor
               //std::cout << "[  " << internal::second_order_tensor[index(atom, i, 0)] << "  " <<
               //                      internal::second_order_tensor[index(atom, i, 1)] << "  " <<
               //                      internal::second_order_tensor[index(atom, i, 2)] << "  " << "]" << std::endl;
               
            }

         }

         // delete memory of neel tensor
         std::vector<double> blank(0);
         internal::neel_tensor.swap(blank);

         return;

      }




/*      void calculate_surface_anisotropy_fields(const int start_index,const int end_index){
      	///======================================================
      	/// 		Subroutine to calculate surface anisotropy fields
      	///
      	///			Version 1.0 Richard Evans 13/09/2011
      	///======================================================

      	// check calling of routine if error checking is activated
      	if(err::check==true){std::cout << "calculate_surface_anisotropy_fields has been called" << std::endl;}

      	for(int atom=start_index;atom<end_index;atom++){
      		// only calculate for surface atoms
      		if(atoms::surface_array[atom]==true){
      			const int imaterial=atoms::type_array[atom];
      			const double Ks=0.5*2.0*mp::material[imaterial].Ks; // note factor two here from differentiation
      			const double S[3]={atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};

      			for(int nn=atoms::nearest_neighbour_list_si[atom];nn<atoms::nearest_neighbour_list_ei[atom];nn++){
      				const double si_dot_eij=(S[0]*atoms::eijx[nn]+S[1]*atoms::eijy[nn]+S[2]*atoms::eijz[nn]);
      				atoms::x_total_spin_field_array[atom]-=Ks*si_dot_eij*atoms::eijx[nn];
      				atoms::y_total_spin_field_array[atom]-=Ks*si_dot_eij*atoms::eijy[nn];
      				atoms::z_total_spin_field_array[atom]-=Ks*si_dot_eij*atoms::eijz[nn];
      			}
      		}
      	}

      	return;
      }



      /// @brief Calculates the surface anisotropy energy for a single spin.
      ///
      /// @section License
      /// Use of this code, either in source or compiled form, is subject to license from the authors.
      /// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
      ///
      /// @section Information
      /// @author  Richard Evans, rfle500@york.ac.uk
      /// @version 1.0
      /// @date    07/02/2011
      ///
      /// @param[in] atom atom number
      /// @param[in] imaterial material of local atom
      /// @param[in] Sx x-spin of local atom
      /// @param[in] Sy y-spin of local atom
      /// @param[in] Sz z-spin of local atom
      /// @return exchange energy
      ///
      /// @internal
      ///	Created:		13/09/2011
      ///	Revision:	  ---
      ///=====================================================================================
      ///
      double spin_surface_anisotropy_energy(const int atom, const int imaterial, const double Sx, const double Sy, const double Sz){

      	double energy=0.0;

      	if(atoms::surface_array[atom]==true && sim::surface_anisotropy==true){
      		const double Ks=mp::material[imaterial].Ks*0.5;
      		for(int nn=atoms::nearest_neighbour_list_si[atom];nn<atoms::nearest_neighbour_list_ei[atom];nn++){
      			const double si_dot_eij=(Sx*atoms::eijx[nn]+Sy*atoms::eijy[nn]+Sz*atoms::eijz[nn]);
      			energy+=Ks*si_dot_eij*si_dot_eij;
      		}
      	}

      	return energy;
      }
*/

   } // end of internal namespace

} // end of anisotropy namespace
