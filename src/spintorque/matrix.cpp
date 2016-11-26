//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>

// Vampire headers
#include "spintorque.hpp"

// Spin Torque headers
#include "internal.hpp"

// Contains functions for matrix operations
namespace st{
   namespace internal{

      //------------------------------------------------
      // Function to set transformation matrix m -> m'
      //------------------------------------------------
      void set_inverse_transformation_matrix(const st::internal::three_vector_t& reference_vector, st::internal::matrix_t& itm){

         // load vector into temporary variables for readability
         const double x = reference_vector.x;
         const double y = reference_vector.y;
         const double z = reference_vector.z;

         // Transform vectors in local coordinate system
         const double D1 = sqrt(y*y+z*z);
         const double D2 = sqrt(x*x+y*y+z*z);
         const double iD1 = 1.0/D1;
         const double iD2 = 1.0/D2;

         // define result matrix
         st::internal::matrix_t tm; // transformation matrix

         tm.xx = D1*iD2;       // D1/D2
         tm.xy = -x*y*iD1*iD2; // -x*y/(D1*D2)
         tm.xz = -x*z*iD1*iD2; // -x*z/(D1*D2)

         tm.yx = 0.0;
         tm.yy = z*iD1;  // z/D1
         tm.yz = -y*iD1; // -y/D1;

         tm.zx = x*iD2;  // x/D2
         tm.zy = y*iD2;  // y/D2
         tm.zz = z*iD2;  // z/D2

         // Set determinants for inverse matrix
         const double det1 = tm.xx*tm.yy*tm.zz + tm.xy*tm.yz*tm.zx + tm.xz*tm.yx*tm.zy;
         const double det2 = tm.yz*tm.zy*tm.xx + tm.xz*tm.zx*tm.yy + tm.xy*tm.yx*tm.zz;

         const double det_T = det1-det2;
         const double inv_det = 1.0/det_T;

         // Calculate inverse transformation matrix
         itm.xx = inv_det*(tm.yy*tm.zz - tm.yz*tm.zy);
         itm.xy = inv_det*(tm.xz*tm.zy - tm.xy*tm.zz);
         itm.xz = inv_det*(tm.xy*tm.yz - tm.xz*tm.yy);

         itm.yx = inv_det*(tm.yz*tm.zx - tm.yx*tm.zz);
         itm.yy = inv_det*(tm.xx*tm.zz - tm.xz*tm.zx);
         itm.yz = inv_det*(tm.xz*tm.yx - tm.xx*tm.yz);

         itm.zx = inv_det*(tm.yx*tm.zy - tm.yy*tm.zx);
         itm.zy = inv_det*(tm.xy*tm.zx - tm.xx*tm.zy);
         itm.zz = inv_det*(tm.xx*tm.yy - tm.xy*tm.yx);

         return;

      }

      //--------------------------------------------------------------------------------------------------
      // Function to transform reference vector (rv) using transformation matrix (tm) moving rv -> result
      //--------------------------------------------------------------------------------------------------
      st::internal::three_vector_t transform_vector(const st::internal::three_vector_t& rv, const st::internal::matrix_t& tm){

         // Declare result vector
         st::internal::three_vector_t result(0.0,0.0,0.0);

         // Calculate transformation
         result.x = tm.xx*rv.x + tm.xy*rv.y + tm.xz*rv.z;
         result.y = tm.yx*rv.x + tm.yy*rv.y + tm.yz*rv.z;
         result.z = tm.zx*rv.x + tm.zy*rv.y + tm.zz*rv.z;

         return result;

      }

      //------------------------------------------------------------------------------------------------
      // Black magic code for Gaussian elimination
      //------------------------------------------------------------------------------------------------

      st::internal::three_vector_t gaussian_elimination(st::internal::matrix_t& M, st::internal::three_vector_t& V){

         // initialise temporary arrays
         double a_array[3][3] = {{M.xx,M.xy,M.xz},{M.yx,M.yy,M.yz},{M.zx,M.zy,M.zz}};
         double b_array[3]    = {V.x,V.y,V.z};

         // initialise result array
         double x[3]={0.0,0.0,0.0};

         // temporary variables
         double ratio,temp;
         int n=3;
         double tem1, tem2, tem3, tem4;
         int p;

         // Rearranging matrix(pivot)
         for(int i=0;i<n;i++){  //column i
            tem1=fabs(a_array[i][i]);
            p=i;
            for(int j=i+1;j<n;j++){
               tem1 = fabs(tem1);
               tem2=fabs(a_array[j][i]);
               if(tem2>tem1){
                  p=j;
                  tem1=a_array[j][i];
               }
            }
            //row exchange in both the matrix
            for(int j=0;j<n;j++){
               tem3=a_array[i][j];
               a_array[i][j]=a_array[p][j];
               a_array[p][j]=tem3;
            }
            tem4=b_array[i];
            b_array[i]=b_array[p];
            b_array[p]=tem4;
         }

         // Gaussian Elimination
         for(int i=0; i<(n-1); i++){
            for(int j=i+1;j<n;j++){
               ratio=a_array[j][i]/a_array[i][i];     //diagonal value should not be zero
               for(int count=i; count<n; count++){
                  a_array[j][count]-=ratio*a_array[i][count];
               }
               b_array[j]-= ratio*b_array[i];
            }
         }

         x[n-1]=b_array[n-1]/a_array[n-1][n-1];
         for(int i=(n-2);i>=0;i--){
            temp=b_array[i];
            for(int j=(i+1);j<n;j++){
               temp-= a_array[i][j]*x[j];
            }
            x[i]=temp/a_array[i][i];
         }

         // copy to three vector
         st::internal::three_vector_t result(x[0],x[1],x[2]);

         return result;

      }

   } // end of internal namespace
} // end of spin torque namespace
