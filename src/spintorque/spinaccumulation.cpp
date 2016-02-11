//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <fstream>
#include <iostream>

// Vampire headers
#include "spintorque.hpp"
#include "vmpi.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Funtion to calculate the spin accumulation and spin torque
      //-----------------------------------------------------------------------------
      void calculate_spin_accumulation(){

         // Zero all shared arrays (essential for parallelisation to work)
         std::fill (st::internal::spin_torque.begin(),st::internal::spin_torque.end(),0.0);

         // Declare resuable temporary variables
         st::internal::matrix_t itm; // inverse transformation matrix
         st::internal::matrix_t M; // general matrix
         st::internal::three_vector_t V(0.0,0.0,0.0); // general 3-vector

         // reference basis vectors
         const st::internal::three_vector_t bx(1.0,0.0,0.0);
         const st::internal::three_vector_t by(0.0,1.0,0.0);
         const st::internal::three_vector_t bz(0.0,0.0,1.0);

         // local basis vectors
         st::internal::three_vector_t b1(1.0,0.0,0.0);
         st::internal::three_vector_t b2(0.0,1.0,0.0);
         st::internal::three_vector_t b3(0.0,0.0,1.0);


       //  std::cout<< st::internal::je << "\t" << st::internal::current_direction << "\t" << micro_cell_thickness << "\t" << micro_cell_size <<   std::endl;

         // set local constants
         const double je = st::internal::je; // current (C/s)
         const double i_muB = 1.0/9.274e-24; // J/T
         const double i_e = 1.0/1.60217662e-19; // electronic charge (Coulombs)
         const double microcell_volume = (st::internal::micro_cell_size *
                                          st::internal::micro_cell_size *
                                          st::internal::micro_cell_thickness)*1.e-30; // m^3


         // loop over all 1D stacks (in parallel)
         for(int stack=0; stack <num_stacks; ++stack){
            // determine starting cell in stack
            const int idx = stack_index[stack];

            // set initial values


            st::internal::sa[3*idx+0] = 0.0;
            st::internal::sa[3*idx+1] = 0.0;
            st::internal::sa[3*idx+2] = 0.0; //10.e6;// t::internal::default_properties.sa_infinity;

            st::internal::j [3*idx+0] = st::internal::initial_beta*je*st::internal::initial_m[0];
            st::internal::j [3*idx+1] = st::internal::initial_beta*je*st::internal::initial_m[1];
            st::internal::j [3*idx+2] = st::internal::initial_beta*je*st::internal::initial_m[2];


      //  std::cout<< st::internal::initial_beta << "\t" << st::internal::j[0] << "\t" << st::internal::j[1]  << "\t" << st::internal::j[2]  << "\t" << std::endl;


            // loop over all cells in stack after first (idx+1)
            for(int cell=idx+1; cell<idx+num_microcells_per_stack; ++cell){

               // calculate cell id's
               const int cellx = 3*cell+0;
               const int celly = 3*cell+1;
               const int cellz = 3*cell+2;

               // calculate previous cell id's
               const int pcellx = 3*(cell-1)+0;
               const int pcelly = 3*(cell-1)+1;
               const int pcellz = 3*(cell-1)+2;

               // copy array values to temporaries for readability
               st::internal::three_vector_t  m(st::internal::m[cellx],  st::internal::m[celly],  st::internal::m[cellz]); // current cell magnetisations
               st::internal::three_vector_t pm(st::internal::m[pcellx], st::internal::m[pcelly], st::internal::m[pcellz]); // previous cell magnetisations

               const double modm = sqrt(m.x*m.x + m.y*m.y + m.z*m.z);
               const double pmodm = sqrt(pm.x*pm.x + pm.y*pm.y + pm.z*pm.z);

               // Check for zero magnetization in normalization
               if(modm > 1.e-8){
                  m.x = m.x/modm;
                  m.y = m.y/modm;
                  m.z = m.z/modm;
               }
               else{
                  m.x = 0.0;
                  m.y = 0.0;
                  m.z = 0.0;
               }
               if(pmodm > 1.e-8){
                  pm.x = pm.x/modm;
                  pm.y = pm.y/modm;
                  pm.z = pm.z/modm;
               }
               else{
                  pm.x = 0.0;
                  pm.y = 0.0;
                  pm.z = 0.0;
               }

               //---------------------------------------------------------------------
               // Step 1 calculate coordinate transformation m -> m' for current cell
               //---------------------------------------------------------------------

               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-8){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m, itm);

                  // Calculate basis vectors
                  b1 = transform_vector(bz, itm);
                  b2 = transform_vector(bx, itm);
                  b3 = transform_vector(by, itm);
               }
               else{
                  // set default rotational frame
                  b1 = bz;
                  b2 = bx;
                  b3 = by;
               }
               //---------------------------------------------------------------------
               // Step 2 determine coefficients for the spin accumulation
               //---------------------------------------------------------------------

               // Initialise temporary constants
               const double Bc = st::internal::beta_cond[cell]; // beta
               const double Bd = st::internal::beta_diff[cell]; // beta_prime
               const double Do = st::internal::diffusion[cell];
               const st::internal::three_vector_t jm0(st::internal::j[pcellx],st::internal::j[pcelly],st::internal::j[pcellz]);

               //  Calculate gradient dsacc/dx
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m.x*m.x - twoDo;     M.xy = preD*m.y*m.x;             M.xz = preD*m.z*m.x;
               M.yx = preD*m.x*m.y;             M.yy = preD*m.y*m.y - twoDo;     M.yz = preD*m.z*m.y;
               M.zx = preD*m.x*m.z;             M.zy = preD*m.y*m.z;             M.zz = preD*m.z*m.z - twoDo;

               V.x = jm0.x - Bc*je*m.x;
               V.y = jm0.y - Bc*je*m.y;
               V.z = jm0.z - Bc*je*m.z;

               const st::internal::three_vector_t divm_0 = gaussian_elimination(M, V);

               // Calculate mp(0), c and d
               const double i_lsdl = 1.0/st::internal::lambda_sdl[cell];
               const double mp_inf = st::internal::sa_infinity[cell];
               const double a = st::internal::a[cell];
               const double b = st::internal::b[cell];
               const double two_a = 2.0*a;
               const double two_b = 2.0*b;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               const st::internal::three_vector_t coeff = gaussian_elimination(M, V);

               const double mp_0 = coeff.x;
               const double c    = coeff.y;
               const double d    = coeff.z;

               //------------------------------------
               // Step 3 calculate spin accumulation
               //------------------------------------

               const double x = st::internal::micro_cell_thickness*1.0e-10; // Convert to metres
               const double cos_bx = cos(b*x);
               const double sin_bx = sin(b*x);
               const double e_xsdl = exp(-x*i_lsdl);
               const double e_ax   = exp(-a*x);
               const double prefac = (2.0*e_ax);

               const double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               const double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               const double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);

               //convert mp and m_perp
               const double sax = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
               const double say = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
               const double saz = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;

               //--------------------------------------------
               // Step 4 calculate the spin current (jm_end)
               //--------------------------------------------
               const double ac_bd = a*c + b*d;
               const double ad_bc = a*d - b*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m.x*divsax + m.y*divsay + m.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m.z*dot;

               const double jmx = Bc*je*m.x - twoDo*pre_jmx;
               const double jmy = Bc*je*m.y - twoDo*pre_jmy;
               const double jmz = Bc*je*m.z - twoDo*pre_jmz;

               // Save values for the spin accumulation
               st::internal::sa[cellx] = sax;
               st::internal::sa[celly] = say;
               st::internal::sa[cellz] = saz;

               // Save values for the spin current
               st::internal::j[cellx] = jmx;
               st::internal::j[celly] = jmy;
               st::internal::j[cellz] = jmz;

               // Calculate spin torque energy for cell (Joules)
               st::internal::spin_torque[cellx] = microcell_volume * st::internal::sd_exchange[cell] * sax * i_e * i_muB;
               st::internal::spin_torque[celly] = microcell_volume * st::internal::sd_exchange[cell] * say * i_e * i_muB;
               st::internal::spin_torque[cellz] = microcell_volume * st::internal::sd_exchange[cell] * saz * i_e * i_muB;


              //--------------------------------------------
              // Step 5 calculate the spin torque of each cell
              //--------------------------------------------

               //convert M of previous cell into basis b1, b2, b3

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm.x;
               V.y = pm.y;
               V.z = pm.z;

               const st::internal::three_vector_t pm_basis = gaussian_elimination(M, V);

               const double pm_b1 = pm_basis.x;
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast

               const double prefac_sc = microcell_volume * st::internal::sd_exchange[cell] * i_e * i_muB;
               const double plus_perp =  (pm_b2*pm_b2 + pm_b3*pm_b3);
               const double minus_perp = (pm_b2*pm_b2 - pm_b3*pm_b3);
	       double a_constant;
               double b_constant;


                if( ( plus_perp <= 1.0e-7 ) ){
                    a_constant = 0.0;
                    b_constant = 0.0; }

                else{
                    a_constant  = prefac_sc*(sa_perp2*pm_b3 - sa_perp3*pm_b2)/plus_perp;
                    b_constant  = prefac_sc*(sa_perp2*pm_b2 + sa_perp3*pm_b3)/plus_perp;
                }

                double SxSp[3], SxSxSp[3];
                st::internal::coeff_ast[cell]  = a_constant;
                st::internal::coeff_nast[cell] = b_constant;

                SxSp[0]=(m.y*pm.z-m.z*pm.y);
                SxSp[1]=(m.z*pm.x-m.x*pm.z);
                SxSp[2]=(m.x*pm.y-m.y*pm.x);

                SxSxSp[0]= (m.y*SxSp[2]-m.z*SxSp[1]);
                SxSxSp[1]= (m.z*SxSp[0]-m.x*SxSp[2]);
                SxSxSp[2]= (m.x*SxSp[1]-m.y*SxSp[0]);

                //calculate directly from J(Sxm)
                st::internal::total_ST[cellx] = prefac_sc*(m.y*saz-m.z*say);
                st::internal::total_ST[celly] = prefac_sc*(m.z*sax-m.x*saz);
                st::internal::total_ST[cellz] = prefac_sc*(m.x*say-m.y*sax);

                st::internal::ast[cellx] = -a_constant*SxSxSp[0];
                st::internal::ast[celly] = -a_constant*SxSxSp[1];
                st::internal::ast[cellz] = -a_constant*SxSxSp[2];

                st::internal::nast[cellx] = b_constant*SxSp[0];
                st::internal::nast[celly] = b_constant*SxSp[1];
                st::internal::nast[cellz] = b_constant*SxSp[2];



            } // end of cell loop
         } // end of stack loop

         // Reduce all microcell spin torques on all nodes
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::spin_torque[0],st::internal::spin_torque.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::total_ST[0],st::internal::total_ST.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif

         st::internal::output_microcell_data();

         return;
      }

   } // end of namespace internal
} // end of namespace st
