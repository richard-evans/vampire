// Standard libraries
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <random>
#include <fstream>

// Vampire header files
#include "internal.hpp"
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "montecarlo.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmath.hpp"
#include "vio.hpp"
#include "internal.hpp"
#include "lsf_mc.hpp"
#include "exchange.hpp"
#include "../simulate/internal.hpp"


namespace montecarlo{

    std::vector<double> mod_S;
    std::vector<double> x_S;
    std::vector<double> y_S;
    std::vector<double> z_S;
    std::vector<double> mod_S_initial;
    std::vector<int> statistics_moves_l;
    std::vector<int> statistics_reject_l;
    bool mc_set=false;

    // *** TESTING ***
    std::ofstream MyFile("spinlength");
    //std::ofstream MyFile("mu_d");
    double simtemp=1;
    std::vector<double> chi;
    double exch;
    std::vector<double> mu_D;
    double mu_D_avg=0;
    double tempcheck=1;
    double exchsum;
    double mu_i;
    int process;
    std::vector <int> mod_S_round;
    int counter=0;
    int loc_counter=0;
    double mod_S_pt;
    double mod_S_fe;
    double spinlength[2]={0,0};
    double ptdir[3];

    void mcinit(){

        // Initialize spin values
        x_S.resize(atoms::num_atoms,0.0);
        y_S.resize(atoms::num_atoms,0.0);
        z_S.resize(atoms::num_atoms,1.0);
        mod_S_initial.resize(atoms::num_atoms,1.0);
        mod_S.resize(atoms::num_atoms,1.0);

        mc_set=true;

        //TESTING ******
        chi.resize(atoms::num_atoms,0.0);
        mu_D.resize(atoms::num_atoms,0.0);
        mod_S_round.resize(401, 0);
        mu_i=mp::material[1].mu_s_SI;
        mod_S_pt=0;
        for(int atom=0; atom<atoms::num_atoms; atom++){
            exch = exchange::single_spin_energy(atom, atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]);
            chi[atom]=mp::material[atoms::type_array[atom]].mu_s_SI/exch;
        }
        for(int i=0; i<3; i++){
            ptdir[i]=0;
        }

    }

    /// Angle move
    /// Move spin within cone near old position
    void mc_angle_transverse(const std::vector<double>& old_spin, std::vector<double>& new_spin, const double angle, const int atom){

        double factorx = mtrandom::grnd();
        double factory = mtrandom::grnd();
        double factorz = mtrandom::grnd();
        double p;
        p = mtrandom::grnd();
        if(p >= 0.5) factorx *= -1;
        p = mtrandom::grnd();
        if(p >= 0.5) factory *= -1;
        p = mtrandom::grnd();
        if(p >= 0.5) factorz *= -1;
        
        new_spin[0] = old_spin[0] + factorx * angle;
        new_spin[1] = old_spin[1] + factory * angle;
        new_spin[2] = old_spin[2] + factorz * angle;


        return;

    }

    int lsf_mc_step(){

        // Enable calling of routine if error checking is activated
        if(err::check==true){std::cout << "sim::lsf_mc has been called" << std::endl;}

        const int num_atoms=atoms::num_atoms;
        
        // Variable declarations
        double accept_reject_probability_exch;
        double accept_reject_probability_l;
        double P;
        bool tf;
        double move_type;
        double k_B=1.38064e-23;
        double k_BT;
        double total_E1;
        double total_E2;
        double delta_E;
        double Eold=0.0;
        double Enew=0.0;
        double DE=0.0;
        double P_t;
        double deltaS;
        double atom;
        double x_S_initial, y_S_initial, z_S_initial;
        double re;
        double defaults[2]={0.1,0.1};
        double default_upper[2]={1.0, 1.0};
        double default_lower[2]={1e-8, 1e-8};
        int testatom = 3;
        double samp=0;
        bool checking = false;

        double statistics_moves = 0.0;
        double statistics_reject = 0.0;
        statistics_moves_l.resize(mp::num_materials,0.0);
        statistics_reject_l.resize(mp::num_materials,0.0);

        if(montecarlo::mc_set==false) montecarlo::mcinit();

        // *** TESTING ***
        //double spinlength[2]={0,0};
        //double mod_S_fe=0;
        //double mod_S_pt=0;
        for (int atom = 0; atom < num_atoms; atom++){

            const int imaterial = atoms::type_array[atom];

            spinlength[imaterial] = sqrt(atoms::x_spin_array[atom]*atoms::x_spin_array[atom] + atoms::y_spin_array[atom]*atoms::y_spin_array[atom] + atoms::z_spin_array[atom]*atoms::z_spin_array[atom]);
            
            if(imaterial==1){
                ptdir[0]+=atoms::x_spin_array[atom]/spinlength[imaterial];
                ptdir[1]+=atoms::y_spin_array[atom]/spinlength[imaterial];
                ptdir[2]+=atoms::z_spin_array[atom]/spinlength[imaterial];
            }

            int spinround = floor((spinlength[imaterial]*100)+0.5);
            if(spinround > 400) spinround=400;
            if(spinround < 0) spinround=0;
            mod_S_round[spinround]+= 1;

            if(imaterial==0){
                mod_S_fe += spinlength[imaterial];
            }
            if(imaterial==1){
                mod_S_pt += spinlength[imaterial];
            }

        }
        if(simtemp!=sim::temperature){
            //MyFile << sim::temperature << " " << mod_S_fe/(atoms::num_atoms/mp::num_materials) << " " << mod_S_pt/(atoms::num_atoms/mp::num_materials) << " " << spinlength[0]/(atoms::num_atoms/mp::num_materials) << " " << spinlength[1]/(atoms::num_atoms/mp::num_materials) << std::endl;
            simtemp=sim::temperature;
        }

        // TESTING***********
        mu_D_avg=0.0;
        //exchsum=0.0;
        for (int atom = 0; atom < num_atoms; atom++){
            exch = exchange::single_spin_energy(atom, atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]);
            exchsum+=exch;
            if(atoms::type_array[atom]==1 && atom==515){
                mu_D[atom]=chi[atom]*exch;
                mu_D_avg+=mu_D[atom];
            }
        }
        if(tempcheck!=sim::temperature){
            //MyFile << sim::temperature << " " << ((mu_D_avg/(atoms::num_atoms/mp::num_materials)))/mu_i << " " << exchsum/(atoms::num_atoms/mp::num_materials) << std::endl;
            if(loc_counter==0) loc_counter=1;
            //MyFile << sim::temperature << " " << exchsum/(loc_counter*atoms::num_atoms) << std::endl;
            MyFile << sim::temperature << " " << mod_S_fe/((atoms::num_atoms/mp::num_materials)*loc_counter) << " " << mod_S_pt/((atoms::num_atoms/mp::num_materials)*loc_counter) << std::endl;
            tempcheck=sim::temperature;
            loc_counter=0;
            exchsum=0;
            mod_S_pt=0;
            mod_S_fe=0;
        }
        spinlength[0]=0;
        spinlength[1]=0;

        for (int i=0; i<atoms::num_atoms; i++){

            atom = int(atoms::num_atoms*mtrandom::grnd());

            accept_reject_probability_exch=mtrandom::grnd();
            accept_reject_probability_l=mtrandom::grnd();
    
            // 50/50 to determine if move is transverse or longitudinal
            tf=false;
            move_type=mtrandom::grnd();
            if (move_type>=0.5){
                tf=false;
            }  
            
            // Move is transverse
            if (tf==false){ 
                
                // TESTING*****
                process=0;

                // add one to number of moves counter
                statistics_moves+=1.0;

                // get material id
                const int imaterial=atoms::type_array[atom];

                // Save old spin position
                internal::Sold[0] = atoms::x_spin_array[atom];
                internal::Sold[1] = atoms::y_spin_array[atom];
                internal::Sold[2] = atoms::z_spin_array[atom];

                // Transverse step
                montecarlo::mc_angle_transverse(internal::Sold, internal::Snew, montecarlo::internal::adaptive_sigma, atom);

                atoms::x_spin_array[atom] = internal::Sold[0];
                atoms::y_spin_array[atom] = internal::Sold[1];
                atoms::z_spin_array[atom] = internal::Sold[2];

                // Calculate current energy
                Eold = sim::calculate_spin_energy(atom);

                atoms::x_spin_array[atom] = internal::Snew[0];
                atoms::y_spin_array[atom] = internal::Snew[1];
                atoms::z_spin_array[atom] = internal::Snew[2];

                // Calculate new energy
                Enew = sim::calculate_spin_energy(atom);

                double oldl=sqrt(internal::Sold[0]*internal::Sold[0]+internal::Sold[1]*internal::Sold[1]+internal::Sold[2]*internal::Sold[2]);
                double newl=sqrt(internal::Snew[0]*internal::Snew[0]+internal::Snew[1]*internal::Snew[1]+internal::Snew[2]*internal::Snew[2]);

                // Calculate difference in Joules/mu_B
                DE = (Enew-Eold)*mp::material[imaterial].mu_s_SI;

                double P_t=exp(-DE/(k_B*sim::temperature));

                double A=sim::internal::lsf_second_order_coefficient[imaterial];
	            double B=sim::internal::lsf_fourth_order_coefficient[imaterial];
	            double C=sim::internal::lsf_sixth_order_coefficient[imaterial];
                double loldE=A*oldl*oldl + B*oldl*oldl*oldl*oldl + C*oldl*oldl*oldl*oldl*oldl*oldl;
                double lnewE=A*newl*newl + B*newl*newl*newl*newl + C*newl*newl*newl*newl*newl*newl;
                double oldexch=exchange::single_spin_energy(atom, internal::Sold[0], internal::Sold[1], internal::Sold[2])*mp::material[imaterial].mu_s_SI;
                double newexch=exchange::single_spin_energy(atom, internal::Snew[0], internal::Snew[1], internal::Snew[2])*mp::material[imaterial].mu_s_SI;
                
                //double exchDE=newexch-oldexch;
                //double lDE=lnewE-loldE;

                //double P_exch=exp(-exchDE/(k_B*sim::temperature));
                //double P_l=exp(-lDE/(k_B*sim::temperature));

                if(imaterial==1 && sim::temperature==samp && checking==true) std::cout << atom << " " << oldl << " " << newl << " " << DE << " " << P_t << " " << loldE << " " << lnewE << " " << oldexch << " " << newexch;
                

                if(DE<0){
                    // Copy new spin position
                    atoms::x_spin_array[atom] = internal::Snew[0];
                    atoms::y_spin_array[atom] = internal::Snew[1];
                    atoms::z_spin_array[atom] = internal::Snew[2];
                    if(imaterial==1 && sim::temperature==samp && checking==true) std::cout << " - accepted DE<0" << std::endl;
                } else if (P_t >= accept_reject_probability_exch){
                    // Copy new spin position
                    atoms::x_spin_array[atom] = internal::Snew[0];
                    atoms::y_spin_array[atom] = internal::Snew[1];
                    atoms::z_spin_array[atom] = internal::Snew[2];
                    if(imaterial==1 && sim::temperature==samp && checking==true) std::cout << " - accepted" << std::endl;
                } else{
                    atoms::x_spin_array[atom] = internal::Sold[0];
                    atoms::y_spin_array[atom] = internal::Sold[1];
                    atoms::z_spin_array[atom] = internal::Sold[2];
                    if(imaterial==1 && sim::temperature==samp && checking==true) std::cout << " - rejected" << std::endl;
                    // add one to rejection counter
                    statistics_reject += 1.0;

                }

            }

        }

        // calculate new adaptive step sigma angle
        if(statistics_moves!=0){
            const double last_rejection_rate = statistics_reject / statistics_moves;
            const double factor = 0.5 / last_rejection_rate;
            montecarlo::internal::adaptive_sigma *= factor;
            // check for excessive range (too small angle takes too long to grow, too large does not improve performance) and truncate
            if (montecarlo::internal::adaptive_sigma > 60.0 || montecarlo::internal::adaptive_sigma < 1e-5) montecarlo::internal::adaptive_sigma = 60.0; // FOR NORMAL RUNS
            //if (montecarlo::internal::adaptive_sigma > 60.0 || montecarlo::internal::adaptive_sigma < 1e-10) montecarlo::internal::adaptive_sigma = 30.0; // FOR SINGLE SPIN
        }
        // Longitudinal sigma
        for(int i=0; i<mp::num_materials; ++i){
            if(statistics_moves_l[i]!=0){
                const double last_l_rejection_rate = statistics_reject_l[i]/statistics_moves_l[i];
                const double factor_l = 0.5 / last_l_rejection_rate;
                montecarlo::internal::adaptive_sigma_l[i] *= factor_l; 
                // check for excessive range
                if (montecarlo::internal::adaptive_sigma_l[i] > default_upper[i] || montecarlo::internal::adaptive_sigma_l[i] < default_lower[i]) montecarlo::internal::adaptive_sigma_l[i]=defaults[i];
            }
        }
        // Save statistics to sim namespace variable
        double l_moves;
        double l_reject;
        for(int i=0; i<mp::num_materials; ++i){
            l_moves += statistics_moves_l[i];
            l_reject += statistics_reject_l[i];
        }
        sim::mc_statistics_moves += statistics_moves+l_moves;
        sim::mc_statistics_reject += statistics_reject+l_reject;

        counter++;
        loc_counter++;
        
        
        // Write spin length sampling data
        if(counter==sim::total_time){
            for (int i=0; i<401; i++){
                //MyFile << i << " " << mod_S_round[i] << std::endl;  
            }
            //std::cout << ptdir[0]/(counter*2662) << " " << ptdir[1]/(counter*2662) << " " << ptdir[2]/(counter*2662) << std::endl;
        }
        
    }

}