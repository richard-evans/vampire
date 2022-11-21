// Program headers
#ifndef LSF_MC_H_
#define LSF_MC_H_

namespace montecarlo{

    // Declare spin vectors
    extern std::vector<double> mod_S;
    extern std::vector<double> x_S;
    extern std::vector<double> y_S;
    extern std::vector<double> z_S;
    extern std::vector<double> mod_S_initial;

    extern bool mc_set;

    extern std::vector<double> mu_container;

    extern int process;
    extern std::vector <int> mod_S_round;

}
#endif /*LSF_MC_H_*/

