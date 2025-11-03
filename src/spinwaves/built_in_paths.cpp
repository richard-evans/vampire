
//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sergiu Ruta 2022. All rights reserved.
//
//   Email: sergiu.ruta@york.ac.uk j.r.hirst@shu.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "spinwaves.hpp"
#include "unitcell.hpp"
#include "internal.hpp"
#include "errors.hpp"
#include "atoms.hpp"
#include "vmpi.hpp"
#include "vio.hpp"


// sw module headers
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

namespace spinwaves {

    namespace internal {

        void path_sc(){
            std::vector<double> dummy_x = {0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5};
            std::vector<double> dummy_y = {0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5};
            std::vector<double> dummy_z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0};
            spinwaves::internal::pathx.swap(dummy_x);
            spinwaves::internal::pathy.swap(dummy_y);
            spinwaves::internal::pathz.swap(dummy_z);
        }

        void path_bcc(){
            std::vector<double> dummy_x = {0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5};
            std::vector<double> dummy_y = {0.0, 1.0, 1.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 1.0, 0.5, 0.5};
            std::vector<double> dummy_z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0};
            spinwaves::internal::pathx.swap(dummy_x);
            spinwaves::internal::pathy.swap(dummy_y);
            spinwaves::internal::pathz.swap(dummy_z);
        }

        void path_fcc(){
            std::vector<double> dummy_x = {0.00, 0.00, 0.00, 0.25, 0.75, 0.00, 0.00, 0.50, 0.50, 0.50, 0.50, 0.00};
            std::vector<double> dummy_y = {0.00, 1.00, 1.00, 1.00, 0.75, 0.00, 0.00, 0.50, 0.50, 1.00, 1.00, 1.00};
            std::vector<double> dummy_z = {0.00, 0.00, 0.00, 0.25, 0.00, 0.00, 0.00, 0.50, 0.50, 0.00, 0.00, 0.00};
            spinwaves::internal::pathx.swap(dummy_x);
            spinwaves::internal::pathy.swap(dummy_y);
            spinwaves::internal::pathz.swap(dummy_z);
        }

        void path_hcp(){
            std::vector<double> dummy_x = {0.0, 0.5, 0.5, 1.0     , 1.0     , 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 1.0     , 1.0     , 0.0, 0.5, 0.5, 1.0     , 1.0};
            std::vector<double> dummy_y = {0.0, 0.5, 0.5, 0.333333, 0.333333, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.333333, 0.333333, 0.0, 0.5, 0.5, 0.333333, 0.333333};
            std::vector<double> dummy_z = {0.0, 0.0, 0.0, 0.0     , 0.0     , 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5     , 0.5     , 0.5, 0.5, 0.0, 0.5     , 0.0};
            spinwaves::internal::pathx.swap(dummy_x);
            spinwaves::internal::pathy.swap(dummy_y);
            spinwaves::internal::pathz.swap(dummy_z);
        }

        void path_mn2au(){
            std::vector<double> dummy_x = {0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.44444444, 0.55555555, 0.0, 0.5, 0.55555555, 0.44444444, 0.0};
            std::vector<double> dummy_y = {0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0       , 0.0       , 0.0, 0.5, 0.44444444, 0.44444444, 0.0};
            std::vector<double> dummy_z = {0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 1.0, 1.0, 1.0       , 0.0       , 0.0, 0.0, 0.0       , 1.0       , 1.0};
            spinwaves::internal::pathx.swap(dummy_x);
            spinwaves::internal::pathy.swap(dummy_y);
            spinwaves::internal::pathz.swap(dummy_z);
        }

        void path_rocksalt(){
            std::vector<double> dummy_x = {0.0, 0.0, 0.0, 0.25,  0.75,  0.0,  0.0,  0.5,  0.5,  0.5,  0.5,  0.0};
            std::vector<double> dummy_y = {0.0, 1.0, 1.0, 1.0 , 0.75 , 0.0 , 0.0 , 0.5 , 0.5 , 1.0 , 1.0 , 1.0};
            std::vector<double> dummy_z = {0.0, 0.0, 0.0, 0.25,  0.0 , 0.0 , 0.0 , 0.5 , 0.5 , 0.0 , 0.0 , 0.0};
            spinwaves::internal::pathx.swap(dummy_x);
            spinwaves::internal::pathy.swap(dummy_y);
            spinwaves::internal::pathz.swap(dummy_z);
        }

        void path_heusler(){
            std::vector<double> dummy_x = {0.00, 0.00, 0.00, 0.25, 0.75, 0.00, 0.00, 0.50, 0.50, 0.50, 0.50, 0.00};
            std::vector<double> dummy_y = {0.00, 1.00, 1.00, 1.00, 0.75, 0.00, 0.00, 0.50, 0.50, 1.00, 1.00, 1.00};
            std::vector<double> dummy_z = {0.00, 0.00, 0.00, 0.25, 0.00, 0.00, 0.00, 0.50, 0.50, 0.00, 0.00, 0.00};
            spinwaves::internal::pathx.swap(dummy_x);
            spinwaves::internal::pathy.swap(dummy_y);
            spinwaves::internal::pathz.swap(dummy_z);
        }



        int gcd(int a, int b){
            a = abs(a);b = abs(b);
            if (a == 0)
                return b;
            if (b == 0)
                return a;
            if (a == b)
                return a;
            if (a > b)
                return gcd(a - b, b);
            return gcd(a, b - a);
        }

    }

}



    // void bcc() {

    //     GAMMA H  :  [0.0, 0.0, 0.0] => [0.0, 1.0, 0.0]
    // H N  :          [0.0, 1.0, 0.0] => [0.5, 0.5, 0.0]
    //     N GAMMA  :  [0.5, 0.5, 0.0] => [0.0, 0.0, 0.0]
    //     GAMMA P  :  [0.0, 0.0, 0.0] => [0.5, 0.5, 0.5]
    //     P H  :      [0.5, 0.5, 0.5] => [0.0, 1.0, 0.0]
    //     P N  :      [0.5, 0.5, 0.5] => [0.5, 0.5, 0.0]

    // }

    // void sc() {

    //     GAMMA X  :  [0.0, 0.0, 0.0] => [0.0, 0.5, 0.0]
    //     X M  :  [0.0, 0.5, 0.0] => [0.5, 0.5, 0.0]
    //     M GAMMA  :  [0.5, 0.5, 0.0] => [0.0, 0.0, 0.0]
    //     GAMMA R  :  [0.0, 0.0, 0.0] => [0.5, 0.5, 0.5]
    //     R X  :  [0.5, 0.5, 0.5] => [0.0, 0.5, 0.0]
    //     R M  :  [0.5, 0.5, 0.5] => [0.5, 0.5, 0.0]

    // }

    // void fcc() {

    //     GAMMA X  :  [0.0, 0.0, 0.0] => [0.0, 1.0, 0.0]
    //     X U  :  [0.0, 1.0, 0.0] => [0.25, 1.0, 0.25]
    //     K GAMMA  :  [0.75, 0.75, 0.0] => [0.0, 0.0, 0.0]
    //     GAMMA L  :  [0.0, 0.0, 0.0] => [0.5, 0.5, 0.5]
    //     L W  :  [0.5, 0.5, 0.5] => [0.5, 1.0, 0.0]
    //     W X  :  [0.5, 1.0, 0.0] => [0.0, 1.0, 0.0]

    // }

    // void hcp() {

    //     GAMMA M  [0 0 0] => [0.5 0.5 0]
    //     M K  :   [0.5 0.5 0] =>[1 0.333333 0]
    //     K GAMMA  [1 0.333333 0] =>[0 0 0 ]
    //     GAMMA A  [0 0 0 ] => [0 0 0.5 ]
    //     A L  :   [0 0 0.5 ] => [0.5 0.5 0.5 ]
    //     L H  : [0.5 0.5 0.5 ] =>[1 0.333333 0.5 ]
    //     H A  : [1 0.333333 0.5 ] =>[0 0 0.5 ]
    //     L M  : [0.5 0.5 0.5 ] =>[0.5 0.5 0]
    //     H K  : [1 0.333333 0.5 ] => [1 0.333333 0]

    // }

    // void mn2au(){

    //     GAMMA X  :  [0.0, 0.0, 0.0] => [0.5, 0.5, 0.0]
    //     X P  :  [0.5, 0.5, 0.0] => [0.5, 0.5, 0.5]
    //     P N  :  [0.5, 0.5, 0.5] => [0.5, 0.0, 0.5]
    //     N GAMMA  :  [0.5, 0.0, 0.5] => [0.0, 0.0, 0.0]
    //     GAMMA M  :  [0.0, 0.0, 0.0] => [0.0, 0.0, 1.0]
    //     M S  :  [0.0, 0.0, 1.0] => [0.444444444444444, 0.0, 1.0]
    //     S_0 GAMMA  :  [0.55555555555555, 0.0, 0.0] => [0.0, 0.0, 0.0]
    //     X R  :  [0.5, 0.5, 0.0] => [0.555555555555555, 0.4444444444444444, 0.0]
    //     G M  :  [0.4444444444444444, 0.4444444444444444, 1.0] => [0.0, 0.0, 1.0]
