//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2025. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
#ifndef VTYPES_H_
#define VTYPES_H_
// C++ standard library headers
#include <vector>
#include <initializer_list>

// Vampire headers
#include "errors.hpp"
#include "vio.hpp"


namespace vtype{

   //---------------------------------
   // struct for storing unit vectors
   //---------------------------------
   class vec_t{

      public:
         double x;
         double y;
         double z;

      // constructors
      vec_t():x(0.0), y(0.0), z(0.0){}; // default
      vec_t(double xi):x(xi), y(xi), z(xi){}; // single value initialisation
      vec_t(double xi, double yi, double zi):x(xi), y(yi), z(zi){}; // 3-value initialisation
      vec_t(std::vector<double> initial_values){
         if(initial_values.size() == 3){
            x = initial_values[0];
            y = initial_values[1];
            z = initial_values[2];
         }
         else{
            std::cerr << "Error: attempting to initialise vtypes::vec_t with a std::vector not containing three values" << std::endl;
            zlog << zTs() << "Error: attempting to initialise vtypes::vec_t with a std::vector not containing three values" << std::endl;
            err::vexit();
         }
      }; // std::vector initialisation
      /*vec_t(std::initializer_list<double> initial_values):{
         if(initial_values.size() == 3){
            x = initial_values[0];
            y = initial_values[1];
            z = initial_values[2];
         }
         else{
            std::cerr << "Error: attempting to initialise vtypes::vec_t with an initializer list {} not containing three values" << std::endl;
            zlog << zTs() << "Error: attempting to initialise vtypes::vec_t with an initializer list {} not containing three values" << std::endl;
            err::vexit();
         }
      }; // std::vector initialisation*/

      // Copy constructor
      vec_t(const vec_t& other) : x(other.x), y(other.y), z(other.z) {}

      // Copy assignment operator
      vec_t& operator=(const vec_t& other) {
          if (this != &other) {
              x = other.x;
              y = other.y;
              z = other.z;
          }
          return *this;
      }

      // Move constructor
      vec_t(vec_t&& other) noexcept : x(other.x), y(other.y), z(other.z) {
          other.x = 0.0;
          other.y = 0.0;
          other.z = 0.0;
      }

      // Move assignment operator
      vec_t& operator=(vec_t&& other) noexcept {
         if (this != &other) {
            x = other.x;
            y = other.y;
            z = other.z;
            other.x = 0.0;
            other.y = 0.0;
            other.z = 0.0;
         }
         return *this;
      }

   };

   //---------------------------------------------
   // simple initialised class for set variables
   //---------------------------------------------
   class set_double_t{

   private:
      double value; // value
      bool setf; // flag specifiying variable has been set

   public:
      // class functions
      // constructor
      set_double_t() : value(0.0), setf(false) { }

      // setting function
      void set(double in_value){
         value = in_value;
         setf = true;
      };

      // get value function
      double get(){ return value; };
      // check if variable is set
      bool is_set(){ return setf; };

   };

}

#endif /*VTYPES_H_*/
