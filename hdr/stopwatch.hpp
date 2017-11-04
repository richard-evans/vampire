#ifndef __STOPWATCH_H__
#define __STOPWATCH_H__

#include <chrono>

class stopwatch_t {
   public:

      typedef std::chrono::high_resolution_clock::time_point tpoint_t;
      tpoint_t t1;
      tpoint_t t2;

      stopwatch_t()
      {};

      ~stopwatch_t() {};

      void start() {
         t1 = std::chrono::high_resolution_clock::now();
      }

      double elapsed_seconds() {
         t2 = std::chrono::high_resolution_clock::now();

         std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
         return time_span.count();
      }

};

#endif
