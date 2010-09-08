#ifndef VCUDA_H_
#define VCUDA_H_

// Includes
#include <vector>

#ifdef CUDA
namespace vcuda{

  extern int LLG(const int);  ///< Local CPU ID

}
#endif

#endif // VCUDA_H_


