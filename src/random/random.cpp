#include "random.hpp"
#include <cmath>

using std::log;
using std::sqrt;

namespace mtrandom
//==========================================================
// Namespace mt random
//==========================================================
{
	
	double x1,x2,w;
	double number1;
	double number2;
	bool logic=false;
	MTRand grnd; // single sequence of random numbers

double gaussian(){
	using namespace mtrandom;

	if(logic==false){
		for(;;){
			x1 = 2.0*grnd() - 1.0;
			x2 = 2.0*grnd() - 1.0;
			w = x1*x1 + x2*x2;
			if(w<1.0) break;
		}
		w=sqrt((-2.0*log(w))/w);
		number1 = w*x1;
		number2 = w*x2;
		logic=true;
		return number1;
	}
	else {
		logic=false;
		return number2;
		}
}
}
