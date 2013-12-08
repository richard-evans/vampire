#include <iostream>
#include <cmath>

double tanf(double x){
  return tanh((x-200.0)/300.0);
}

int main(){

  std::cout << 1501 << std::endl;
  for(int i=0; i<1501; i++){
    std::cout << i << "\t" << (tanf(double(i))-tanf(0.0))*0.5 +1.0 << std::endl;
  }

  return 0;

}
