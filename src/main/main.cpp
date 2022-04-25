#include "lib.h"
#include <iostream>
#include <cmath>







int main()
{
  SimpleBasis basis({ 0.0, 0.33, 0.66, 1 });

  for(int i = 0; i < 100; ++i)
    std::cout << basis.get_value(0, 0.01*i) << " ";
  std::cout << "\n";

  
  PQCalc cl(basis, [](double x) {return sin(x); }, [](double x) { return x*(1 - x); },  1.0);

  cl.find_sol();

  //std::cout << integrate([](double x) {return sin(x); }, 0, 1, 1e-2) << "\n";
  return 0;
}
