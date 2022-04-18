#include "lib.h"
#include <iostream>
#include <cmath>


int main()
{
  auto la = [](double x) { return sin(x); };
  std::cout << Riman_integral(sin, 0, 1, 1e-2) << "\n";
  return 0;
}