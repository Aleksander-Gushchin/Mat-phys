#include "lib.h"
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>






int main()
{
  const int size = 1001;
  std::vector<double> vec(size, 0.0);
  std::iota(vec.begin(), vec.end(), 0);
  std::for_each(vec.begin(), vec.end(), [&](double& x) {x /= (size - 1); });
  SimpleBasis basis(std::move(vec));

  auto p = [](double x) {return 2 + 0.5 * sin(2 * x); };
  double q = 1.5;
  auto f = [](double x) { return x * (1 - x); };

  PQCalc calc(basis, p, f, q);
  
  PQCalc cl(basis, [](double x) {return sin(x); }, [](double x) { return x*(1 - x); },  1.0);

  std::vector<double> res = calc.find_sol();

  for (double v: res) {
    std::cout << v << std::endl;
  }

  return 0;
}
