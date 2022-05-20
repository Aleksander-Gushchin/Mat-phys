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

  //for (int j = 0; j < 4; ++j) {
  //  for (int i = 0; i <= 10; ++i)
  //    std::cout << basis.get_deriv(j, 0.1 * i) << " ";
  //  std::cout << "\n";

  //  for (int i = 0; i <= 10; ++i)
  //    std::cout << basis.get_value(j, 0.1 * i) << " ";
  //  std::cout << "\n";
  //  std::cout << "\n";
  //  std::cout << "\n";
  //}
  auto p = [](double x) {return 2 + 0.5 * sin(2 * x); };
  double q = 1.5;
  auto f = [](double x) { return x * (1 - x); };

  PQCalc calc(basis, p, f, q);
  
  PQCalc cl(basis, [](double x) {return sin(x); }, [](double x) { return x*(1 - x); },  1.0);

  std::vector<double> res = calc.find_sol();

  for (double v: res) {
    std::cout << v << std::endl;
  }

  //std::cout << integrate([](double x) {return sin(x); }, 0, 1, 1e-2) << "\n";

  /*TridMatrix test(4);
  test.setXi1(-1./2.);
  test.setXi2(-1./4.);
  test.set_A(0, 1.);
  test.set_C(0, 10.);
  test.set_B(0, -5.);
  test.set_A(1, 1.);
  test.set_C(1, -5.);
  test.set_B(1, 2.);
  std::vector<double> fi(4);
  fi[0] = -5./2.;
  fi[1] = -18.;
  fi[2] = -40.;
  fi[3] = -27./4.;

  for(double v: test.trid_matrix_alg(fi)) {
    std::cout << v << std::endl;
  }*/
  return 0;
}
