#pragma once
#include <cmath>
#include <functional>

using function = double(double);

double Riman_integral(function func, double x0, double x1, double h) {
  double res = 0;

  for (double x = x0; x < x1; x+=h ) {
    res += func(x + h / 2);
  }

  return res * h;
}

int foo() { return 10; }
