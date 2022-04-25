#include "lib.h"

double TridMatrix::get_A(int i)
{
  return data[i];
}

double TridMatrix::get_B(int i)
{
  return data[size + i];
}

double TridMatrix::get_C(int i)
{
  return data[2 * size + i + 1];
}

TridMatrix::TridMatrix(int _size) : size(_size)
{
  data = std::vector<double>(3 * size, 0); // a[], b[], c[]
}

void TridMatrix::set(int i, int j, double val)
{
  data[(j - i + 1) * size + j] = val;
}

TridMatrix::~TridMatrix()
{
}

std::vector<double> TridMatrix::trid_matrix_alg(std::vector<double>& vec)
{
  return std::vector<double>();
}

double integrate(std::function<double(double)> func, double x0, double x1, double h)
{
  double res = 0;

  for (double x = x0; x < x1; x += h) {
    res += func(x + h / 2);
  }

  return res * h;
}

double integrate(std::function<double(double)> func, double x0, double x1, int n)
{
  double res = 0;

  double h = (x1 - x0) / n;

  for (double x = x0; x < x1; x += h) {
    res += func(x + h / 2);
  }

  return res * h;
}

BaseBasis::BaseBasis(const std::vector<double>& vec) : range(vec)  {}

BaseBasis::BaseBasis(std::vector<double>&& vec) : range(std::move(vec)) {}

BaseBasis::BaseBasis(BaseBasis&& basis) : range(std::move(basis.range))
{
}

BaseBasis::BaseBasis(const BaseBasis& basis) : range(basis.range)
{
}

int BaseBasis::size() const
{
  return range.size();
}

BaseBasis::~BaseBasis()
{
}
