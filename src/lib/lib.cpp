#include "lib.h"
#include <stdexcept>

void TridMatrix::set_A(int i, double val) {
  data[i] = val;
}

void TridMatrix::set_B(int i, double val) {
  data[2 * (size - 2) + i] = val;
}

void TridMatrix::set_C(int i, double val) {
  data[i + (size - 2)] = val;
}

double TridMatrix::get_A(int i)
{
  return data[i];
}

double TridMatrix::get_B(int i)
{
  return data[2 * (size - 2) + i];
}

double TridMatrix::get_C(int i)
{
  return data[(size - 2) + i];
}

TridMatrix::TridMatrix(int _size) : size(_size)
{
  data = std::vector<double>(3 * (size - 2), 0); // a[], b[], c[]
}

void TridMatrix::set(int i, int j, double val)
{
  data[(j - i + 1) * size + j] = val;
}

void TridMatrix::setXi1(double val) {
  xi1 = val;
}

void TridMatrix::setXi2(double val) {
  xi2 = val;
}

TridMatrix::~TridMatrix()
{
}


// подразумеваем, что все данные внесены учитыая знаки, даже учтены в vec
std::vector<double> TridMatrix::trid_matrix_alg(std::vector<double>& vec)
{
  if (vec.size() != size) {
    throw std::invalid_argument("size of vec != size of matrix");
  }
  mu1 = vec.front();
  mu2 = vec.back();

  std::vector<double> alpha(size - 1);
  std::vector<double> beta(size - 1);
  alpha[0] = xi1;
  beta[0] = mu1;
  // прямой ход прогонки
  for (int i = 1; i < size - 1; i++) {
    double A = get_A(i - 1);
    double B = get_B(i - 1);
    double C = -get_C(i - 1);
    alpha[i] = B / (C - A * alpha[i - 1]);
    beta[i] = (-vec[i] + A * beta[i - 1]) / (C - A * alpha[i - 1]);
  }

  std::vector<double> V(size);
  V[size - 1] = (-xi2 * beta.back() - mu2) / (xi2 * alpha.back() - 1);
  // обратный ход прогонки
  for (int i = size - 2; i >= 0; i--) {
    V[i] = alpha[i] * V[i + 1] + beta[i];
  }

  return V;
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
  if (!func)
    throw std::invalid_argument("empty lambda");

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
