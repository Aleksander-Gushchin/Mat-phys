#pragma once
#include <cmath>
#include <vector>
#include <functional>
#include <iostream>
#include <stdexcept>


double integrate(std::function<double(double)> func, double x0, double x1, double h);
double integrate(std::function<double(double)> func, double x0, double x1, int n);

class TridMatrix {
private:
  std::vector<double> data;
  int size;

  double get_A(int i);
  double get_B(int i);
  double get_C(int i);
public:
  TridMatrix(int _size);
  void set(int i, int j, double val); 
  ~TridMatrix();

  std::vector<double> trid_matrix_alg(std::vector<double>& vec); // метод прогонки
};

class BaseBasis {
protected:
  std::vector<double> range;
public:
  BaseBasis(const std::vector<double>& vec);
  BaseBasis(std::vector<double>&& vec);
  BaseBasis(BaseBasis&& basis);
  BaseBasis(const BaseBasis& basis);
  int size() const;
  virtual std::pair<double, double> get_range(int i) const = 0;
  virtual double get_value(int i, double x) const = 0;
  virtual double get_deriv(int i, double x) const = 0;

  double get_func(const std::vector<double>& c_vec, double x) {
    if (range.size() != c_vec.size())
      throw std::invalid_argument("Wrong size of c");

    int size = range.size();
    double res = 0;

    for (int i = 0; i < size; ++i) {
      res += get_value(i, x);
    }

    return res;
  }

  virtual ~BaseBasis();
};


class SimpleBasis : public BaseBasis {
public:
  SimpleBasis(const std::vector<double>& vec) : BaseBasis(vec) {};
  SimpleBasis(std::vector<double>&& vec) : BaseBasis(std::move(vec)) {};
  SimpleBasis(SimpleBasis&& basis) : BaseBasis(std::move(basis)) {};
  SimpleBasis(const SimpleBasis& basis) : BaseBasis(basis) {};

  //return x[i-1], x[i+1] if posible
  std::pair<double, double> get_range(int i) const {
    return std::make_pair(range[std::max(0, i - 1)], range[std::min(i + 1, size() - 1)]);
  }

  double get_value(int i, double x)  const override  {
    if (i < 0 || i  > size() - 1)
      return 0;
    if (i - 1 >= 0 && x >= range[i - 1] && x <= range[i])
      return (x - range[i - 1]) / (range[i] - range[i - 1]);
    if (i + 1 <= size() - 1 && x >= range[i] && x <= range[i + 1])
      return (range[i + 1] - x) / (range[i + 1] - range[i]);
    return 0;

  };

  double get_deriv(int i, double x) const override {
    if (i < 0 || i > size() - 1)
      return 0;
    if (i - 1 >= 0 && x >= range[i - 1] && x <= range[i])
      return 1 / (range[i] - range[i - 1]);
    if (i + 1 <= size() - 1 && x >= range[i] && x <= range[i + 1])
      return -1 / (range[i + 1] - range[i]);
    return 0;
  }

};


using function = std::function<double(double)>;

template<typename T>
class BaseCalc {
protected:
  T basis;
public:
  BaseCalc() = delete;
  BaseCalc(const T& _basis) : basis(_basis) {}
  BaseCalc(T&& _basis) : basis(std::move(_basis)) {};

  virtual std::pair<double, double> get_first_a() = 0;
  virtual std::pair<double, double> get_last_a() = 0;
  virtual double get_a(int i, int j) = 0;

  virtual double get_l(int i) = 0;

  std::vector<double> find_sol() {
    //
    int size = basis.size();
    std::vector<double> l_vec(size);
    TridMatrix matrix(size);
    //
    auto a_first = get_first_a();
    auto a_last = get_last_a();

    matrix.set(0, 0, a_first.first);
    matrix.set(0, 1, a_first.second);

    for (int i = 1; i < size - 1; ++i) {
      matrix.set(i, i - 1, get_a(i, i - 1));
      matrix.set(i, i, get_a(i, i));
      matrix.set(i, i + 1, get_a(i, i + 1));
    }

    matrix.set(size - 1, size - 2, a_last.first);
    matrix.set(size - 1, size - 1, a_last.second);

    for (int i = 0; i < size; ++i) {
      l_vec[i] = get_l(i);
    }

    return matrix.trid_matrix_alg(l_vec);
  };
};

template<typename T>
class PQCalc : public BaseCalc<T> {
private:
  function p;
  function f;
  double q;
public:
  PQCalc(const T& _basis, function _p, function _f, double _q) : BaseCalc(_basis), p(_p), q(_q), f(_f) {};
  PQCalc(T&& _basis, function _p, function _f, double _q) : BaseCalc(std::move(_basis)), p(_p), q(_q), f(_f) {};

  double get_a(int i, int j) override {
    auto [first, last] = basis.get_range(i);
    auto first_part = [&](double x) {
      return p(x) * (basis.get_deriv(i, x) * basis.get_deriv(j, x));
    };

    auto second_part = [&](double x) {
      return q * (basis.get_value(i, x) * basis.get_value(j, x));
    };
    
    return integrate(first_part, first, last, 100) + integrate(second_part, first, last, 100);
  }

  std::pair<double, double> get_first_a() override {
    return std::make_pair(get_a(0, 0) + p(0), get_a(0, 1));
  }

  std::pair<double, double> get_last_a() override {
    int size = basis.size();
    return std::make_pair(get_a(size -1, size - 2), get_a(size - 1, size - 1));
  }

  double get_l(int i) override {
    auto [first, last] = basis.get_range(i);
    auto tmp_la = [=](double x) {
      return f(x) * basis.get_value(i, x);
    };

    return integrate(tmp_la, first, last, 100);
  }
};



