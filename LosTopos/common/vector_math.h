#ifndef VECTOR_MATH_H
#define VECTOR_MATH_H

#include <algorithm>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

template <class T>
T sum(std::vector<T>& data) {
  T result = 0;
  for (unsigned int i = 0; i < data.size(); ++i) result += data[i];
  return result;
}

template <class T>
void copy(const std::vector<T>& src, std::vector<T>& dest) {
  std::copy(src.begin(), src.end(), dest.begin());
}

template <class T>
T dot(const std::vector<T>& a, const std::vector<T>& b) {
  return Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>(a.data(),
                                                               a.size())
      .dot(Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>(b.data(),
                                                                 a.size()));
}

template <class T>
void scale(T factor, std::vector<T>& data) {
  Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(data.data(), data.size()) *=
      factor;
}

template <class T>
void add_scaled(T alpha, const std::vector<T>& x,
                std::vector<T>& y) {  // y = y + alpha*x
  Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(y.data(), x.size()) +=
      Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>(x.data(),
                                                            x.size()) *
      alpha;
}

#endif
