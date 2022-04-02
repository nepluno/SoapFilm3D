#include <Eigen/Core>
#include <Eigen/Dense>

#include <dense_matrix.h>

#include <cstring>

void DenseMatrix::clear(void) {
  m = n = 0;
  value.clear();
}

void DenseMatrix::set_zero(void) {
  Eigen::Map<Eigen::VectorXd>(value.data(), value.size()).setZero();
}

void DenseMatrix::resize(int m_, int n_) {
  m = m_;
  n = n_;
  value.resize(m * n, 0);
}

void DenseMatrix::apply(const double *x, double *y) const {
  assert(x && y);
  Eigen::Map<Eigen::VectorXd>(y, m) =
      Eigen::Map<const Eigen::MatrixXd>(value.data(), m, n) *
      Eigen::Map<const Eigen::VectorXd>(x, n);
}

void DenseMatrix::apply_and_subtract(const double *x, const double *y,
                                     double *z) const {
  assert(x && y);
  Eigen::Map<Eigen::VectorXd>(z, m) =
      Eigen::Map<const Eigen::VectorXd>(y, m) -
      Eigen::Map<const Eigen::MatrixXd>(value.data(), m, n) *
          Eigen::Map<const Eigen::VectorXd>(x, n);
}

void DenseMatrix::apply_transpose(const double *x, double *y) const {
  assert(x && y);
  Eigen::Map<Eigen::VectorXd>(y, n) =
      Eigen::Map<const Eigen::MatrixXd>(value.data(), m, n).transpose() *
      Eigen::Map<const Eigen::VectorXd>(x, m);
}

void DenseMatrix::apply_transpose_and_subtract(const double *x, const double *y,
                                               double *z) const {
  assert(x && y);
  Eigen::Map<Eigen::VectorXd>(z, n) =
      Eigen::Map<const Eigen::VectorXd>(y, n) -
      Eigen::Map<const Eigen::MatrixXd>(value.data(), m, n).transpose() *
          Eigen::Map<const Eigen::VectorXd>(x, m);
}

void DenseMatrix::write_matlab(std::ostream &output,
                               const char *variable_name) const {
  output << variable_name << "=[";
  std::streamsize old_precision = output.precision();
  output.precision(18);
  for (int i = 0; i < m; ++i) {
    if (i > 0) output << " ";
    for (int j = 0; j < n - 1; ++j) output << value[i + j * m] << " ";
    output << value[i + (n - 1) * m];
    if (i < m - 1)
      output << std::endl;
    else
      output << "];" << std::endl;
  }
  output.precision(old_precision);
}

void transpose(const DenseMatrix &A, DenseMatrix &Atranspose) {
  Atranspose.resize(A.n, A.m);
  for (int j = 0; j < A.n; ++j)
    for (int i = 0; i < A.m; ++i) {
      Atranspose(j, i) = A(i, j);
    }
}

void multiply(const DenseMatrix &A, const DenseMatrix &B, DenseMatrix &C) {
  assert(A.n == B.m);
  C.resize(A.m, B.n);
  Eigen::Map<Eigen::MatrixXd>(C.value.data(), A.m, B.n) =
      Eigen::Map<const Eigen::MatrixXd>(A.value.data(), A.m, A.n) *
      Eigen::Map<const Eigen::MatrixXd>(B.value.data(), B.m, B.n);
}

void multiply_with_transpose(const DenseMatrix &A, DenseMatrix &ATA) {
  ATA.resize(A.n, A.n);
  Eigen::Map<Eigen::MatrixXd>(ATA.value.data(), A.n, A.n) =
      Eigen::Map<const Eigen::MatrixXd>(A.value.data(), A.m, A.n).transpose() *
      Eigen::Map<const Eigen::MatrixXd>(A.value.data(), A.m, A.n);
}
