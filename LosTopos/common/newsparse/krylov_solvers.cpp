#include <Eigen/Core>
#include <Eigen/Dense>
#include <krylov_solvers.h>

#include <cassert>

//============================================================================
KrylovSolverStatus CG_Solver::solve(const LinearOperator &A, const double *rhs,
                                    double *result,
                                    const LinearOperator *preconditioner,
                                    bool use_given_initial_guess) {
  const int n = A.m;
  assert(A.n == n);
  assert(preconditioner == 0 ||
         (preconditioner->m == n && preconditioner->n == n));
  if ((int)s.size() != n) {
    r.resize(n);
    z.resize(n);
    s.resize(n);
  }
  // convergence tolerance
  double tol =
      tolerance_factor *
      Eigen::Map<const Eigen::VectorXd>(rhs, n).lpNorm<Eigen::Infinity>();
  // initial guess
  if (use_given_initial_guess) {
    A.apply_and_subtract(result, rhs, &r[0]);
  } else {
    Eigen::Map<Eigen::VectorXd>(result, n).setZero();
    Eigen::Map<Eigen::VectorXd>(r.data(), n) =
        Eigen::Map<const Eigen::VectorXd>(rhs, n);
  }
  // check instant convergence
  iteration = 0;
  residual_norm =
      Eigen::Map<const Eigen::VectorXd>(r.data(), n).lpNorm<Eigen::Infinity>();
  if (residual_norm == 0) return status = KRYLOV_CONVERGED;
  // set up CG
  double rho;
  if (preconditioner)
    preconditioner->apply(r, z);
  else
    z = r;
  rho = Eigen::Map<const Eigen::VectorXd>(r.data(), n)
            .dot(Eigen::Map<const Eigen::VectorXd>(z.data(), n));
  if (rho <= 0 || rho != rho) return status = KRYLOV_BREAKDOWN;
  s = z;
  // and iterate
  for (iteration = 1; iteration < max_iterations; ++iteration) {
    double alpha;
    A.apply(s, z);  // reusing z=A*s
    double sz = Eigen::Map<const Eigen::VectorXd>(s.data(), n)
                    .dot(Eigen::Map<const Eigen::VectorXd>(z.data(), n));
    if (sz <= 0 || sz != sz) return status = KRYLOV_BREAKDOWN;
    alpha = rho / sz;
    Eigen::Map<Eigen::VectorXd>(result, n) +=
        Eigen::Map<const Eigen::VectorXd>(s.data(), n) * alpha;
    Eigen::Map<Eigen::VectorXd>(r.data(), n) +=
        -Eigen::Map<const Eigen::VectorXd>(z.data(), n) * alpha;
    residual_norm = Eigen::Map<const Eigen::VectorXd>(r.data(), n)
                        .lpNorm<Eigen::Infinity>();
    if (residual_norm <= tol) return status = KRYLOV_CONVERGED;
    if (preconditioner)
      preconditioner->apply(r, z);
    else
      z = r;
    double rho_new = Eigen::Map<const Eigen::VectorXd>(r.data(), n)
                         .dot(Eigen::Map<const Eigen::VectorXd>(z.data(), n));
    if (rho_new <= 0 || rho_new != rho_new) return status = KRYLOV_BREAKDOWN;
    double beta = rho_new / rho;
    Eigen::Map<Eigen::VectorXd>(z.data(), n) +=
        Eigen::Map<const Eigen::VectorXd>(s.data(), n) * beta;
    s.swap(z);  // s=beta*s+z
    rho = rho_new;
  }
  return status = KRYLOV_EXCEEDED_MAX_ITERATIONS;
}

//============================================================================
KrylovSolverStatus MINRES_CR_Solver::solve(const LinearOperator &A,
                                           const double *rhs, double *result,
                                           const LinearOperator *preconditioner,
                                           bool use_given_initial_guess) {
  const int n = A.m;
  assert(A.n == n);
  assert(preconditioner == 0 ||
         (preconditioner->m == n && preconditioner->n == n));
  if ((int)s.size() != n) {
    r.resize(n);
    z.resize(n);
    q.resize(n);
    s.resize(n);
    t.resize(n);
  }
  // convergence tolerance
  double tol =
      tolerance_factor *
      Eigen::Map<const Eigen::VectorXd>(rhs, n).lpNorm<Eigen::Infinity>();
  // initial guess
  if (use_given_initial_guess) {
    A.apply_and_subtract(result, rhs, &r[0]);
  } else {
    Eigen::Map<Eigen::VectorXd>(result, n).setZero();
    Eigen::Map<Eigen::VectorXd>(r.data(), n) =
        Eigen::Map<const Eigen::VectorXd>(rhs, n);
  }
  // check instant convergence
  iteration = 0;
  residual_norm =
      Eigen::Map<const Eigen::VectorXd>(r.data(), n).lpNorm<Eigen::Infinity>();
  if (residual_norm == 0) return status = KRYLOV_CONVERGED;
  // set up CR
  double rho;
  if (preconditioner)
    preconditioner->apply(r, z);
  else
    s = r;
  A.apply(s, t);
  rho = Eigen::Map<const Eigen::VectorXd>(r.data(), n)
            .dot(Eigen::Map<const Eigen::VectorXd>(t.data(), n));
  if (rho == 0 || rho != rho) return status = KRYLOV_BREAKDOWN;
  // and iterate
  for (iteration = 1; iteration < max_iterations; ++iteration) {
    double alpha;
    double tt = Eigen::Map<const Eigen::VectorXd>(t.data(), n)
                    .dot(Eigen::Map<const Eigen::VectorXd>(t.data(), n));
    if (tt == 0 || tt != tt) return status = KRYLOV_BREAKDOWN;
    alpha = rho / tt;
    Eigen::Map<Eigen::VectorXd>(result, n) +=
        Eigen::Map<const Eigen::VectorXd>(s.data(), n) * alpha;
    Eigen::Map<Eigen::VectorXd>(r.data(), n) +=
        -Eigen::Map<const Eigen::VectorXd>(t.data(), n) * alpha;

    residual_norm = Eigen::Map<const Eigen::VectorXd>(r.data(), n)
                        .lpNorm<Eigen::Infinity>();
    if (residual_norm <= tol) return KRYLOV_CONVERGED;
    if (preconditioner)
      preconditioner->apply(r, z);
    else
      z = r;
    A.apply(z, q);
    double rho_new = Eigen::Map<const Eigen::VectorXd>(r.data(), n)
                         .dot(Eigen::Map<const Eigen::VectorXd>(q.data(), n));
    if (rho_new == 0 || rho_new != rho_new) return KRYLOV_BREAKDOWN;
    double beta = rho_new / rho;
    Eigen::Map<Eigen::VectorXd>(z.data(), n) +=
        Eigen::Map<const Eigen::VectorXd>(s.data(), n) * beta;
    s.swap(z);  // s=beta*s+z
    Eigen::Map<Eigen::VectorXd>(q.data(), n) +=
        Eigen::Map<const Eigen::VectorXd>(t.data(), n) * beta;
    t.swap(q);  // t=beta*t+q
    rho = rho_new;
  }
  return KRYLOV_EXCEEDED_MAX_ITERATIONS;
}

//============================================================================
KrylovSolverStatus CGNR_Solver::solve(const LinearOperator &A,
                                      const double *rhs, double *result,
                                      const LinearOperator *preconditioner,
                                      bool use_given_initial_guess) {
  const int m = A.m, n = A.n;
  assert(preconditioner == 0 ||
         (preconditioner->m == n && preconditioner->n == n));
  if ((int)s.size() != n) {
    r.resize(n);
    z.resize(n);
    s.resize(n);
    u.resize(m);
  }
  // convergence tolerance
  A.apply_transpose(rhs, &r[0]);  // form A^T*rhs in r
  double tol =
      tolerance_factor *
      Eigen::Map<const Eigen::VectorXd>(r.data(), n).lpNorm<Eigen::Infinity>();
  // initial guess
  if (use_given_initial_guess) {
    A.apply_and_subtract(result, rhs, &u[0]);
    A.apply_transpose(u, r);
  } else {
    Eigen::Map<Eigen::VectorXd>(result, n).setZero();
  }
  // check instant convergence
  iteration = 0;
  residual_norm =
      Eigen::Map<const Eigen::VectorXd>(r.data(), n).lpNorm<Eigen::Infinity>();
  if (residual_norm == 0) return status = KRYLOV_CONVERGED;
  // set up CG
  double rho;
  if (preconditioner)
    preconditioner->apply(r, z);
  else
    z = r;
  rho = Eigen::Map<const Eigen::VectorXd>(r.data(), n)
            .dot(Eigen::Map<const Eigen::VectorXd>(z.data(), n));
  if (rho <= 0 || rho != rho) return status = KRYLOV_BREAKDOWN;
  s = z;
  // and iterate
  for (iteration = 1; iteration < max_iterations; ++iteration) {
    double alpha;
    A.apply(s, u);
    A.apply_transpose(u, z);
    double sz = Eigen::Map<const Eigen::VectorXd>(u.data(), n)
                    .dot(Eigen::Map<const Eigen::VectorXd>(u.data(), n));
    if (sz <= 0 || sz != sz) return status = KRYLOV_BREAKDOWN;
    alpha = rho / sz;
    Eigen::Map<Eigen::VectorXd>(result, n) +=
        Eigen::Map<const Eigen::VectorXd>(s.data(), n) * alpha;
    Eigen::Map<Eigen::VectorXd>(r.data(), n) +=
        -Eigen::Map<const Eigen::VectorXd>(z.data(), n) * alpha;
    residual_norm = Eigen::Map<const Eigen::VectorXd>(r.data(), n)
                        .lpNorm<Eigen::Infinity>();
    if (residual_norm <= tol) return status = KRYLOV_CONVERGED;
    if (preconditioner)
      preconditioner->apply(r, z);
    else
      z = r;
    double rho_new = Eigen::Map<const Eigen::VectorXd>(r.data(), n)
                         .dot(Eigen::Map<const Eigen::VectorXd>(z.data(), n));
    if (rho_new <= 0 || rho_new != rho_new) return status = KRYLOV_BREAKDOWN;
    double beta = rho_new / rho;
    Eigen::Map<Eigen::VectorXd>(z.data(), n) +=
        Eigen::Map<const Eigen::VectorXd>(s.data(), n) * beta;
    s.swap(z);  // s=beta*s+z
    rho = rho_new;

    //      if ( iteration % 5000 == 0 )
    //      {
    //         std::cout << "CGNR_Solver --- residual_norm: " << residual_norm
    //         << std::endl;
    //      }
  }
  return status = KRYLOV_EXCEEDED_MAX_ITERATIONS;
}
