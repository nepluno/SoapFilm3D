#pragma once
/** @file BarycentricTaylor.hpp
 * @brief Implements the Barycentric kernel with cartesian Taylor expansions.
 *
 * K(t,s) = 1 / s-t
 */

#include <numeric>

#include "fmmtl/Expansion.hpp"

#include "Barycentric.kern"
#include "Util/GradedPolynomial.hpp"


/** BarycentricTaylor
 * @tparam P  The order of the polynomial expansion
 */
template <unsigned P>
class BarycentricTaylor
    : public fmmtl::Expansion<Barycentric, BarycentricTaylor<P> >
{
 public:
  FMMTL_IMPORT_KERNEL_TRAITS(Barycentric);

  //! Point type
  typedef Vec<1,double> point_type;

  //! Multipole expansion type
  typedef fmmtl::GradedPolynomial<double,1,P> multipole_type;
  //! Local expansion type
  typedef fmmtl::GradedPolynomial<double,1,P> local_type;

  //! A multiindex ranging from 0 to P
  static constexpr fmmtl::multiindex<P> k{};


  void init_multipole(multipole_type& M, const point_type&, unsigned) const {
    M.fill(0);
  }
  void init_local(local_type& L, const point_type&, unsigned) const {
    L.fill(0);
  }


  // Define the S2M
  void S2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    //Mt += pow(center-source, k) * charge;

    // For all k < P,  M[k] = pow(r,k)/k! * charge
    double r = center[0] - source[0];
    double rk = charge;               // rk = r^k / k! * charge
    M[0] += rk;                       // k == 0
    for (unsigned k = 1; k <= P; ++k)
      M[k] += (rk *= r / k);
  }

  // Define the M2M
  void M2M(const multipole_type& Ms,
                 multipole_type& Mt,
           const point_type& translation) const {
    Mt += pow(translation, k) * Ms;

    /*
    double r = translation[0];

    Mt[0] += Ms[0];
    for (unsigned n = 1; n <= P; ++n) {
      double rk = 1;                       // rk = r^k / k!
      double& Mtn = (Mt[n] += Ms[n]);      // k == 0
      for (unsigned k = 1; k <= n; ++k)
        Mtn += (rk *= r / k) * Ms[n-k];
    }
    */
  }

  // Define the M2L
  void M2L(const multipole_type& M,
                 local_type& L,
           const point_type& translation) const {
    /*
    // Precompute the factors (-1)^k k! / r^{k+1}
    fmmtl::GradedPolynomial<double,1,2*P> tmp(-translation, translation[0]);
    for (unsigned k = 0; k <= 2*P; ++k)  tmp[k] = 1.0 / tmp[k];

    // Multiply
    for (unsigned n = 0; n <= P; ++n)
      L[n] = std::inner_product(M.begin(), M.end(), tmp.begin()+n, L[n]);
    */

    // L[n] += sum_{k < P} (-1)^{n+k} (n+k)! / r^{n+k+1} * M[k]
    double r   = 1.0 / translation[0];
    double rn  = r;        // r^{n+1}     n!   (-1)^(n)
    double rnk = rn;       // r^{n+k+1} (n+k)! (-1)^(n+k)

    double& Ln = (L[0] += rnk * M[0]);     // n == 0, k == 0
    for (unsigned k = 1; k <= P; ++k)
      Ln += (rnk *= -r * k) * M[k];

    for (unsigned n = 1; n <= P; ++n) {
      rnk = (rn *= -r * n);
      double& Ln = (L[n] += rnk * M[0]);   // k == 0
      for (unsigned k = 1; k <= P; ++k)
        Ln += (rnk *= -r * (n+k)) * M[k];
    }
  }

  // Define the L2L
  void L2L(const local_type& Ls,
                 local_type& Lt,
           const point_type& translation) const {
    Lt += pow(translation, k) / Ls;

    /*
    double r = translation[0];

    // For all n < P, Lt[n] += sum_{k >= n} 1/(k-n)! r^{k-n} Ls[k]
    // TODO: Unroll & Compile-time
    for (unsigned n = 0; n < P; ++n) {
      double rk = 1;                       // rk = r^k / k!
      double& Ltn = (Lt[n] += Ls[n]);      // k == 0
      for (unsigned k = 1; k <= P-n; ++k)
        Ltn += (rk *= r / k) * Ls[n+k];
    }
    Lt[P] += Ls[P];
    */
  }

  // Define the L2T
  void L2T(const local_type& L, const point_type& center,
           const target_type& target, result_type& result) const {
    // result += inner_prod(pow(target-center,k), L);

    // result += sum_k pow(r,k)/k! * L[k]
    double r = target[0] - center[0];
    double rk = 1;                         // rk = r^k / k!
    result += L[0];                        // k == 0
    for (unsigned k = 1; k <= P; ++k)
      result += (rk *= r / k) * L[k];
  }
};
