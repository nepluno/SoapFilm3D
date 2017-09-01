#pragma once

#include <cmath>

#include "fmmtl/numeric/Complex.hpp"

/** TODO: Documentation
 */
template <typename point_type,
          typename multipole_type,
          typename local_type>
struct SphericalMultipole3D {
  //! real_type == either float or double
  typedef typename point_type::value_type  real_type;
  //! complex_type
  typedef std::complex<real_type>          complex_type;

  //! (-1)^n
  inline static constexpr real_type neg1pow(int n) {
    return ((n & 1) ? -1 : 1);
  }

  /** Kernel S2M operation
   * M += Op(s) * c where M is the multipole and s is the source
   */
  template <typename charge_type>
  inline static
  void S2M(int P, const point_type& translation, const charge_type& charge,
           multipole_type& M) {
    real_type rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
//    complex_type Z[P*(P+1)/2];   // Avoid initialization?
    complex_type * Z = (complex_type *)alloca(sizeof (complex_type) * (P*(P+1)/2));
    evalZ(rho, theta, phi, P, Z);
    int nm = 0;   // n*(n+1)/2+m
    for (int n = 0; n < P; ++n) {
      for (int m = 0; m <= n; ++m, ++nm) {
        M[nm] += neg1pow(m) * conj(Z[nm]) * charge;
      }
    }
  }

  /** Kernel M2M operator
   * M_t += Op(M_s) where M_t is the target and M_s is the source
   */
  inline static
  void M2M(int P,
           const multipole_type& Msource,
           multipole_type& Mtarget,
           const point_type& translation) {
    real_type rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
//    complex_type Z[P*(P+1)/2];
    complex_type * Z = (complex_type *)alloca(sizeof (complex_type) * (P*(P+1)/2));
    evalZ(rho, theta, phi, P, Z);
    int nm = 0;   // n*(n+1)/2+m
    for (int n = 0; n != P; ++n) {
      for (int m = 0; m <= n; ++m, ++nm) {
        auto& M = Mtarget[nm];

        for (int j = 0; j <= n; ++j) {
          // Compute the offset for Y_j
          const auto Zj = Z + j*(j+1)/2;
          // Compute the offset for M_{n-j}
          const auto Mnj = Msource.begin() + (n-j)*(n-j+1)/2;

          // All k with -j <= k <= 0 and 0 <= m-k <= n-j
          // Thus, k >= -j and k >= -n+j+m
          int k = std::max(-j, -n+j+m);
          // Thus, k <= 0 and k <= m
          for ( ; k <= 0; ++k) {
            // k is negative and m-k is positive
            M += Zj[-k] * Mnj[m-k];
          }

          // All k with 0 < k <= j and 0 <= m-k <= n-j
          // Thus, k <= j and k <= m
          int end = std::min(j, m);
          for ( ; k <= end; ++k) {
            // k is positive and m-k is positive
            M += neg1pow(k) * conj(Zj[k]) * Mnj[m-k];
          }

          // All k with 0 <= k < j and -(n-j) <= m-k <= 0
          // Thus, k <= j and k <= n+m-j
          end = std::min(j, n+m-j);
          for ( ; k <= end; ++k) {
            // k is positive and m-k is negative
            M += conj(neg1pow(m) * Zj[k] * Mnj[k-m]);
          }
        }
      }
    }
  }

  /** Kernel M2L operation
   * L += Op(M)
   *
   * @param[in] Msource The multpole expansion source
   * @param[in,out] Ltarget The local expansion target
   * @param[in] translation The vector from source to target
   * @pre translation obeys the multipole-acceptance criteria
   * @pre Msource includes the influence of all points within its box
   */
  inline static
  void M2L(int P,
           const multipole_type& Msource,
           local_type& Ltarget,
           const point_type& translation) {
    real_type rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
//    complex_type W[P*(2*P+1)];
    complex_type * W = (complex_type *)alloca(sizeof (complex_type) * (P*(2*P+1)));
    evalW(rho, theta, phi, 2*P, W);
    int nm = 0;    // n*(n+1)/2 + m
    for (int n = 0; n != P; ++n) {
      for (int m = 0; m <= n; ++m, ++nm) {
        auto& L = Ltarget[nm];

        for (int j = 0; j != P; ++j) {
          // Compute the offset for M_j
          const auto Mj = Msource.begin() + j*(j+1)/2;
          // Compute the offset for W_{j+n}
          const auto Wjn = W + (j+n)*(j+n+1)/2;

          // All k with -j <= k <= 0 and -(j+n) <= k-m <= 0
          // Thus, k >= -j and k >= m-n-j
          int k = -j;
          // Thus, k <= 0 and k <= m
          for ( ; k <= 0; ++k) {
            // k is negative and k-m is negative
            L += conj(neg1pow(m) * Wjn[m-k] * Mj[-k]);
          }

          // All k with 0 <= k <= j and -(j+n) <= k-m <= 0
          // Thus, k <= j and k <= m
          int end = std::min(j, m);
          for ( ; k <= end; ++k) {
            // k is positive and k-m is negative
            L += neg1pow(k-m) * conj(Wjn[m-k]) * Mj[k];
          }

          // All k with 0 <= k <= j and 0 <= k-m <= j+n
          // Thus, k <= j and k <= m+n+j
          for ( ; k <= j; ++k) {
            // k is positive and k-m is positive
            L += Wjn[k-m] * Mj[k];
          }
        }
      }
    }
  }

  /** Kernel L2L operator
   * L_t += Op(L_s) where L_t is the target and L_s is the source
   *
   * @param[in] source The local source at the parent level
   * @param[in,out] target The local target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Lsource includes the influence of all points outside its box
   */
  inline static
  void L2L(int P,
           const local_type& Lsource,
           local_type& Ltarget,
           const point_type& translation) {
    real_type rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
//    complex_type Z[P*(P+1)/2];
    complex_type * Z = (complex_type *)alloca(sizeof (complex_type) * (P*(P+1)/2));
    evalZ(rho, theta, phi, P, Z);
    int nm = 0;    // n*(n+1)/2 + m
    for (int n = 0; n != P; ++n) {
      for (int m = 0; m <= n; ++m, ++nm) {
        auto& L = Ltarget[nm];

        for (int j = n; j != P; ++j) {
          // Compute the offset for L_j
          const auto Lj = Lsource.begin() + j*(j+1)/2;
          // Compute the offset for Z_{j-n}
          const auto Zjn = Z + (j-n)*(j-n+1)/2;

          // All k with -j <= k <= 0 and -(j-n) <= k-m <= 0
          // Thus, k >= -j and k >= n+m-j
          int k = n+m-j;
          // Thus, k <= 0 and k <= m
          for ( ; k <= 0; ++k) {
            // k is negative and k-m is negative
            L += conj(neg1pow(m) * Zjn[m-k] * Lj[-k]);
          }

          // All k with 0 <= k <= j and -(j-n) <= k-m <= 0
          // Thus, k <= j and k <= m
          int end = std::min(j, m);
          for ( ; k <= end; ++k) {
            // k is positive and k-m is negative
            L += neg1pow(k-m) * conj(Zjn[m-k]) * Lj[k];
          }

          // All k with 0 <= k <= j and 0 <= k-m <= j-n
          // Thus, k <= j and k <= m-n+j
          end = std::min(j, m-n+j);
          for ( ; k <= end; ++k) {
            // k is positive and k-m is positive
            L += Zjn[k-m] * Lj[k];
          }
        }
      }
    }
  }


  /** Spherical to cartesian coordinates */
  inline static
  point_type sph2cart(real_type rho, real_type theta, real_type phi,
                      const point_type& s) {
    using std::cos;
    using std::sin;
    const real_type st = sin(theta);
    const real_type ct = cos(theta);
    const real_type sp = sin(phi);
    const real_type cp = cos(phi);
    return point_type(s[0]*st*cp + s[1]*ct*cp/rho - s[2]*sp/(rho*st),
                      s[0]*st*sp + s[1]*ct*sp/rho + s[2]*cp/(rho*st),
                      s[0]*ct    - s[1]*st/rho);
  }

  /** Cartesian to spherical coordinates */
  inline static
  void cart2sph(real_type& r, real_type& theta, real_type& phi,
                const point_type& x) {
    using std::acos;
    using std::atan2;
    r = norm_2(x);
    theta = acos(x[2] / (r + 1e-100));
    phi = atan2(x[1], x[0]);
  }

  /** Computes the function
   * Z[n*(n+1)/2+m]
   *        = Z_n^m
   *        = i^-|m| A_n^m (-1)^n rho^n Y_n^m(theta, phi)
   *        = i^-|m| (-1)^n/sqrt((n+m)!(n-m)!) (-1)^n rho^n Y_n^m(theta, phi)
   *        = rho^n/(n+|m|)! P_n^|m|(cos theta) exp(i m phi) i^-|m|
   * for all 0 <= n < P and all 0 <= m <= n.
   *
   * @param[in] rho    The radial distance, rho > 0.
   * @param[in] theta  The inclination or polar angle, 0 <= theta <= pi.
   * @param[in] phi    The azimuthal angle, 0 <= phi <= 2pi.
   * @param[in] P      The maximum degree spherical harmonic to compute, P > 0.
   * @param[out] Y,dY  The output arrays to store Znm and its theta-derivative.
   *                     Each has storage for P*(P+1)/2 elements.
   *
   * Note this uses the spherical harmonics definition
   * Y_n^m(\theta, \phi) =
   *        \sqrt{\frac{(n-|m|)!}{(n+|m|)!}} P_n^{|m|}(\cos \theta) e^{i m \phi)
   * which has symmetries
   *   Y_n^m(\pi - \theta, \pi + \phi) = (-1)^n Y_n^m(\theta, \phi)
   *   Y_n^{-m} = (-1)^m conj(Y_n^m)
   *
   * Note the symmetry in the result of this function
   *    Z_n^{-m} = (-1)^m conj(Z_n^m)
   * where conj denotes complex conjugation.
   *
   * These are not the spherical harmonics, but are the spherical
   * harmonics with the prefactor (often denoted A_n^m) included. These are useful
   * for computing multipole and local expansions in an FMM.
   */
  inline static
  void evalZ(real_type rho, real_type theta, real_type phi, int P,
             complex_type* Y, complex_type* dY = nullptr) {
    typedef real_type     real;
    typedef complex_type  complex;
    using std::cos;
    using std::sin;

    const real    ct = cos(theta);
    const real    st = sin(theta);
    const complex ei = complex(sin(phi),-cos(phi)); // exp(i*phi) i^-1

    int m = 0;
    int nm;
    real    Pmm = 1;                                // Init Legendre Pmm(ct)
    real   rhom = 1;                                // Init rho^n / (n+m)!
    complex eim = 1;                                // Init exp(i*m*phi) i^-m
    while (true) {                                  // For all 0 <= m < P
      // n == m
      nm = m*(m+1)/2 + m;                           //  Index of Znm for m > 0
      Y[nm] = rhom * Pmm * eim;                     //  Znm for m > 0
      if (dY)
        dY[nm] = m*ct/st * Y[nm];                   //  Theta derivative of Znm

      // n == m+1
      int n = m + 1;
      if (n == P) return;                           //  Done! m == P-1

      real Pnm  = ct * (2*m+1) * Pmm;               //  P_{m+1}^m(x) = x(2m+1)Pmm
      real rhon = rhom * rho / (n+m);               //  rho^n / (n+m)!
      nm += n;                                      //  n*(n+1)/2 + m
      Y[nm] = rhon * Pnm * eim;                     //  Znm for m > 0
      if (dY)
        dY[nm] = (n*ct-(n+m)*Pmm/Pnm)/st * Y[nm];   //  Theta derivative of Znm

      // m+1 < n < P
      real Pn1m = Pmm;                              //  P_{n-1}^m
      while (++n != P) {
        real Pn2m = Pn1m;                           //   P_{n-2}^m
        Pn1m = Pnm;                                 //   P_{n-1}^m
        Pnm = (ct*(2*n-1)*Pn1m-(n+m-1)*Pn2m)/(n-m); //   P_n^m recurrence
        rhon *= rho / (n + m);                      //   rho^n / (n+m)!

        nm += n;                                    //   n*(n+1)/2 + m
        Y[nm] = rhon * Pnm * eim;                   //   Znm for m > 0
        if (dY)
          dY[nm] = (n*ct-(n+m)*Pn1m/Pnm)/st * Y[nm];//   Theta derivative of Znm
      }

      ++m;                                          //  Increment m

      rhom *= rho / (2*m*(2*m-1));                  //  rho^m / (2m)!
      Pmm  *= -st * (2*m-1);                        //  P_{m+1}^{m+1} recurrence
      eim  *= ei;                                   //  exp(i*m*phi) i^-m
    }                                               // End loop over m in Znm
  }


  /** Computes the function
   * W[n*(n+1)/2+m]
   *         = W_n^m
   *         = i^|m| / A_n^m rho^{-n-1} Y_n^m(theta, phi)
   *         = i^|m| (-1)^n sqrt((n+m)!(n-m)!) rho^{-n-1} Y_n^m(theta, phi)
   *         = (-1)^n rho^{-n-1} (n-|m|)! P_n^|m|(cos theta) exp(i m phi) i^|m|
   * for all 0 <= n < P and all 0 <= m <= n.
   *
   * @param[in] rho    The radial distance, rho > 0.
   * @param[in] theta  The inclination or polar angle, 0 <= theta <= pi.
   * @param[in] phi    The azimuthal angle, 0 <= phi <= 2pi.
   * @param[in] P      The maximum degree spherical harmonic to compute, P > 0.
   * @param[out] Y,dY  The output arrays to store Wnm and its theta-derivative.
   *                     Each has storage for P*(P+1)/2 elements.
   *
   * Note this uses the spherical harmonics definition
   * Y_n^m(\theta, \phi) =
   *        \sqrt{\frac{(n-|m|)!}{(n+|m|)!}} P_n^{|m|}(\cos \theta) e^{i m \phi)
   * which has symmetries
   *   Y_n^m(\pi - \theta, \pi + \phi) = (-1)^n Y_n^m(\theta, \phi)
   *   Y_n^{-m} = (-1)^m conj(Y_n^m)
   *
   * Note the symmetry in the result of this function
   *    W_n^{-m} = (-1)^m conj(W_n^m)
   * where conj denotes complex conjugation.
   *
   * These are not the spherical harmonics, but are the spherical
   * harmonics with the prefactor (often denoted A_n^m) included. These are useful
   * for computing multipole and local expansions in an FMM.
   */
  inline static
  void evalW(real_type rho, real_type theta, real_type phi, int P,
             complex_type* Y) {
    typedef real_type     real;
    typedef complex_type  complex;
    using std::cos;
    using std::sin;

    const real    ct = cos(theta);
    const real    st = sin(theta);
    const complex ei = complex(-sin(phi),cos(phi));// exp(i*phi) i

    rho = 1 / rho;
    int m = 0;
    int nm;
    real    Pmm = 1;                               // Init Legendre Pmm(ct)
    real   rhom = rho;                             // Init -1^n rho^{-n-1} (n-m)!
    complex eim = 1;                               // Init exp(i*m*phi) i^-m
    while (true) {                                 // For all 0 <= m < P
      // n == m
      nm = m*(m+1)/2 + m;                          //  Index of Wnm for m > 0
      Y[nm] = rhom * Pmm * eim;                    //  Wnm for m > 0

      // n == m+1
      int n = m+1;
      if (n == P) return;                          //  Done! m == P-1

      real Pnm  = ct * (2*m+1) * Pmm;              //  P_{m+1}^m(x) = x(2m+1)Pmm
      real rhon = rhom * -rho;                     //  -1^n rho^{-n-1} (n-m)!
      nm += n;                                     //  n*(n+1)/2 + m
      Y[nm] = rhon * Pnm * eim;                    //  Wnm for m > 0

      // m+1 < n < P
      real Pn1m = Pmm;                              //  P_{n-1}^m
      while (++n != P) {
        real Pn2m = Pn1m;                           //   P_{n-2}^m
        Pn1m = Pnm;                                 //   P_{n-1}^m
        Pnm = (ct*(2*n-1)*Pn1m-(n+m-1)*Pn2m)/(n-m); //   P_n^m recurrence
        rhon *= -rho * (n - m);                     //   -1^n rho^{-n-1} (n-m)!

        nm += n;                                    //   n*(n+1)/2 + m
        Y[nm] = rhon * Pnm * eim;                   //   Wnm for m > 0
      }

      ++m;                                          //  Increment m

      rhom *= -rho;                                 //  -1^m rho^{-m-1} (n-m)!
      Pmm  *= -st * (2*m-1);                        //  P_{m+1}^{m+1} recurrence
      eim  *= ei;                                   //  exp(i*m*phi) i^m
    }                                               // End loop over m in Wnm
  }

};
