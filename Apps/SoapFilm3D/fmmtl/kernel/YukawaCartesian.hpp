#pragma once
/** @file YukawaCartesian.hpp
 * @brief Implements the Yukawa kernel with cartesian expansions
 *
 * K(t,s) = exp(-Kappa*|t-s|) / |t-s|                           // Potential
 * K(t,s) = -(Kappa*|t-s|+1) exp(-Kappa*|t-s|) (t-s) / |t-s|^3  // Force
 */

#include <complex>
#include <vector>
#include <cassert>

#include "Yukawa.kern"

#include "fmmtl/Expansion.hpp"
#include "fmmtl/numeric/Vec.hpp"

class YukawaCartesian
    : public fmmtl::Expansion<YukawaKernel, YukawaCartesian> {
 protected:
  typedef double real;
  typedef std::complex<real> complex;
  typedef std::vector<real> real_vec;
  struct IndexCache;

  //! Expansion order
  const int P;
  //! number of multipole terms
  unsigned MTERMS;

  //! I, J, K arrays
  std::vector<unsigned> I, J, K;
  //! indices of multipole terms
  std::vector<unsigned> index;
  //! factorial cache
  std::vector<double> fact;

 protected:
  //! Store all possible index combinations returned from setIndex
  // TODO: Static/Analytic computation
  struct IndexCache {
   private:
    unsigned setIndex(unsigned i, unsigned j, unsigned k) const {
      unsigned II = 0;
      for (unsigned ii = 0; ii < i; ++ii)
        for (unsigned jj = 1; jj < P_+2-ii; ++jj)
          II += jj;
      for (unsigned jj = P_+2-j; jj < P_+2; ++jj)
        II += jj-i;
      return II + k;
    }

    // use 1D vector & index into it
    std::vector<unsigned> indices;
    unsigned P_;

   public:
    IndexCache(unsigned P)
        : P_(P) {
      indices = std::vector<unsigned>((P_+1)*(P_+1)*(P_+1),0);
      for (unsigned i=0; i<P_+1; i++) {
        for (unsigned j=0; j<P_+1-i; j++) {
          for (unsigned k=0; k<P_+1-i-j; k++) {
            indices[i*(P_+1)*(P_+1)+j*(P_+1) + k] = setIndex(i,j,k);
          }
        }
      }
    }

    inline unsigned operator()(unsigned i, unsigned j, unsigned k) const {
      return indices[i*(P_+1)*(P_+1)+j*(P_+1)+k];
    }
  };

  IndexCache index_cache;

 public:
  //! Point type to use for the trees
  typedef Vec<3,real> point_type;

  //! Multipole expansion type
  typedef std::vector<real> multipole_type;
  //! Local expansion type
  typedef std::vector<real> local_type;

  //! Default constructor - use delegating constructor
  YukawaCartesian()
      : YukawaCartesian(4,0.125) {
  }
  //! Constructor
  YukawaCartesian(int p, double _kappa)
      : Expansion(YukawaKernel(_kappa)),
        P(p), MTERMS((P+1)*(P+2)*(P+3)/6), fact(2*P), index_cache(P) {
    //kappa = _kappa;  // Sets the YukawaKernel kappa

    I = std::vector<unsigned>(MTERMS,0);
    J = std::vector<unsigned>(MTERMS,0);
    K = std::vector<unsigned>(MTERMS,0);
    index = std::vector<unsigned>(MTERMS,0);

    // generate n -> (i,j,k) arrays
    unsigned idx=0;
    for (int ii=0; ii<P+1; ii++) {
      for (int jj=0; jj<P+1-ii; jj++) {
        for (int kk=0; kk<P+1-ii-jj; kk++) {
          index[idx] = index_cache(ii,jj,kk);
          I[idx] = ii;
          J[idx] = jj;
          K[idx] = kk;
          idx++;
        }
      }
    }

    // generate factorials
    fact[0] = 1;
    for (int i=1; i<2*P; i++) fact[i] = i*fact[i-1];
  }

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M, const point_type&, unsigned) const {
    M = std::vector<real>((P+1)*(P+2)*(P+3)/6, 0);
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L, const point_type&, unsigned) const {
    L = std::vector<real>((P+1)*(P+2)*(P+3)/6, 0);
  }

  /** Kernel S2M operation
   * M += Op(s) * c where M is the multipole and s is the source
   *
   * @param[in] source The source to accumulate into the multipole
   * @param[in] charge The source's corresponding charge
   * @param[in] center The center of the box containing the multipole expansion
   * @param[in,out] M The multipole expansion to accumulate into
   */
  void S2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    point_type dX = center - source;

    for (unsigned i = 0; i < MTERMS; ++i) {
      double fact_term = fact[I[i]]*fact[J[i]]*fact[K[i]];
      M[i] += charge * pow(dX[0],I[i])
                     * pow(dX[1],J[i])
                     * pow(dX[2],K[i]) / fact_term;
    }
  }

  /** Kernel M2M operator
   * M_t += Op(M_s) where M_t is the target and M_s is the source
   *
   * @param[in] source The multipole source at the child level
   * @param[in,out] target The multipole target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Msource includes the influence of all points within its box
   */
  void M2M(const multipole_type& Msource,
                 multipole_type& Mtarget,
           const point_type& translation) const {
    const unsigned MTERMS = (P+1)*(P+2)*(P+3)/6;
    const auto dX = translation;

    for (unsigned i=0; i<MTERMS; ++i) {
      unsigned n[3] = {I[i], J[i], K[i]};

      for (unsigned ii=0; ii<n[0]+1; ii++) {
        for (unsigned jj=0; jj<n[1]+1; jj++) {
          for (unsigned kk=0; kk<n[2]+1; kk++) {

            unsigned Midx = index_cache(ii,jj,kk);
            double fact_term = fact[n[0]-ii]*fact[n[1]-jj]*fact[n[2]-kk];
            Mtarget[i] += Msource[Midx]*pow(dX[0],n[0]- ii)*pow(dX[1],n[1]-jj)*pow(dX[2],n[2]-kk)/fact_term;
          }
        }
      }
    }
  }

  /** Kernel M2T operation
   * r += Op(M, t) where M is the multipole and r is the result
   *
   * @param[in] M The multpole expansion
   * @param[in] center The center of the box with the multipole expansion
   * @param[in] target The target to evaluate the multipole at
   * @param[in,out] result The target's corresponding result to accumulate into
   * @pre M includes the influence of all sources within its box
   */
  void M2T(const multipole_type& M, const point_type& center,
           const target_type& target, result_type& result) const {
    std::vector<real> a_aux(MTERMS,0), ax_aux(MTERMS,0), ay_aux(MTERMS,0), az_aux(MTERMS,0);

    point_type dX = target - center;

    getCoeff(a_aux, ax_aux, ay_aux, az_aux, dX);

    // loop over tarSize
    for (unsigned j=0; j<MTERMS; j++) {
      double fact_term = fact[I[j]]*fact[J[j]]*fact[K[j]];

      result[0] +=  a_aux[j]*M[j]*fact_term;
      result[1] += ax_aux[j]*M[j]*fact_term;
      result[2] += ay_aux[j]*M[j]*fact_term;
      result[3] += az_aux[j]*M[j]*fact_term;
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
  void M2L(const multipole_type& Msource,
                 local_type& Ltarget,
           const point_type& translation) const {
    std::vector<real> a_aux(MTERMS,0), ax_aux(MTERMS,0), ay_aux(MTERMS,0), az_aux(MTERMS,0);

    getCoeff(a_aux, ax_aux, ay_aux, az_aux, translation);

    // get rid of factorial terms from getCoeff
    for (unsigned i=0; i<MTERMS; i++) {
      a_aux[i] *= fact[I[i]]*fact[J[i]]*fact[K[i]];
    }

    for (int ik=0; ik<P+1; ik++) {
      for (int jk=0; jk<P+1-ik; jk++) {
        for (int kk=0; kk<P+1-ik-jk; kk++) {
          // k = (ik, jk, kk)
          int Lk  = index_cache(ik,jk,kk);
          //
          for (int in=0; in<P+1-ik; in++) {
            // printf("in: %d\n",in);
            for (int jn=0; jn<P+1-jk; jn++) {
              for (int kn=0; kn<P+1-kk; kn++) {

                if (in+jn+kn > P) continue;
                if (in+ik+jn+jk+kn+kk > P) continue;
                //     n = (in, jn, kn)
                // (n+k) = (ik+in, jk+jn, kk+kn)
                // L_k   = \sum_{n=0}^{p-k} a_aux[getIndex(n+k)]*M[getIndex(n)]
                //
                int Mn  = index_cache(in,jn,kn);
                int npk = index_cache(ik+in,jk+jn,kk+kn);

                Ltarget[Lk] += a_aux[npk] * Msource[Mn];
              }
            }
          } // end inner loop
        }
      }
    } // end outer loop
  }

  /** Kernel L2L operation
   * L_t += Op(L_s) where L_t is the target and L_s is the source
   *
   * @param[in] source The local source at the parent level
   * @param[in,out] target The local target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Lsource includes the influence of all points outside its box
   */
  void L2L(const local_type& Lsource,
                 local_type& Ltarget,
           const point_type& translation) const {
    for (unsigned i=0; i<MTERMS; i++) {
      // n = (I[i], J[i], K[i])
      unsigned n[3] = {I[i],J[i],K[i]};
      for (int ii=I[i]; ii<P+1; ii++) {
        for (int jj=J[i]; jj<P+1-ii; jj++) {
          for (int kk=K[i]; kk<P+1-ii-jj; kk++) {
            // k = (ii, jj, kk)
            int Lk = index_cache(ii,jj,kk);

            double fact_term = fact[ii-n[0]]*fact[jj-n[1]]*fact[kk-n[2]];
            Ltarget[i] += Lsource[Lk]*pow(translation[0],ii-n[0])*pow(translation[1],jj-n[1])*pow(translation[2],kk-n[2]) / fact_term;
          }
        }
      }
    }
  }

  /** Kernel L2T operation
   * r += Op(L, t) where L is the local expansion and r is the result
   *
   * @param[in] L The local expansion
   * @param[in] center The center of the box with the local expansion
   * @param[in] target The target of this L2T operation
   * @param[in] result The result to accumulate into
   * @pre L includes the influence of all sources outside its box
   */
  void L2T(const local_type& L, const point_type& center,
           const target_type& target, result_type& result) const {
    auto dx = target-center;

    for (unsigned i=0; i<MTERMS; i++) {
      unsigned k[3] = {I[i], J[i], K[i]};
      // calculate potential
      double phi = L[i]*pow(dx[0],k[0])*pow(dx[1],k[1])*pow(dx[2],k[2]) / (fact[k[0]]*fact[k[1]]*fact[k[2]]);
      result[0] += phi;

      double inv[3];
      for (int i=0; i<3; i++) inv[i] = (fabs(dx[i]) < 1e-12) ? 0 : 1. / dx[i];

      // using potential, calculate derivatives
      result[1] += phi * k[0] * inv[0];
      result[2] += phi * k[1] * inv[1];
      result[3] += phi * k[2] * inv[2];
    }
  }

 protected:
  // get coefficients given by d^n / dx^n f(x)
  void getCoeff(real_vec& a, real_vec& ax, real_vec& ay, real_vec& az,
                const point_type& dX) const
  {
    real_vec b(a.size(),0);
    real dx = dX[0], dy = dX[1], dz = dX[2];

    auto R2 = norm_2_sq(dX);
    auto R  = sqrt(R2);

    int i,j,k,I,Im1x,Im2x,Im1y,Im2y,Im1z,Im2z;
    real C,C1,C2,Cb, R2_1;

    R2_1 = 1/R2;

    // First coefficient
    b[0] = exp(-kappa*R);
    a[0] = b[0]/R;

    // Two indices = 0
    I = index_cache(1,0,0); // setIndex(P,1,0,0);
    b[I]   = -kappa * (dx*a[0]); // 1,0,0
    b[P+1] = -kappa * (dy*a[0]); // 0,1,0
    b[1]   = -kappa * (dz*a[0]); // 0,0,1

    a[I]   = -R2_1*dx*(kappa*b[0]+a[0]);
    a[P+1] = -R2_1*dy*(kappa*b[0]+a[0]);
    a[1]   = -R2_1*dz*(kappa*b[0]+a[0]);

    ax[0]  = a[I];
    ay[0]  = a[P+1];
    az[0]  = a[1];

    for (i=2; i<P+1; i++)
    {
      Cb   = -kappa/i;
      C    = R2_1/i;
      I    = index_cache(i,0,0); // setIndex(P,i,0,0);
      Im1x = index_cache(i-1,0,0); // setIndex(P,i-1,0,0);
      Im2x = index_cache(i-2,0,0); // setIndex(P,i-2,0,0);
      b[I] = Cb * (dx*a[Im1x] + a[Im2x]);
      a[I] = C * ( -kappa*(dx*b[Im1x] + b[Im2x]) -(2*i-1)*dx*a[Im1x] - (i-1)*a[Im2x] );
      ax[Im1x] = a[I]*i;

      I    = index_cache(0,i,0); // setIndex(P,0,i,0);
      Im1y = I-(P+2-i);
      Im2y = Im1y-(P+2-i+1);
      b[I] = Cb * (dy*a[Im1y] + a[Im2y]);
      a[I] = C * ( -kappa*(dy*b[Im1y] + b[Im2y]) -(2*i-1)*dy*a[Im1y] - (i-1)*a[Im2y] );
      ay[Im1y] = a[I]*i;

      I   = i;
      Im1z = I-1;
      Im2z = I-2;
      b[I] = Cb * (dz*a[Im1z] + a[Im2z]);
      a[I] = C * ( -kappa*(dz*b[Im1z] + b[Im2z]) -(2*i-1)*dz*a[Im1z] - (i-1)*a[Im2z] );
      az[Im1z] = a[I]*i;
    }

    // One index = 0, one = 1 other >=1
    Cb   = -kappa/2;
    C    = R2_1/2.;
    I    = index_cache(1,1,0); //setIndex(P,1,1,0);
    Im1x = P+1;
    Im1y = I-P;
    b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y]);
    a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]) - 3*(dx*a[Im1x]+dy*a[Im1y]) );
    ax[Im1x] = a[I];
    ay[Im1y] = a[I];
    I    = index_cache(1,0,1); // setIndex(P,1,0,1);
    Im1x = 1;
    Im1z = I-1;
    b[I] = Cb * (dx*a[Im1x] + dz*a[Im1z]);
    a[I] = C * ( -kappa*(dx*b[Im1x]+dz*b[Im1z]) - 3*(dx*a[Im1x]+dz*a[Im1z]) );
    ax[Im1x] = a[I];
    az[Im1z] = a[I];
    I    = index_cache(0,1,1); // setIndex(P,0,1,1);
    Im1y = I-(P+1);
    Im1z = I-1;
    b[I] = Cb * (dy*a[Im1y] + dz*a[Im1z]);
    a[I] = C * ( -kappa*(dy*b[Im1y]+dz*b[Im1z]) - 3*(dy*a[Im1y]+dz*a[Im1z]) );
    ay[Im1y] = a[I];
    az[Im1z] = a[I];

    for (i=2; i<P; i++)
    {
      Cb   = -kappa/(i+1);
      C    = R2_1/(1+i);
      C1   = 1+2*i;
      I    = index_cache(1,i,0); // setIndex(P,1,i,0);
      Im1x = index_cache(0,i,0); // setIndex(P,0,i,0);
      Im1y = I-(P+1-i);
      Im2y = Im1y-(P+2-i);
      b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + a[Im2y]);
      a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+b[Im2y]) - C1*(dx*a[Im1x]+dy*a[Im1y]) - (1+i-1)*(a[Im2y]) );
      ax[Im1x] = a[I];
      ay[Im1y] = a[I]*i;

      I    = index_cache(1,0,i); // setIndex(P,1,0,i);
      Im1x = index_cache(0,0,i); // setIndex(P,0,0,i);
      Im1z = I-1;
      Im2z = I-2;
      b[I] = Cb * (dx*a[Im1x] + dz*a[Im1z] + a[Im2z]);
      a[I] = C * ( -kappa*(dx*b[Im1x]+dz*b[Im1z]+b[Im2z]) - C1*(dx*a[Im1x]+dz*a[Im1z]) - (1+i-1)*(a[Im2z]) );
      ax[Im1x] = a[I];
      az[Im1z] = a[I]*i;

      I    = index_cache(0,1,i); // setIndex(P,0,1,i);
      Im1y = I-(P+1);
      Im1z = I-1;
      Im2z = I-2;
      b[I] = Cb * (dy*a[Im1y] + dz*a[Im1z] + a[Im2z]);
      a[I] = C * ( -kappa*(dy*b[Im1y]+dz*b[Im1z]+b[Im2z]) - C1*(dy*a[Im1y]+dz*a[Im1z]) - (1+i-1)*(a[Im2z]) );
      ay[Im1y] = a[I];
      az[Im1z] = a[I]*i;

      I    = index_cache(i,1,0); // setIndex(P,i,1,0);
      Im1y = I-(P+1-i);
      Im1x = index_cache(i-1,1,0); // setIndex(P,i-1,1,0);
      Im2x = index_cache(i-2,1,0); // setIndex(P,i-2,1,0);
      b[I] = Cb * (dy*a[Im1y] + dx*a[Im1x] + a[Im2x]);
      a[I] = C * ( -kappa*(dy*b[Im1y]+dx*b[Im1x]+b[Im2x]) - C1*(dy*a[Im1y]+dx*a[Im1x]) - (1+i-1)*(a[Im2x]) );
      ax[Im1x] = a[I]*i;
      ay[Im1y] = a[I];

      I    = index_cache(i,0,1); // setIndex(P,i,0,1);
      Im1z = I-1;
      Im1x = index_cache(i-1,0,1); // setIndex(P,i-1,0,1);
      Im2x = index_cache(i-2,0,1); // setIndex(P,i-2,0,1);
      b[I] = Cb * (dz*a[Im1z] + dx*a[Im1x] + a[Im2x]);
      a[I] = C * ( -kappa*(dz*b[Im1z]+dx*b[Im1x]+b[Im2x]) - C1*(dz*a[Im1z]+dx*a[Im1x]) - (1+i-1)*(a[Im2x]) );
      ax[Im1x] = a[I]*i;
      az[Im1z] = a[I];

      I    = index_cache(0,i,1); // setIndex(P,0,i,1);
      Im1z = I-1;
      Im1y = I-(P+2-i);
      Im2y = Im1y-(P+3-i);
      b[I] = Cb * (dz*a[Im1z] + dy*a[Im1y] + a[Im2y]);
      a[I] = C * ( -kappa*(dz*b[Im1z]+dy*b[Im1y]+b[Im2y]) - C1*(dz*a[Im1z]+dy*a[Im1y]) - (1+i-1)*(a[Im2y]) );
      ay[Im1y] = a[I]*i;
      az[Im1z] = a[I];
    }

    // One index 0, others >=2
    for (i=2; i<P+1; i++)
    {
      for (j=2; j<P+1-i; j++)
      {
        Cb   = -kappa/(i+j);
        C    = R2_1/(i+j);
        C1   = 2*(i+j)-1;
        I    = index_cache(i,j,0); // setIndex(P,i,j,0);
        Im1x = index_cache(i-1,j,0); // setIndex(P,i-1,j,0);
        Im2x = index_cache(i-2,j,0); // setIndex(P,i-2,j,0);
        Im1y = I-(P+2-j-i);
        Im2y = Im1y-(P+3-j-i);
        b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + a[Im2x] + a[Im2y]);
        a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+b[Im2x]+b[Im2y]) - C1*(dx*a[Im1x]+dy*a[Im1y]) -(i+j-1)*(a[Im2x]+a[Im2y]) );
        ax[Im1x] = a[I]*i;
        ay[Im1y] = a[I]*j;

        I    = index_cache(i,0,j); // setIndex(P,i,0,j);
        Im1x = index_cache(i-1,0,j); // setIndex(P,i-1,0,j);
        Im2x = index_cache(i-2,0,j); // setIndex(P,i-2,0,j);
        Im1z = I-1;
        Im2z = I-2;
        b[I] = Cb * (dx*a[Im1x] + dz*a[Im1z] + a[Im2x] + a[Im2z]);
        a[I] = C * ( -kappa*(dx*b[Im1x]+dz*b[Im1z]+b[Im2x]+b[Im2z]) - C1*(dx*a[Im1x]+dz*a[Im1z]) -(i+j-1)*(a[Im2x]+a[Im2z]) );
        ax[Im1x] = a[I]*i;
        az[Im1z] = a[I]*j;

        I    = index_cache(0,i,j); // setIndex(P,0,i,j);
        Im1y = I-(P+2-i);
        Im2y = Im1y-(P+3-i);
        Im1z = I-1;
        Im2z = I-2;
        b[I] = Cb * (dy*a[Im1y] + dz*a[Im1z] + a[Im2y] + a[Im2z]);
        a[I] = C * ( -kappa*(dy*b[Im1y]+dz*b[Im1z]+b[Im2y]+b[Im2z]) - C1*(dy*a[Im1y]+dz*a[Im1z]) -(i+j-1)*(a[Im2y]+a[Im2z]) );
        ay[Im1y] = a[I]*i;
        az[Im1z] = a[I]*j;
      }
    }

    if (P>2)
    {
      // Two index = 1, other>=1
      C    = R2_1/3;
      Cb   = -kappa/3;
      I    = index_cache(1,1,1); // setIndex(P,1,1,1);
      Im1x = index_cache(0,1,1); // setIndex(P,0,1,1);
      Im1y = index_cache(1,0,1); // setIndex(P,1,0,1);
      Im1y = I-(P);
      Im1z = I-1;
      b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z]);
      a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]) - 5*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) );
      ax[Im1x] = a[I];
      ay[Im1y] = a[I];
      az[Im1z] = a[I];
      for (i=2; i<P-1; i++)
      {
        Cb   = -kappa/(2+i);
        C    = R2_1/(i+2);
        C1   = 2*i+3;
        I    = index_cache(i,1,1); // setIndex(P,i,1,1);
        Im1x = index_cache(i-1,1,1); // setIndex(P,i-1,1,1);
        Im1y = I-(P+1-i);
        Im1z = I-1;
        Im2x = index_cache(i-2,1,1); // setIndex(P,i-2,1,1);
        b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x]);
        a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]) - C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - (i+1)*(a[Im2x]) );
        ax[Im1x] = a[I]*i;
        ay[Im1y] = a[I];
        az[Im1z] = a[I];

        I    = index_cache(1,i,1); // setIndex(P,1,i,1);
        Im1x = index_cache(0,i,1); // setIndex(P,0,i,1);
        Im1y = I-(P+1-i);
        Im2y = Im1y-(P+2-i);
        Im1z = I-1 ;
        b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2y]);
        a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2y]) - C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - (i+1)*(a[Im2y]) );
        ax[Im1x] = a[I];
        ay[Im1y] = a[I]*i;
        az[Im1z] = a[I];

        I    = index_cache(1,1,i); // setIndex(P,1,1,i);
        Im1x = index_cache(0,1,i); // setIndex(P,0,1,i);
        Im1y = I-(P);
        Im1z = I-1;
        Im2z = I-2;
        b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2z]);
        a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2z]) - C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - (i+1)*(a[Im2z]) );
        ax[Im1x] = a[I];
        ay[Im1y] = a[I];
        az[Im1z] = a[I]*i;
      }
    }

    // One index = 1, others >=2
    if (P>4)
    {
      for (i=2; i<P-2; i++)
      {
        for (j=2; j<P-i; j++)
        {
          Cb = -kappa/(1+i+j);
          C  =  R2_1/(1+i+j);
          C1 = -(2.*(i+j)+1);
          C2 = (i+j);
          I    = index_cache(1,i,j); // setIndex(P,1,i,j);
          Im1x = index_cache(0,i,j); // setIndex(P,0,i,j);
          Im1y = I-(P+1-i);
          Im2y = Im1y-(P+2-i);
          Im1z = I-1;
          Im2z = I-2;
          b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2y] + a[Im2z]);
          a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2y]+b[Im2z]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2y]+a[Im2z]) );
          ax[Im1x] = a[I];
          ay[Im1y] = a[I]*i;
          az[Im1z] = a[I]*j;

          I    = index_cache(i,1,j); // setIndex(P,i,1,j);
          Im1x = index_cache(i-1,1,j); // setIndex(P,i-1,1,j);
          Im1y = I-(P+1-i);
          Im2x = index_cache(i-2,1,j); // setIndex(P,i-2,1,j);
          Im1z = I-1;
          Im2z = I-2;
          b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x] + a[Im2z]);
          a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]+b[Im2z]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2x]+a[Im2z]) );
          ax[Im1x] = a[I]*i;
          ay[Im1y] = a[I];
          az[Im1z] = a[I]*j;

          I    = index_cache(i,j,1); // setIndex(P,i,j,1);
          Im1x = index_cache(i-1,j,1); // setIndex(P,i-1,j,1);
          Im2x = index_cache(i-2,j,1); // setIndex(P,i-2,j,1);
          Im1y = I-(P+2-i-j);
          Im2y = Im1y-(P+3-i-j);
          Im1z = I-1;
          b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x] + a[Im2y]);
          a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]+b[Im2y]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2x]+a[Im2y]) );
          ax[Im1x] = a[I]*i;
          ay[Im1y] = a[I]*j;
          az[Im1z] = a[I];
        }
      }
    }

    // All indices >= 2
    if (P>5)
    {
      for (i=2;i<P-3;i++)
      {
        for (j=2;j<P-1-i;j++)
        {
          for (k=2;k<P+1-i-j;k++)
          {
            Cb = -kappa/(i+j+k);
            C  = R2_1/(i+j+k);
            C1 = -(2.*(i+j+k)-1);
            C2 = i+j+k-1.;
            I    = index_cache(i,j,k); // setIndex(P,i,j,k);
            Im1x = index_cache(i-1,j,k); // setIndex(P,i-1,j,k);
            Im2x = index_cache(i-2,j,k); // setIndex(P,i-2,j,k);
            Im1y = I-(P+2-i-j);
            Im2y = Im1y-(P+3-i-j);
            Im1z = I-1;
            Im2z = I-2;
            b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x] + a[Im2y] + a[Im2z]);
            a[I] = C * ( -kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]+b[Im2y]+b[Im2z]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2x]+a[Im2y]+a[Im2z]) );
            ax[Im1x] = a[I]*i;
            ay[Im1y] = a[I]*j;
            az[Im1z] = a[I]*k;
          }
        }
      }
    }
  }
};
