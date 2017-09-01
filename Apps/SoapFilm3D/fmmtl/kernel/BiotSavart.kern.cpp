#ifndef FMMTL_BIOTSAVART_KERN
#define FMMTL_BIOTSAVART_KERN
//#define FMMTL_KERNEL

#include "fmmtl/Kernel.hpp"

#include "fmmtl/numeric/Vec.hpp"

struct BiotSavart
    : public fmmtl::Kernel<BiotSavart> {
  typedef Vec<3,double> source_type;
  typedef Vec<3,double> target_type;
  typedef Vec<3,double> charge_type;
  typedef Vec<3,double> result_type;

  struct kernel_value_type {
    Vec<3,double> v;
    FMMTL_INLINE
    kernel_value_type(const Vec<3,double>& _v) : v(_v) {}

    FMMTL_INLINE
    result_type operator*(const charge_type& c) const {
      return cross(v, c);
    }
  };

  FMMTL_INLINE
  kernel_value_type operator()(const target_type& t,
                               const source_type& s) const {
    Vec<3,double> dist = s - t;            //   Vector from target to source
    double R2 = norm_2_sq(dist);              //   R^2
    double invR2 = 1.0 / R2;               //   1 / R^2
    if (R2 < 1e-20) invR2 = 0;             //   Exclude self interaction
    dist *= invR2 * std::sqrt(invR2);
    return kernel_value_type(dist);
  }

  FMMTL_INLINE
  kernel_value_type transpose(const kernel_value_type& kts) const {
    return kernel_value_type(-kts.v);
  }
};
FMMTL_KERNEL_EXTRAS(BiotSavart);


struct RosenheadMoore
    : public fmmtl::Kernel<RosenheadMoore> {
  typedef Vec<3,double> source_type;
  typedef Vec<3,double> target_type;
  typedef Vec<3,double> charge_type;
  typedef Vec<3,double> result_type;

  struct kernel_value_type {
    Vec<3,double> v;
    FMMTL_INLINE
    kernel_value_type(const Vec<3,double>& _v) : v(_v) {}

    FMMTL_INLINE
    result_type operator*(const charge_type& c) const {
      return cross(v, c);
    }
  };

  double aSq;

  FMMTL_INLINE
  RosenheadMoore(double _a = 1) : aSq(_a * _a) {}

  FMMTL_INLINE
  kernel_value_type operator()(const target_type& t,
                               const source_type& s) const {
    Vec<3,double> dist = s - t;            //   Vector from target to source
    double R2 = aSq + norm_2_sq(dist);
    dist /= R2 * std::sqrt(R2);
    return kernel_value_type(dist);
  }

  FMMTL_INLINE
  kernel_value_type transpose(const kernel_value_type& kts) const {
    return kernel_value_type(-kts.v);
  }
};
FMMTL_KERNEL_EXTRAS(RosenheadMoore);


#endif
