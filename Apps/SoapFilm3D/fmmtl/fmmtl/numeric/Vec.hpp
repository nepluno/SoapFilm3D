#pragma once
/** @file Vec.hpp
 * @brief A small N-dimensional numerical vector type that works on CPU and GPU.
 *
 * TODO: Place in fmmtl namespace.
 */

#include <iostream>
#include <cmath>

#include "fmmtl/config.hpp"

#define for_i    for(std::size_t i = 0; i != N; ++i)
#define for_I(K) for(std::size_t i = K; i != N; ++i)

/** @class Vec
 * @brief Class representing ND points and vectors.
 */
template <std::size_t N, typename T>
struct Vec {
  FMMTL_STATIC_ASSERT(N > 0, "Vec<N,T> needs N >= 1");

  T elem[N];

  typedef T               value_type;
  typedef T&              reference;
  typedef const T&        const_reference;
  typedef T*              iterator;
  typedef const T*        const_iterator;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;

  // CONSTRUCTORS

  // TODO: require zero-initialization?
  FMMTL_INLINE Vec() {
    for_i elem[i] = T(0);
    //std::fill(this->begin(), this->end(), value_type());
  }
  // TODO: Force 0-initialization to get POD and trivial semantics?
  //FMMTL_INLINE Vec() = default;
  FMMTL_INLINE explicit Vec(const value_type& b) {
    for_i elem[i] = b;
    //std::fill(this->begin(), this->end(), b);
  }
  FMMTL_INLINE Vec(const value_type& b0,
                   const value_type& b1) {
    FMMTL_STATIC_ASSERT(N >= 2, "Too many arguments to Vec constructor");
    elem[0] = b0; elem[1] = b1;
    for_I(2) elem[i] = T(0);
    //std::fill(this->begin() + 2, this->end(), value_type());
  }
  FMMTL_INLINE Vec(const value_type& b0,
                   const value_type& b1,
                   const value_type& b2) {
    FMMTL_STATIC_ASSERT(N >= 3,  "Too many arguments to Vec constructor");
    elem[0] = b0; elem[1] = b1; elem[2] = b2;
    for_I(3) elem[i] = T(0);
    //std::fill(this->begin() + 3, this->end(), value_type());
  }
  FMMTL_INLINE Vec(const value_type& b0,
                   const value_type& b1,
                   const value_type& b2,
                   const value_type& b3) {
    FMMTL_STATIC_ASSERT(N >= 4,  "Too many arguments to Vec constructor");
    elem[0] = b0; elem[1] = b1; elem[2] = b2; elem[3] = b3;
    for_I(4) elem[i] = T(0);
    //std::fill(this->begin() + 4, this->end(), value_type());
  }
  template <typename U>
  FMMTL_INLINE explicit Vec(const Vec<N,U>& v) {
    for_i elem[i] = v[i];
    //std::copy(v.begin(), v.end(), this->begin());
  }
  /*
  template <typename Generator>
  FMMTL_INLINE explicit Vec(const Generator& gen) {
    for_i elem[i] = gen(i);
  }
  */
  template <typename U, typename OP>
  FMMTL_INLINE explicit Vec(const Vec<N,U>& v, OP op) {
    for_i elem[i] = op(v[i]);
  }
  template <typename U1, typename OP, typename U2>
  FMMTL_INLINE explicit Vec(const Vec<N,U1>& v1, const Vec<N,U2>& v2, OP op) {
    for_i elem[i] = op(v1[i], v2[i]);
  }
  template <typename U1, typename OP, typename U2>
  FMMTL_INLINE explicit Vec(const U1& v1, const Vec<N,U2>& v2, OP op) {
    for_i elem[i] = op(v1, v2[i]);
  }
  template <typename U1, typename OP, typename U2>
  FMMTL_INLINE explicit Vec(const Vec<N,U1>& v1, const U2& v2, OP op) {
    for_i elem[i] = op(v1[i], v2);
  }

  // COMPARATORS

  FMMTL_INLINE bool operator==(const Vec& b) const {
    for_i { if (!(elem[i] == b[i])) return false; }
    return true;
    //return std::equal(this->begin(), this->end(), b.begin());
  }
  FMMTL_INLINE bool operator!=(const Vec& b) const {
    return !(*this == b);
  }

  // MODIFIERS

  /** Add scalar @a b to this Vec */
  template <typename U>
  FMMTL_INLINE Vec& operator+=(const U& b) {
    for_i elem[i] += b;
    return *this;
  }
  /** Subtract scalar @a b from this Vec */
  template <typename U>
  FMMTL_INLINE Vec& operator-=(const U& b) {
    for_i elem[i] -= b;
    return *this;
  }
  /** Scale this Vec up by scalar @a b */
  template <typename U>
  FMMTL_INLINE Vec& operator*=(const U& b) {
    for_i elem[i] *= b;
    return *this;
  }
  /** Scale this Vec down by scalar @a b */
  template <typename U>
  FMMTL_INLINE Vec& operator/=(const U& b) {
    for_i elem[i] /= b;
    return *this;
  }
  /** Add Vec @a b to this Vec */
  template <typename U>
  FMMTL_INLINE Vec& operator+=(const Vec<N,U>& b) {
    for_i elem[i] += b[i];
    return *this;
  }
  /** Subtract Vec @a b from this Vec */
  template <typename U>
  FMMTL_INLINE Vec& operator-=(const Vec<N,U>& b) {
    for_i elem[i] -= b[i];
    return *this;
  }
  /** Scale this Vec up by factors in @a b */
  template <typename U>
  FMMTL_INLINE Vec& operator*=(const Vec<N,U>& b) {
    for_i elem[i] *= b[i];
    return *this;
  }
  /** Scale this Vec down by factors in @a b */
  template <typename U>
  FMMTL_INLINE Vec& operator/=(const Vec<N,U>& b) {
    for_i elem[i] /= b[i];
    return *this;
  }

  // ACCESSORS

  FMMTL_INLINE reference       operator[](size_type i)       { return elem[i]; }
  FMMTL_INLINE const_reference operator[](size_type i) const { return elem[i]; }

  FMMTL_INLINE T*       data()       { return elem; }
  FMMTL_INLINE const T* data() const { return elem; }

  FMMTL_INLINE reference       front()       { return elem[0]; }
  FMMTL_INLINE const_reference front() const { return elem[0]; }
  FMMTL_INLINE reference        back()       { return elem[N-1]; }
  FMMTL_INLINE const_reference  back() const { return elem[N-1]; }

  FMMTL_INLINE static size_type     size() { return N; }
  FMMTL_INLINE static size_type max_size() { return N; }
  FMMTL_INLINE static bool         empty() { return false; }

  // ITERATORS

  FMMTL_INLINE iterator        begin()       { return elem; }
  FMMTL_INLINE const_iterator  begin() const { return elem; }
  FMMTL_INLINE const_iterator cbegin() const { return elem; }

  FMMTL_INLINE iterator          end()       { return elem+N; }
  FMMTL_INLINE const_iterator    end() const { return elem+N; }
  FMMTL_INLINE const_iterator   cend() const { return elem+N; }
};

// OPERATORS

/** Write a Vec to an output stream */
template <std::size_t N, typename T>
inline std::ostream& operator<<(std::ostream& s, const Vec<N,T>& a) {
  s << a[0];
  for (unsigned i = 1; i != a.size(); ++i) s << " " << a[i];
  return s;
}
/** Read a Vec from an input stream */
template <std::size_t N, typename T>
inline std::istream& operator>>(std::istream& s, Vec<N,T>& a) {
  for_i s >> a[i];
  return s;
}

/** Compute cross product of two 3D Vecs */
template <typename T>
FMMTL_INLINE Vec<3,T> cross(const Vec<3,T>& a, const Vec<3,T>& b) {
  return Vec<3,T>(a[1]*b[2] - a[2]*b[1],
                  a[2]*b[0] - a[0]*b[2],
                  a[0]*b[1] - a[1]*b[0]);
}

// ARITHMETIC UNARY OPERATORS

/** Unary negation: Return -@a a */
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,T> operator-(Vec<N,T> a) {
  for_i a[i] = -a[i];
  return a;
}
/** Unary plus: Return @a a. ("+a" should work if "-a" works.) */
template <std::size_t N, typename T>
FMMTL_INLINE const Vec<N,T>& operator+(const Vec<N,T>& a) {
  return a;
}

// ARITHEMTIC BINARY OPERATORS

#if !defined(__CUDACC__)
/** This is not being compiled with CUDA, but with C++11 compatible compiler.
 * Type promotion with SFINAE selection can be accomplished with
 *   template <typename T, typename U>
 *   using prod_type = decltype(std::declval<T>() * std::declval<U>());
 * but we use a C++03 style for compatibility with the nvcc version below.
 */
#  define FMMTL_BINARY_PROMOTE_DECLARE(NAME, OP)                            \
  template <typename T, typename U,                                         \
            typename R = decltype(std::declval<T>() OP std::declval<U>())>  \
  struct NAME##_op {                                                        \
    typedef R type;                                                         \
    FMMTL_INLINE R operator()(const T& a, const U& b) const {               \
      return a OP b;                                                        \
    }                                                                       \
  }

#  define FMMTL_UNARY_PROMOTE_DECLARE(NAME, OP)                             \
  template <typename T,                                                     \
            typename R = decltype(OP(std::declval<T>()))>                   \
  struct NAME##_op {                                                        \
    typedef R type;                                                         \
    FMMTL_INLINE R operator()(const T& a) const {                           \
      return OP(a);                                                         \
    }                                                                       \
  }
#else
/** This is being compiled by CUDA, which does not have decltype.
 * Instead, a simple fix is to disallow type promotion and cross our fingers.
 * TODO: Improve.
 */
#  define FMMTL_BINARY_PROMOTE_DECLARE(NAME, OP)                            \
  template <typename T, typename U>                                         \
  struct NAME##_op {};                                                      \
  template <typename T>                                                     \
  struct NAME##_op<T,T> {                                                   \
    typedef T type;                                                         \
    FMMTL_INLINE T operator()(const T& a, const T& b) const {               \
      return a OP b;                                                        \
    }                                                                       \
  }

#  define FMMTL_UNARY_PROMOTE_DECLARE(NAME, OP)                             \
  template <typename T>                                                     \
  struct NAME##_op {                                                        \
    typedef T type;                                                         \
    FMMTL_INLINE T operator()(const T& a) const {                           \
      return OP(a);                                                         \
    }                                                                       \
  }
#endif

namespace fmmtl {

FMMTL_BINARY_PROMOTE_DECLARE(sum,+);
FMMTL_BINARY_PROMOTE_DECLARE(diff,-);
FMMTL_BINARY_PROMOTE_DECLARE(prod,*);
FMMTL_BINARY_PROMOTE_DECLARE(div,/);

}

// TODO: namespace Vec and all operators

template <std::size_t N, typename T, typename U>
FMMTL_INLINE Vec<N,typename fmmtl::sum_op<T,U>::type>
operator+(const Vec<N,T>& a, const Vec<N,U>& b) {
  return Vec<N,typename fmmtl::sum_op<T,U>::type>(a, b, fmmtl::sum_op<T,U>());
}
template <std::size_t N, typename T, typename U>
FMMTL_INLINE Vec<N,typename fmmtl::sum_op<T,U>::type>
operator+(const Vec<N,T>& a, const U& b) {
  return Vec<N,typename fmmtl::sum_op<T,U>::type>(a, b, fmmtl::sum_op<T,U>());
}
template <std::size_t N, typename T, typename U>
FMMTL_INLINE Vec<N,typename fmmtl::sum_op<T,U>::type>
operator+(const T& a, const Vec<N,U>& b) {
  return Vec<N,typename fmmtl::sum_op<T,U>::type>(a, b, fmmtl::sum_op<T,U>());
}
template <std::size_t N, typename T, typename U>
FMMTL_INLINE Vec<N,typename fmmtl::diff_op<T,U>::type>
operator-(const Vec<N,T>& a, const Vec<N,U>& b) {
  return Vec<N,typename fmmtl::diff_op<T,U>::type>(a, b, fmmtl::diff_op<T,U>());
}
template <std::size_t N, typename T, typename U>
FMMTL_INLINE Vec<N,typename fmmtl::diff_op<T,U>::type>
operator-(const Vec<N,T>& a, const U& b) {
  return Vec<N,typename fmmtl::diff_op<T,U>::type>(a, b, fmmtl::diff_op<T,U>());
}
template <std::size_t N, typename T, typename U>
FMMTL_INLINE Vec<N,typename fmmtl::diff_op<T,U>::type>
operator-(const T& a, const Vec<N,U>& b) {
  return Vec<N,typename fmmtl::diff_op<T,U>::type>(a, b, fmmtl::diff_op<T,U>());
}
template <std::size_t N, typename T, typename U>
FMMTL_INLINE Vec<N,typename fmmtl::prod_op<T,U>::type>
operator*(const Vec<N,T>& a, const Vec<N,U>& b) {
  return Vec<N,typename fmmtl::prod_op<T,U>::type>(a, b, fmmtl::prod_op<T,U>());
}
template <std::size_t N, typename T, typename U>
FMMTL_INLINE Vec<N,typename fmmtl::prod_op<T,U>::type>
operator*(const Vec<N,T>& a, const U& b) {
  return Vec<N,typename fmmtl::prod_op<T,U>::type>(a, b, fmmtl::prod_op<T,U>());
}
template <std::size_t N, typename T, typename U>
FMMTL_INLINE Vec<N,typename fmmtl::prod_op<T,U>::type>
operator*(const T& a, const Vec<N,U>& b) {
  return Vec<N,typename fmmtl::prod_op<T,U>::type>(a, b, fmmtl::prod_op<T,U>());
}
template <std::size_t N, typename T, typename U>
FMMTL_INLINE Vec<N,typename fmmtl::div_op<T,U>::type>
operator/(const Vec<N,T>& a, const Vec<N,U>& b) {
  return Vec<N,typename fmmtl::div_op<T,U>::type>(a, b, fmmtl::div_op<T,U>());
}
template <std::size_t N, typename T, typename U>
FMMTL_INLINE Vec<N,typename fmmtl::div_op<T,U>::type>
operator/(const Vec<N,T>& a, const U& b) {
  return Vec<N,typename fmmtl::div_op<T,U>::type>(a, b, fmmtl::div_op<T,U>());
}
template <std::size_t N, typename T, typename U>
FMMTL_INLINE Vec<N,typename fmmtl::div_op<T,U>::type>
operator/(const T& a, const Vec<N,U>& b) {
  return Vec<N,typename fmmtl::div_op<T,U>::type>(a, b, fmmtl::div_op<T,U>());
}

// ELEMENTWISE OPERATORS

//namespace fmmtl {

FMMTL_UNARY_PROMOTE_DECLARE(abs, abs);
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,typename abs_op<T>::type>
abs(const Vec<N,T>& a) {
  return Vec<N,typename abs_op<T>::type>(a, abs_op<T>());
}

FMMTL_UNARY_PROMOTE_DECLARE(sqrt, sqrt);
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,typename sqrt_op<T>::type>
sqrt(const Vec<N,T>& a) {
  return Vec<N,typename sqrt_op<T>::type>(a, sqrt_op<T>());
}

FMMTL_UNARY_PROMOTE_DECLARE(conj, conj);
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,typename conj_op<T>::type>
conj(const Vec<N,T>& a) {
  return Vec<N,typename conj_op<T>::type>(a, conj_op<T>());
}

FMMTL_UNARY_PROMOTE_DECLARE(real, real);
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,typename real_op<T>::type>
real(const Vec<N,T>& a) {
  return Vec<N,typename real_op<T>::type>(a, real_op<T>());
}

FMMTL_UNARY_PROMOTE_DECLARE(imag, imag);
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,typename imag_op<T>::type>
imag(const Vec<N,T>& a) {
  return Vec<N,typename imag_op<T>::type>(a, imag_op<T>());
}


//} // end namespace fmmtl


// Compliance with std::
namespace std {

template <std::size_t I, std::size_t N, typename T>
typename Vec<N,T>::reference
get(Vec<N,T>& a) {
  FMMTL_STATIC_ASSERT(I < N, "I must be less than N.");
  return a[I];
}

template <std::size_t I, std::size_t N, typename T>
typename Vec<N,T>::const_reference
get(const Vec<N,T>& a) {
  FMMTL_STATIC_ASSERT(I < N, "I must be less than N.");
  return a[I];
}

template <std::size_t N, typename T>
void
swap(Vec<N,T>& a, Vec<N,T>& b) {
  std::swap_ranges(a.begin(), a.end(), b.begin());
}

} // end namespace std


// META OPERATIONS

#include "fmmtl/meta/dimension.hpp"

namespace fmmtl {

template <std::size_t N, typename T>
struct dimension<Vec<N,T> > {
  const static std::size_t value = N;
};

} // end namespace fmmtl


#undef for_i
#undef FMMTL_BINARY_PROMOTE_DECLARE
#undef FMMTL_UNARY_PROMOTE_DECLARE

#include "fmmtl/numeric/norm.hpp"
#include "fmmtl/numeric/random.hpp"
