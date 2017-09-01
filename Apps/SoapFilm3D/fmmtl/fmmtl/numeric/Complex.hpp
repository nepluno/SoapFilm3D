#pragma once
/*
 *  Copyright 2008-2013 NVIDIA Corporation
 *  Copyright 2013 Filipe RNC Maia
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/*! \file complex.h
 *  \brief Complex numbers
 */

#include <cmath>
#include <complex>
#include <sstream>

#include "fmmtl/config.hpp"

namespace fmmtl {

/*! \p complex is equivalent to <tt>std::complex</tt>. It is functionally
 *  equivalent to it, but can also be used in device code which <tt>std::complex</tt> currently cannot.
 *
 *  \tparam T The type used to hold the real and imaginary parts. Should be <tt>float</tt>
 *  or <tt>double</tt>. Others types are not supported.
 *
 */
template <typename T>
struct complex {
 public:

  /*! \p value_type is the type of \p complex's real and imaginary parts.
   */
  typedef T value_type;

  /* --- Constructors --- */

  /*! Construct a complex number from its real and imaginary parts.
   *
   *  \param re The real part of the number.
   *  \param im The imaginary part of the number.
   */
  FMMTL_INLINE
  complex(const T& re = T(), const T& im = T()) {
    real(re); imag(im);
  }

  /*! This copy constructor copies from a \p complex with a type that
   *  is convertible to this \p complex \c value_type.
   *
   *  \param z The \p complex to copy from.
   *
   *  \tparam X is convertible to \c value_type.
   */
  template <typename X>
  FMMTL_INLINE
  complex(const complex<X>& z) {
    real(static_cast<T>(z.real())); imag(static_cast<T>(z.imag()));
  }

  /*! This copy constructor copies from a <tt>std::complex</tt> with a type that
   *  is convertible to this \p complex \c value_type.
   *
   *  \param z The \p complex to copy from.
   *
   *  \tparam X is convertible to \c value_type.
   */
  template <typename X>
  FMMTL_INLINE
  complex(const std::complex<X>& z) {
    real(static_cast<T>(z.real())); imag(static_cast<T>(z.imag()));
  }

  /* --- Compound Assignment Operators --- */

  /*! Adds a \p complex to this \p complex and
   *  assigns the result to this \p complex.
   *
   *  \param z The \p complex to be Added.
   */
  FMMTL_INLINE
  complex<T>& operator+=(const complex<T>& z) {
    real(real()+z.real());
    imag(imag()+z.imag());
    return *this;
  }

  /*! Subtracts a \p complex from this \p complex and
   *  assigns the result to this \p complex.
   *
   *  \param z The \p complex to be subtracted.
   */
  FMMTL_INLINE
  complex<T>& operator-=(const complex<T>& z) {
    real(real()-z.real());
    imag(imag()-z.imag());
    return *this;
  }

  /*! Multiplies this \p complex by another \p complex and
   *  assigns the result to this \p complex.
   *
   *  \param z The \p complex to be multiplied.
   */
  FMMTL_INLINE
  complex<T>& operator*=(const complex<T>& z) {
    return *this = *this * z;
  }

  /*! Divides this \p complex by another \p complex and
   *  assigns the result to this \p complex.
   *
   *  \param z The \p complex to be divided.
   */
  FMMTL_INLINE
  complex<T>& operator/=(const complex<T> z) {
    return *this = *this / z;
  }

  /* --- Getter functions ---
   * The volatile ones are there to help for example
   * with certain reductions optimizations
   */

  /*! Returns the real part of this \p complex.
   */
  FMMTL_INLINE const T& real() const volatile { return m_data[0]; }

  /*! Returns the imaginary part of this \p complex.
   */
  FMMTL_INLINE const T& imag() const volatile { return m_data[1]; }

  /*! Returns the real part of this \p complex.
   */
  FMMTL_INLINE const T& real() const { return m_data[0]; }

  /*! Returns the imaginary part of this \p complex.
   */
  FMMTL_INLINE const T& imag() const { return m_data[1]; }

  /* --- Setter functions ---
   * The volatile ones are there to help for example
   * with certain reductions optimizations
   */

  /*! Sets the real part of this \p complex.
   *
   *  \param re The new real part of this \p complex.
   */
  FMMTL_INLINE void real(const T& re) volatile { m_data[0] = re; }

  /*! Sets the imaginary part of this \p complex.
   *
   *  \param im The new imaginary part of this \p complex.e
   */
  FMMTL_INLINE void imag(const T& im) volatile { m_data[1] = im; }

  /*! Sets the real part of this \p complex.
   *
   *  \param re The new real part of this \p complex.
   */
  FMMTL_INLINE void real(const T& re) { m_data[0] = re; }

  /*! Sets the imaginary part of this \p complex.
   *
   *  \param im The new imaginary part of this \p complex.
   */
  FMMTL_INLINE void imag(const T& im) { m_data[1] = im; }

  /* --- Casting functions --- */

  /*! Casts this \p complex to a <tt>std::complex</tt> of the same type.
   */
  inline operator std::complex<T>() const {
    return std::complex<T>(real(),imag());
  }

 private:
  T m_data[2];
};


/* --- Equality Operators --- */

/*! Returns true if two \p complex numbers are equal and false otherwise.
 *
 *  \param lhs The first \p complex.
 *  \param rhs The second \p complex.
 */
template <typename T>
FMMTL_INLINE
bool operator==(const complex<T>& lhs, const complex<T>& rhs) {
  return lhs.real() == rhs.real() && lhs.imag() == rhs.imag();
}

/*! Returns true if the imaginary part of the  \p complex number is zero and the real part is equal to the scalar. Returns false otherwise.
 *
 *  \param lhs The scalar.
 *  \param rhs The \p complex.
 */
template <typename T>
FMMTL_INLINE
bool operator==(const T& lhs, const complex<T>& rhs) {
  return lhs == rhs.real() && T(0) == rhs.imag();
}

/*! Returns true if the imaginary part of the  \p complex number is zero and the real part is equal to the scalar. Returns false otherwise.
 *
 *  \param lhs The \p complex.
 *  \param rhs The scalar.
 */
template <typename T>
FMMTL_INLINE
bool operator==(const complex<T>& lhs, const T& rhs) {
  return lhs.real() == rhs && lhs.imag() == T(0);
}

/*! Returns true if two \p complex numbers are different and false otherwise.
 *
 *  \param lhs The first \p complex.
 *  \param rhs The second \p complex.
 */
template <typename T>
FMMTL_INLINE
bool operator!=(const complex<T>& lhs, const complex<T>& rhs) {
  return !(lhs == rhs);
}

/*! Returns true if the imaginary part of the  \p complex number is not zero or the real part is different from the scalar. Returns false otherwise.
 *
 *  \param lhs The scalar.
 *  \param rhs The \p complex.
 */
template <typename T>
FMMTL_INLINE
bool operator!=(const T& lhs, const complex<T>& rhs) {
  return !(lhs == rhs);
}

/*! Returns true if the imaginary part of the \p complex number is not zero or the real part is different from the scalar. Returns false otherwise.
 *
 *  \param lhs The \p complex.
 *  \param rhs The scalar.
 */
template <typename T>
FMMTL_INLINE
bool operator!=(const complex<T>& lhs, const T& rhs) {
  return !(lhs == rhs);
}


/* --- Unary Arithmetic operators --- */

/*! Unary plus, returns its \p complex argument.
 *
 *  \param rhs The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
const complex<T>& operator+(const complex<T>& rhs) {
  return rhs;
}

/*! Unary minus, returns the additive inverse (negation) of its \p complex argument.
 *
 *  \param rhs The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> operator-(const complex<T>& rhs) {
  return complex<T>(-rhs.real(), -rhs.imag());
}


/* --- Binary Arithmetic operators --- */

/*! Multiplies two \p complex numbers.
 *
 *  \param lhs The first \p complex.
 *  \param rhs The second \p complex.
 */
template <typename T>
FMMTL_INLINE
complex<T> operator*(const complex<T>& lhs, const complex<T>& rhs) {
  return complex<T>(lhs.real()*rhs.real() - lhs.imag()*rhs.imag(),
                    lhs.real()*rhs.imag() + lhs.imag()*rhs.real());
}

/*! Multiplies a \p complex number by a scalar.
 *
 *  \param lhs The \p complex.
 *  \param rhs The scalar.
 */
template <typename T>
FMMTL_INLINE
complex<T> operator*(const complex<T>& lhs, const T& rhs) {
  return complex<T>(lhs.real()*rhs, lhs.imag()*rhs);
}

/*! Multiplies a scalr by a \p complex number.
 *
 *  \param lhs The scalar.
 *  \param rhs The \p complex.
 */
template <typename T>
FMMTL_INLINE
complex<T> operator*(const T& lhs, const complex<T>& rhs) {
  return complex<T>(lhs*rhs.real(), lhs*rhs.imag());
}

/*! Divides two \p complex numbers.
 *
 *  \param lhs The numerator (dividend).
 *  \param rhs The denomimator (divisor).
 */
template <typename T>
FMMTL_INLINE
complex<T> operator/(const complex<T>& lhs, const complex<T>& rhs) {
  using std::abs;
  // XXX: Revisit
  T s = T(1) / (abs(rhs.real()) + abs(rhs.imag()));
  T ars = lhs.real() * s;
  T ais = lhs.imag() * s;
  T brs = rhs.real() * s;
  T bis = rhs.imag() * s;
  s = T(1) / ((brs * brs) + (bis * bis));
  return complex<T>(((ars * brs) + (ais * bis)) * s,
                    ((ais * brs) - (ars * bis)) * s);
}

/*! Divides a \p complex number by a scalar.
 *
 *  \param lhs The complex numerator (dividend).
 *  \param rhs The scalar denomimator (divisor).
 */
template <typename T>
FMMTL_INLINE
complex<T> operator/(const complex<T>& lhs, const T& rhs) {
  return complex<T>(lhs.real()/rhs, lhs.imag()/rhs);
}

/*! Divides a scalar by a \p complex number.
 *
 *  \param lhs The scalar numerator (dividend).
 *  \param rhs The complex denomimator (divisor).
 */
template <typename T>
FMMTL_INLINE
complex<T> operator/(const T& lhs, const complex<T>& rhs) {
  using std::abs;
  // XXX: Revisit
  T s = T(1) / (abs(rhs.real()) + abs(rhs.imag()));
  T ars = lhs.real() * s;
  T brs = rhs.real() * s;
  T bis = rhs.imag() * s;
  s = T(1) / ((brs * brs) + (bis * bis));
  return complex<T>((ars * brs) *  s,
                    (ars * bis) * -s);
}

/*! Adds two \p complex numbers.
 *
 *  \param lhs The first \p complex.
 *  \param rhs The second \p complex.
 */
template <typename T>
FMMTL_INLINE
complex<T> operator+(const complex<T>& lhs, const complex<T>& rhs) {
  return complex<T>(lhs.real()+rhs.real(), lhs.imag()+rhs.imag());
}

/*! Adds a scalar to a \p complex number.
 *
 *  \param lhs The \p complex.
 *  \param rhs The scalar.
 */
template <typename T>
FMMTL_INLINE
complex<T> operator+(const complex<T>& lhs, const T& rhs) {
  return complex<T>(lhs.real()+rhs, lhs.imag());
}

/*! Adds a \p complex number to a scalar.
 *
 *  \param lhs The scalar.
 *  \param rhs The \p complex.
 */
template <typename T>
FMMTL_INLINE
complex<T> operator+(const T& lhs, const complex<T>& rhs) {
  return complex<T>(lhs+rhs.real(), rhs.imag());
}

/*! Subtracts two \p complex numbers.
 *
 *  \param lhs The first \p complex (minuend).
 *  \param rhs The second \p complex (subtrahend).
 */
template <typename T>
FMMTL_INLINE
complex<T> operator-(const complex<T>& lhs, const complex<T>& rhs) {
  return complex<T>(lhs.real()-rhs.real(), lhs.imag()-rhs.imag());
}

/*! Subtracts a scalar from a \p complex number.
 *
 *  \param lhs The \p complex (minuend).
 *  \param rhs The scalar (subtrahend).
 */
template <typename T>
FMMTL_INLINE
complex<T> operator-(const complex<T>& lhs, const T& rhs) {
  return complex<T>(lhs.real()-rhs, lhs.imag());
}

/*! Subtracts a \p complex number from a scalar.
 *
 *  \param lhs The scalar (minuend).
 *  \param rhs The \p complex (subtrahend).
 */
template <typename T>
FMMTL_INLINE
complex<T> operator-(const T& lhs, const complex<T>& rhs) {
  return complex<T>(lhs-rhs.real(), -rhs.imag());
}



/* --- General Functions --- */

/*! Returns the magnitude (also known as absolute value) of a \p complex.
 *
 *  \param z The \p complex from which to calculate the absolute value.
 */
template <typename T>
FMMTL_INLINE
T abs(const complex<T>& z) {
  using std::sqrt;
  return sqrt(norm(z));   // XXX std::hypot on C++11
}

/*! Returns the phase angle (also known as argument) in radians of a \p complex.
 *
 *  \param z The \p complex from which to calculate the phase angle.
 */
template <typename T>
FMMTL_INLINE
T arg(const complex<T>& z) {
  using std::atan2;
  return atan2(z.imag(), z.real());
}

/*! Returns the square of the magnitude of a \p complex.
 *
 *  \param z The \p complex from which to calculate the norm.
 */
template <typename T>
FMMTL_INLINE
T norm(const complex<T>& z) {
  // XXX: Under/Overflow?
  return z.real()*z.real() + z.imag()*z.imag();
}

/*! Returns the complex conjugate of a \p complex.
 *
 *  \param z The \p complex from which to calculate the complex conjugate.
 */
template <typename T>
FMMTL_INLINE
complex<T> conj(const complex<T>& z) {
  return complex<T>(z.real(), -z.imag());
}

/*! Returns a \p complex with the specified magnitude and phase.
 *
 *  \param m The magnitude of the returned \p complex.
 *  \param theta The phase of the returned \p complex in radians.
 */
template <typename T>
FMMTL_INLINE
complex<T> polar(const T& m, const T& theta = 0) {
  using std::cos;
  using std::sin;
  return complex<T>(m * cos(theta), m * sin(theta));
}

/*! Returns the projection of a \p complex on the Riemann sphere.
 *  For all finite \p complex it returns the argument. For \p complexs
 *  with a non finite part returns (INFINITY,+/-0) where the sign of
 *  the zero matches the sign of the imaginary part of the argument.
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> proj(const T& z);


/* --- Exponential Functions --- */

/*! Returns the complex exponential of a \p complex number.
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> exp(const complex<T>& z) {
  // XXX: Under/Overflow
  using std::exp;
  return polar(exp(z.real()), z.imag());
}

/*! Returns the complex natural logarithm of a \p complex number.
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> log(const complex<T>& z) {
  using std::log;
  return complex<T>(log(abs(z)), arg(z));
}

/*! Returns the complex base 10 logarithm of a \p complex number.
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> log10(const complex<T>& z) {
  return log(z) / T(2.30258509299404568402);
}


/* --- Power Functions --- */

/*! Returns a \p complex number raised to another.
 *
 *  \param x The base.
 *  \param y The exponent.
 */
template <typename T>
FMMTL_INLINE
complex<T> pow(const complex<T>& x, const complex<T>& y) {
  return exp(log(x) * y);
}

/*! Returns a \p complex number raised to a scalar.
 *
 *  \param x The \p complex base.
 *  \param y The scalar exponent.
 */
template <typename T>
FMMTL_INLINE
complex<T> pow(const complex<T>& x, const T& y) {
  return exp(log(x) * y);
}

/*! Returns a scalar raised to a \p complex number.
 *
 *  \param x The scalar base.
 *  \param y The \p complex exponent.
 */
template <typename T>
FMMTL_INLINE
complex<T> pow(const T& x, const complex<T>& y) {
  using std::log;
  return exp(log(x) * y);
}

/*! Returns the complex square root of a \p complex number.
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> sqrt(const complex<T>& z) {
  using std::sqrt;
  return polar(sqrt(abs(z)), arg(z)/T(2));
}


/* --- Trigonometric Functions --- */

/*! Returns the complex cosine of a \p complex number.
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> cos(const complex<T>& z) {
  using std::sin; using std::sinh;
  using std::cos; using std::cosh;
  return complex<T>( cos(z.real()) * cosh(z.imag()),
                    -sin(z.real()) * sinh(z.imag()));
}

/*! Returns the complex sine of a \p complex number.
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> sin(const complex<T>& z) {
  using std::sin; using std::sinh;
  using std::cos; using std::cosh;
  return complex<T>(sin(z.real()) * cosh(z.imag()),
                    cos(z.real()) * sinh(z.imag()));
}

/*! Returns the complex tangent of a \p complex number.
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> tan(const complex<T>& z) {
  return sin(z) / cos(z);
}


/* --- Hyperbolic Functions --- */

/*! Returns the complex hyperbolic cosine of a \p complex number.
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> cosh(const complex<T>& z) {
  using std::sin; using std::sinh;
  using std::cos; using std::cosh;
  return complex<T>(cosh(z.real()) * cos(z.imag()),
                    sinh(z.real()) * sin(z.imag()));
}

/*! Returns the complex hyperbolic sine of a \p complex number.
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> sinh(const complex<T>& z) {
  using std::sin; using std::sinh;
  using std::cos; using std::cosh;
  return complex<T>(sinh(z.real()) * cos(z.imag()),
                    cosh(z.real()) * sin(z.imag()));
}

/*! Returns the complex hyperbolic tangent of a \p complex number.
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> tanh(const complex<T>& z) {
  const complex<T> r = exp(T(2) * z);
  return (r - T(1)) / (r + T(1));
}


/* --- Inverse Trigonometric Functions --- */

/*! Returns the complex arc cosine of a \p complex number.
 *
 *  The range of the real part of the result is [0, Pi] and
 *  the range of the imaginary part is [-inf, +inf]
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> acos(const complex<T>& z) {
  const complex<T> ret = asin(z);
  const T pi_half = T(3.1415926535897932384626433832795) / T(2);
  return complex<T>(pi_half - ret.real(), -ret.imag());
}

/*! Returns the complex arc sine of a \p complex number.
 *
 *  The range of the real part of the result is [-Pi/2, Pi/2] and
 *  the range of the imaginary part is [-inf, +inf]
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> asin(const complex<T>& z) {
  // asin(z) = -i*asinh(i*z)
  const complex<T> ret = asinh(complex<T>(-z.imag(), z.real()));
  return complex<T>(ret.imag(), -ret.real());
}

/*! Returns the complex arc tangent of a \p complex number.
 *
 *  The range of the real part of the result is [-Pi/2, Pi/2] and
 *  the range of the imaginary part is [-inf, +inf]
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> atan(const complex<T>& z) {
  // asin(z) = -i*atanh(i*z)
  const complex<T> ret = atanh(complex<T>(-z.imag(), z.real()));
  return complex<T>(ret.imag(), -ret.real());
}


/* --- Inverse Hyperbolic Functions --- */

/*! Returns the complex inverse hyperbolic cosine of a \p complex number.
 *
 *  The range of the real part of the result is [0, +inf] and
 *  the range of the imaginary part is [-Pi, Pi]
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> acosh(const complex<T>& z) {
  complex<T> ret((z.real()-z.imag()) * (z.real()+z.imag()) - T(1),
                 T(2) * z.real() * z.imag());
  ret = sqrt(ret);
  if (z.real() < T(0))
    ret = -ret;
  ret += z;
  ret = log(ret);
  if (ret.real() < T(0))
    ret = -ret;
  return ret;
}

/*! Returns the complex inverse hyperbolic sine of a \p complex number.
 *
 *  The range of the real part of the result is [-inf, +inf] and
 *  the range of the imaginary part is [-Pi/2, Pi/2]
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> asinh(const complex<T>& z) {
  return log(sqrt(z*z+T(1.0))+z);
}

/*! Returns the complex inverse hyperbolic tangent of a \p complex number.
 *
 *  The range of the real part of the result is [-inf, +inf] and
 *  the range of the imaginary part is [-Pi/2, Pi/2]
 *
 *  \param z The \p complex argument.
 */
template <typename T>
FMMTL_INLINE
complex<T> atanh(const complex<T>& z) {
  T imag2 = z.imag()*z.imag();
  T n = T(1) + z.real();
  n = imag2 + n*n;

  T d = T(1) - z.real();
  d = imag2 + d*d;
  d = T(1) - z.real()*z.real() - imag2;
  using std::log;
  using std::atan2;
  return complex<T>(T(0.25)*(log(n) - log(d)), T(0.5)*atan2(T(2)*z.imag(),d));
}


/* --- Stream Operators --- */

/*! Writes to an output stream a \p complex number in the form (real,imaginary).
 *
 *  \param os The output stream.
 *  \param z The \p complex number to output.
 */
template <typename T, class charT, class traits>
std::basic_ostream<charT,traits>&
operator<<(std::basic_ostream<charT,traits>& os, const complex<T>& z) {
  return os << '(' << z.real() << ',' << z.imag() << ')';
}

/*! Reads a \p complex number from an input stream.
 *  The recognized formats are:
 * - real
 * - (real)
 * - (real, imaginary)
 *
 * The values read must be convertible to the \p complex's \c value_type
 *
 *  \param is The input stream.
 *  \param z The \p complex number to set.
 */
template <typename T, typename charT, class traits>
std::basic_istream<charT,traits>&
operator>>(std::basic_istream<charT,traits>& is, complex<T>& z) {
  T re, im;

  charT ch;
  is >> ch;

  if(ch == '(') {
    is >> re >> ch;
    if (ch == ',') {
      is >> im >> ch;
      if (ch == ')') {
	      z = complex<T>(re, im);
	    } else {
	      is.setstate(std::ios_base::failbit);
	    }
    } else if (ch == ')') {
      z = re;
    } else {
      is.setstate(std::ios_base::failbit);
    }
  } else {
    is.putback(ch);
    is >> re;
    z = re;
  }
  return is;
}

} // end namespace fmmtl
