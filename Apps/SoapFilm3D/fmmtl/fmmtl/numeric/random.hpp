#pragma once

#include <limits>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_same.hpp>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/numeric/Complex.hpp"

namespace fmmtl {

static boost::random::mt19937 default_generator;

using boost::enable_if;
using boost::is_integral;
using boost::is_floating_point;

template <typename T, class Enable = void>
struct random;

template <typename T>
struct random<T, typename enable_if<is_integral<T> >::type> {
  typedef T result_type;
  result_type operator()(T a, T b) const { return get(a,b); }
  result_type operator()()         const { return get();    }

  static result_type get(T a, T b) {
    boost::random::uniform_int_distribution<T> dist(a, b);
    return dist(default_generator);
  }
  static result_type get() {
    return get(T(0), std::numeric_limits<T>::max());
  }
};

template <typename T>
struct random<T, typename enable_if<is_floating_point<T> >::type> {
  typedef T result_type;
  result_type operator()(T a, T b) const { return get(a,b); }
  result_type operator()()         const { return get();    }

  static result_type get(T a, T b) {
    boost::random::uniform_real_distribution<T> dist(a, b);
    return dist(default_generator);
  }
  static result_type get() {
    return get(T(0), T(1));
  }
};

template <typename T>
struct random<complex<T> > {
  typedef complex<T> result_type;
  result_type operator()(T a, T b) const { return get(a,b); }
  result_type operator()()         const { return get();    }

  static result_type get(T a, T b) {
    return complex<T>(random<T>::get(a,b), random<T>::get(a,b));
  }
  static result_type get() {
    return get(T(0), T(1));
  }
};

template <typename T>
struct random<std::complex<T> > {
  typedef std::complex<T> result_type;
  result_type operator()(T a, T b) const { return get(a,b); }
  result_type operator()()         const { return get();    }

  static result_type get(T a, T b) {
    return std::complex<T>(random<T>::get(a,b), random<T>::get(a,b));
  }
  static result_type get() {
    return get(T(0), T(1));
  }
};

template <std::size_t N, typename T>
struct random<Vec<N,T> > {
  typedef Vec<N,T> result_type;
  result_type operator()(T a, T b) const { return get(a,b); }
  result_type operator()()         const { return get();    }

  static result_type get(T a, T b) {
    Vec<N,T> v;
    for (std::size_t i = 0; i != N; ++i)
      v[i] = fmmtl::random<T>::get(a, b);
    return v;
  }
  static result_type get() {
    return get(T(0), T(1));
  }
};


class random_n {
  template <typename T>
  struct generator {
    T operator()(const std::size_t&) const { return random<T>::get(); }
  };
  std::size_t N;

 public:
  random_n(const std::size_t& _N) : N(_N) {}

  template <typename Container>
  operator Container() const {
    typedef typename Container::value_type value_type;
    return Container(
        boost::make_transform_iterator(
            boost::make_counting_iterator(std::size_t(0)),
            generator<value_type>()),
        boost::make_transform_iterator(
            boost::make_counting_iterator(std::size_t(N)),
            generator<value_type>()));
  }
};


} // end namespace fmmtl
