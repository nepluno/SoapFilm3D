#pragma once
/** @file math.hpp
 * A collection of constexpr math operations.
 */

namespace fmmtl {

// Greatest common divisor
template <typename T>
constexpr T gcd(const T& x, const T& y) {
  return ((x%y) == 0) ? y : gcd(y,x%y);
}

// Sum
template <typename T>
constexpr const T& sum(const T& n) {
  return n;
}
template <typename T, typename... More>
constexpr T sum(const T& n, More... ns) {
  return n + sum(ns...);
}

// Product
template <typename T>
constexpr const T& product(const T& n) {
  return n;
}
template <typename T, typename... More>
constexpr T product(const T& n, More... ns) {
  return n * product(ns...);
}

// Factorial
// n! = n * (n-1) * (n-2) * ... * 1
constexpr double factorial(std::size_t n) {
	return (n <= 1) ? 1 : n * factorial(n-1);
};

// Combinatorial Permutation
// (n)_k = n! / (n-k)! = n * (n-1) * (n-2) * ... * (n-k+1)
constexpr double permutation(std::size_t n, std::size_t k) {
	return (k == 0) ? 1 : n * permutation(n-1,k-1);
};

// Combinatorial Combination -- Binomial Coefficient
// n choose k = n! / (k! (n-k)!)
constexpr double combination_impl(std::size_t n, std::size_t k) {
  return (k == 0) ? 1 : n * combination_impl(n-1,k-1) / k;
}
constexpr double combination(std::size_t n, std::size_t k) {
  return (n/2 > k) ? combination_impl(n,n-k) : combination_impl(n,k);
};

} // end namespace fmmtl
