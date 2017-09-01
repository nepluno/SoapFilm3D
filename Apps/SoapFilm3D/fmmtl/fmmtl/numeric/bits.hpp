#pragma once

#include <type_traits>
#include <limits>

namespace fmmtl {

/** Smears the bits in c into the low bits by steps of N
 *
 * Example: N = 3
 *          0000010000000000 -> 0000010010010010
 */
template <unsigned N, typename T>
inline constexpr T smear_low_(T c) {
  static_assert(std::is_unsigned<T>::value, "Smear requires unsigned type");
  return 0 < N && N < std::numeric_limits<T>::digits
      ? smear_low_<2*N,T>(c | (c >> N))
      : c;
  /* Equivalent to (unrolled)
     for (unsigned i = N; i < std::numeric_limits<T>::digits; i *= 2)
       c |= (c >> i);
     return c;
   */
}

/** Smears the bits in c into the high bits by steps of N
 *
 * Example: N = 3
 *          0000000000001000 -> 1001001001001000
 */
template <unsigned N, typename T>
inline constexpr T smear_high_(T c) {
  static_assert(std::is_unsigned<T>::value, "Smear requires unsigned type");
  return 0 < N && N < std::numeric_limits<T>::digits
      ? smear_high_<2*N,T>(c | (c << N))
      : c;
  /* Equivalent to (unrolled)
     for (unsigned i = N; i < std::numeric_limits<T>::digits; i *= 2)
       c |= (c << i);
     return c;
  */
}

/** Returns the smallest power of two greater than the input value
 */
template <typename T>
inline constexpr T next_pow_2(T c) {
  return smear_low_<1>(c) + 1;
}

/** Returns the smallest power of two greater than or equal to the input value
 * @note Undefined for 0.
 */
template <typename T>
inline constexpr T ceil_pow_2(T c) {
  return smear_low_<1>(c - 1) + 1;
}

/** Returns the greatest power of two less than the input value
 * @note Undefined for 0 and 1.
 */
template <typename T>
inline constexpr T prev_pow_2(T c) {
  return (smear_low_<1>(c - 1) >> 1) + 1;
}

/** Returns the greatest power of two less than or equal to the input value.
 * @note Undefined for 0.
 */
template <typename T>
inline constexpr T floor_pow_2(T c) {
  return (smear_low_<1>(c) >> 1) + 1;
}

/** Recursive template helper function for spread_bits_
 */
template <unsigned N,            // Zeros stride
          typename T,            // Underlying unsigned type
          unsigned I>            // Iteration counter (Group size)
inline constexpr T spread_helper_(T c) {
  // c = (c | (c << (I*(N-1)))) & 0b(I ones and I*(N-1) zeros repeating)
  return I > 0
      ? spread_helper_<N,T,I/2>((c | (c << (I*N-I))) & smear_high_<I*N>((T(1) << I)-1))
      : c;
}

/** Inserts N-1 zeros in between the first B/N bits of c, where B is the number
 * of bits in the type T.
 *
 * Example: N = 3
 *        0000000000101101 -> 1000001001000001
 *        0000000000****** -> *00*00*00*00*00*
 *
 * @tparam T models an unsigned integer type.
 * @pre 0 <= c < std::pow(2, B/N), where B = std::numeric_limits<T>::digits
 */
template <unsigned N, typename T>
inline constexpr T spread_bits_(T c) {
  static_assert(std::is_unsigned<T>::value, "Spread requires unsigned type");
  return N > 1
      ? spread_helper_<N,T,prev_pow_2(std::numeric_limits<T>::digits/N)>(c)
      : c;
}

template <unsigned N,             // Zeros stride
          typename T,             // Underlying type
          unsigned I>             // Iteration counter
inline constexpr T compact_helper_(T c) {
  return I*N < std::numeric_limits<T>::digits
      ? compact_helper_<N,T,2*I>((c | (c >> (I*N-I))) & smear_high_<2*I*N>((T(1) << (2*I))-1))
      : c;
}

/** Compacts every Nth bit into the first B/N bits of c, where B is the number
 * of bits in the type T.
 *
 * Example: N = 3
 *        1000001001000001 -> 0000000000101101
 *        *00*00*00*00*00* -> 0000000000******
 *
 * @tparam T models an unsigned integer type
 */
template <unsigned N, typename T>
inline constexpr T compact_bits_(T c) {
  static_assert(std::is_unsigned<T>::value, "Compact requires unsigned type");
  return N > 1
      ? compact_helper_<N,T,1>(c & smear_high_<N>(T(1)))
      : c;
}

} // end namespace fmmtl
