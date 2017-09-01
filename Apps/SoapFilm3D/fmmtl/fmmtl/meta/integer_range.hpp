#pragma once

/** Index range generation
 * TODO: Use C++14 std::integer_sequence
 */
template <class T, T... Ints>
struct integer_sequence {};
template <std::size_t... Ints>
using index_sequence = integer_sequence<std::size_t, Ints...>;

template <class T, std::size_t N, T... Is>
struct generate_integer_sequence {
  using type = typename generate_integer_sequence<T, N-1, N-1, Is...>::type;
};
template <class T, T... Is>
struct generate_integer_sequence<T, 0, Is...> {
  using type = integer_sequence<T, Is...>;
};

template <class T, T N>
using make_integer_sequence = typename generate_integer_sequence<T, N>::type;
template <std::size_t N>
using make_index_sequence = make_integer_sequence<std::size_t, N>;
