#pragma once

#include <type_traits>

namespace fmmtl {

///////////
// Aliases
template <bool N>
using bool_t = std::integral_constant<bool, N>;
template <char N>
using char_t = std::integral_constant<char, N>;
template <unsigned char N>
using uchar_t = std::integral_constant<unsigned char, N>;
template <short N>
using short_t = std::integral_constant<short, N>;
template <unsigned short N>
using ushort_t = std::integral_constant<unsigned short, N>;
template <int N>
using int_t = std::integral_constant<int, N>;
template <unsigned int N>
using uint_t = std::integral_constant<unsigned int, N>;
template <long N>
using long_t = std::integral_constant<long, N>;
template <unsigned long N>
using ulong_t = std::integral_constant<unsigned long, N>;
template <long long N>
using longlong_t = std::integral_constant<long long, N>;
template <unsigned long long N>
using ulonglong_t = std::integral_constant<unsigned long long, N>;

template <std::size_t N>
using index_t = std::integral_constant<std::size_t, N>;


// TODO: Use C++14 std::integer_sequence
template <typename T, T... N>
struct integer_sequence {
  using type = integer_sequence;
  using value_type = T;
  static constexpr std::size_t size() noexcept {
    return sizeof...(N);
  }
};

template <std::size_t... N>
using index_sequence = integer_sequence<std::size_t, N...>;

//////////////////
// make_*
namespace detail {

// Glue two sets of integer_sequence together
template <typename I1, typename I2>
struct glue_integer_sequence;

template <typename T, T... N1, T... N2>
struct glue_integer_sequence<integer_sequence<T, N1...>,
                             integer_sequence<T, N2...>> {
  using type = integer_sequence<T, N1..., (sizeof...(N1) + N2)...>;
};

template <typename T, std::size_t N>
struct make_integer_sequence_
    : glue_integer_sequence<typename make_integer_sequence_<T,    N/2>::type,
                            typename make_integer_sequence_<T,N - N/2>::type> {
};

template <typename T>
struct make_integer_sequence_<T, 0> {
  using type = integer_sequence<T>;
};

template <typename T>
struct make_integer_sequence_<T, 1> {
  using type = integer_sequence<T, 0>;
};

template <typename T, typename Seq, T Begin>
struct make_integer_range_;

template <typename T, T... Is, T Begin>
struct make_integer_range_<T, integer_sequence<T, Is...>, Begin> {
  using type = integer_sequence<T, Begin + Is...>;
};

} // end namespace detail


// Generate integer_sequence [0,N) in O(log(N)) time
template <typename T, T N>
using make_integer_sequence = typename detail::make_integer_sequence_<T,N>::type;

template <std::size_t N>
using make_index_sequence = make_integer_sequence<std::size_t, N>;

/* Similar to make_integer_sequence<>, except it goes from [Begin, End) */
// TODO: Strided and O(log(End-Begin))
template <typename T, T Begin, T End>
using make_integer_range =
    typename detail::make_integer_range_<T,
                                         make_integer_sequence<T, End-Begin>,
                                         Begin>::type;

/* Similar to make_index_sequence, except it goes from [Begin, End) */
template <std::size_t Begin, std::size_t End>
using make_index_range = make_integer_range<std::size_t, Begin, End>;

///////////////
// Utilities //
///////////////

///////////////
// intseq_size
template <typename Seq>
struct intseq_size {};

template <typename T, T... rest>
struct intseq_size<integer_sequence<T, rest...>>
    : std::integral_constant<std::size_t, sizeof...(rest)> {
};

//////////////
// intseq_cat
template <typename... Seqs>
struct intseq_cat;

template <>
struct intseq_cat<> {
  using type = integer_sequence<std::size_t>;
};

template <typename T, T... N1>
struct intseq_cat<integer_sequence<T, N1...>> {
  using type = integer_sequence<T, N1...>;
};

template <typename T, T... N1, T... N2>
struct intseq_cat<integer_sequence<T, N1...>,
                  integer_sequence<T, N2...>> {
  using type = integer_sequence<T, N1..., N2...>;
};

template <typename T, T... N1, T... N2, typename... Seqs>
struct intseq_cat<integer_sequence<T, N1...>,
                  integer_sequence<T, N2...>,
                  Seqs...>
    : intseq_cat<integer_sequence<T, N1..., N2...>, Seqs...> {
};

template <typename... Seq>
using intseq_cat_t = typename intseq_cat<Seq...>::type;

///////////////
// intseq_repeat
template <std::size_t N, typename T = std::size_t, T i = 0>
struct intseq_repeat
    : intseq_cat<typename intseq_repeat<    N/2, T, i>::type,
                 typename intseq_repeat<N - N/2, T, i>::type> {
};

template <typename T, T i>
struct intseq_repeat<0, T, i> {
  using type = integer_sequence<T>;
};

template <typename T, T i>
struct intseq_repeat<1, T, i> {
  using type = integer_sequence<T, i>;
};

template <std::size_t N, typename T = std::size_t, T i = 0>
using intseq_repeat_t = typename intseq_repeat<N, T, i>::type;

template <std::size_t N, std::size_t i = 0>
using idxseq_repeat_t = typename intseq_repeat<N, std::size_t, i>::type;

/////////////////////////////////////////////
// intseq_element
namespace detail {

template <std::size_t N, typename S = make_index_sequence<N>>
struct at;

template <std::size_t N, std::size_t... ignore>
struct at<N, index_sequence<ignore...>> {
  template <typename Nth>
  static constexpr Nth apply(decltype(ignore,(void*)0)..., Nth nth, ...)
  { return nth; }
};

} // end namespace detail

template <std::size_t N, typename... Xs>
constexpr auto at(Xs... xs) -> decltype(*detail::at<N>::apply(&xs...)) {
  return *detail::at<N>::apply(&xs...);
}

template <std::size_t N, typename Seq>
struct intseq_element {};

template <std::size_t N, typename T, T... I>
struct intseq_element<N, integer_sequence<T, I...>> {
  using type = std::integral_constant<T, at<N>(I...)>;
};

template<std::size_t N, typename Seq>
using intseq_element_t = typename intseq_element<N, Seq>::type;

////////////////
// intseq_front
template <typename Seq>
struct intseq_front {};

template <typename T, T i, T... rest>
struct intseq_front<integer_sequence<T, i, rest...>>
    : std::integral_constant<T, i> {
};

/////////////////////
// intseq_push_front
template <typename Seq, typename Int>
struct intseq_push_front {};

template <typename T, T... rest, typename U, U i>
struct intseq_push_front<integer_sequence<T, rest...>,
                         std::integral_constant<U, i>> {
  using type = integer_sequence<T, i, rest...>;
};

template <typename Seq, typename Int>
using intseq_push_front_t = typename intseq_push_front<Seq, Int>::type;

////////////////////
// intseq_pop_front
template <typename Seq>
struct intseq_pop_front {};

template <typename T, T i, T... rest>
struct intseq_pop_front<integer_sequence<T, i, rest...>> {
  using type = integer_sequence<T, rest...>;
};

template <typename Seq>
using intseq_pop_front_t = typename intseq_pop_front<Seq>::type;

///////////////
// intseq_back
template <typename Seq>
struct intseq_back {};

template <typename T, T i, T... rest>
struct intseq_back<integer_sequence<T, i, rest...>>
    : intseq_element<sizeof...(rest), integer_sequence<T, i, rest...>> {
};

////////////////////
// intseq_push_back
template <typename Seq, typename Int>
struct intseq_push_back {};

template <typename T, T... rest, typename U, U i>
struct intseq_push_back<integer_sequence<T, rest...>,
                        std::integral_constant<U, i>> {
  using type = integer_sequence<T, rest..., i>;
};

template<typename Seq, typename Int>
using intseq_push_back_t = typename intseq_push_back<Seq, Int>::type;

///////////////////
// intseq_pop_back
// Not provided.
// It cannot be made to meet the complexity guarantees one would expect.

///////////////////
// intseq_is_empty
template <typename Seq>
struct intseq_is_empty : std::false_type {
};

template <typename T>
struct intseq_is_empty<integer_sequence<T>> : std::true_type {
};

} // end namespace fmmtl
