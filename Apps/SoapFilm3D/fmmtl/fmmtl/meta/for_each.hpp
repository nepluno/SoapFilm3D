#pragma once
/** @file for_each.hpp
 * @brief A short implementation of boost::mpl::for_each.
 * This is optimized for default item construction and uses a non-recursive
 * implementation for integer_sequences (avoid template recursion depth errors).
 *
 * TODO: Improve by sequence dispatching on static iterator traversal
 */

#include "integer_sequence.hpp"

namespace fmmtl {

// Declarations
template <typename T>
struct size;
template <typename T>
struct begin;
template <typename T>
struct end;
template <typename T>
struct deref;
template <typename T>
struct next;

template <typename T>
using size_t = typename size<T>::type;
template <typename T>
using begin_t = typename begin<T>::type;
template <typename T>
using end_t = typename end<T>::type;
template <typename T>
using deref_t = typename deref<T>::type;
template <typename T>
using next_t = typename next<T>::type;

// for_each
template <typename First, typename Last, typename F>
constexpr
typename std::enable_if<std::is_same<First,Last>::value>::type
for_each_impl(F&&) {}

template <typename First, typename Last, typename F>
constexpr
typename std::enable_if<!std::is_same<First,Last>::value>::type
for_each_impl(F&& f) {
  using i    = typename fmmtl::deref<First>::type;
  using Next = typename fmmtl::next<First>::type;

  std::forward<F>(f)(i{}), void(), for_each_impl<Next,Last>(std::forward<F>(f));
}

/** Static for_each over a static forward iterator range.
 * @tparam ...
 * @tparam F models a function object with signature
 *    ignored-return f(typename deref<First>::type{})
 *
 *
 * Applies @a f to a default constructed type in the range [First, Last).
 */
template <typename First, typename Last, typename F>
constexpr
void
for_each(First, Last, F&& f) {
  return for_each_impl<First,Last>(std::forward<F>(f));
}

/** Static for_each for a forward sequence.
 * @tparam Sequence models a static forward sequence with:
 *  typename fmmtl::begin<Sequence>::type   --  Begin static iterator
 *  typename fmmtl::end<Sequence>::type     --  End static iterator
 *  typename fmmtl::deref<Iterator>::type   --  Dereferenced static iterator
 *  typename fmmtl::next<Iterator>::type    --  Increment static iterator
 * @tparam F models a function object with signature
 *    ignored-return f(typename deref<typename begin<Sequence>::type>::type{})
 *
 * Applies @f to a default-constructed type in the range [Begin, End).
 */
template <typename Sequence, typename F>
constexpr
void
for_each(Sequence, F&& f) {
  using First = typename fmmtl::begin<Sequence>::type;
  using Last  = typename fmmtl::end<Sequence>::type;

  return for_each_impl<First,Last>(std::forward<F>(f));
}

/** Static for_each for an integer_sequence.
 *
 * Applies @f to a std::integral_constant<T,i> for each i in I...
 */
template <typename T, T... I, typename F>
constexpr
void
for_each(integer_sequence<T,I...>, F&& f) {
  using eat = int[];
  (void) eat {0,(void(std::forward<F>(f)(std::integral_constant<T,I>{})),0)...};
}

}
