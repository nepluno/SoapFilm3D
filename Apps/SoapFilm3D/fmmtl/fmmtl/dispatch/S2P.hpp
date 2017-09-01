#pragma once
/** @file S2P.hpp
 * @brief Dispatch methods for the S2P stage
 *
 */

#include <type_traits>

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

/** Default behavior is no-op, but check for convertibility at compile-time */
template <typename Expansion, bool _ = ExpansionTraits<Expansion>::has_S2P>
struct S2P_Functor;

/** No S2P operator, default to no-op */
template <typename E>
struct S2P_Functor<E, false> {
  S2P_Functor(const E&) {}
  inline const typename ExpansionTraits<E>::source_type&
  operator()(const typename ExpansionTraits<E>::source_type& source) const {
    // No S2P operator, make sure the types are convertible
    FMMTL_STATIC_ASSERT(
        std::is_convertible<typename ExpansionTraits<E>::source_type,
                            typename ExpansionTraits<E>::point_type>::value,
        "source_type is not convertible to point_type and no S2P provided!");

    return source;
  }
};

template <typename E>
struct S2P_Functor<E, true> {
  S2P_Functor(const E& e) : e_(e) {}
  inline typename ExpansionTraits<E>::point_type
  operator()(const typename ExpansionTraits<E>::source_type& source) const {
    return e_.S2P(source);
  }
 private:
  const E& e_;
};

template <typename Expansion>
inline S2P_Functor<Expansion> S2P(const Expansion& e) {
  return S2P_Functor<Expansion>(e);
}
