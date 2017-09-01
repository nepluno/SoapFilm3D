#pragma once
/** @file T2P.hpp
 * @brief Dispatch methods for the T2P stage
 *
 */

#include <type_traits>

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

/** Default behavior is no-op, but check for convertibility at compile-time */
template <typename Expansion, bool _ = ExpansionTraits<Expansion>::has_T2P>
struct T2P_Functor;

/** No T2P operator, default to no-op */
template <typename E>
struct T2P_Functor<E, false> {
  T2P_Functor(const E&) {}
  inline const typename ExpansionTraits<E>::target_type&
  operator()(const typename ExpansionTraits<E>::target_type& target) const {
    // No T2P operator, make sure the types are convertible
    FMMTL_STATIC_ASSERT(
        std::is_convertible<typename ExpansionTraits<E>::target_type,
                            typename ExpansionTraits<E>::point_type>::value,
        "target_type is not convertible to point_type and no T2P provided!");

    return target;
  }
};

template <typename E>
struct T2P_Functor<E, true> {
  T2P_Functor(const E& e) : e_(e) {}
  inline typename ExpansionTraits<E>::point_type
  operator()(const typename ExpansionTraits<E>::target_type& target) const {
    return e_.T2P(target);
  }
 private:
  const E& e_;
};

template <typename Expansion>
inline T2P_Functor<Expansion> T2P(const Expansion& e) {
  return T2P_Functor<Expansion>(e);
}
