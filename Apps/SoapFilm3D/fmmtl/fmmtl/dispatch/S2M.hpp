#pragma once
/** @file S2M.hpp
 * @brief Dispatch methods for the S2M stage
 *
 */

#include <type_traits>

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

/** Default behavior gives a warning -- using non-existent method */
template <bool has_s2m>
struct S2M_Helper {
  template <typename... Args>
  inline static void apply(Args&&...) {
    std::cerr << "WARNING: Expansion does not have a correct S2M!\n";
  }
  template <typename... Args>
  inline static void eval(Args&&...) {
    apply();
  }
};

/** Expansion has an S2M method (scalar/vector) to dispatch to */
template <>
struct S2M_Helper<true> {
  /** The Expansion provides a vector S2M accumulator. */
  template <typename Expansion, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_vector_S2M>::type
  apply(const Expansion& K,
        SourceIter s_begin, SourceIter s_end, ChargeIter c_begin,
        const typename Expansion::point_type& center,
        typename Expansion::multipole_type& M) {
    K.S2M(s_begin, s_end, c_begin, center, M);
  }

  /** The Expansion provides a scalar S2M accumulator. */
  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_S2M &
                          !ExpansionTraits<Expansion>::has_vector_S2M>::type
  apply(const Expansion& K,
        const typename Expansion::source_type& source,
        const typename Expansion::charge_type& charge,
        const typename Expansion::point_type& center,
        typename Expansion::multipole_type& M) {
    K.S2M(source, charge, center, M);
  }

  /** The Expansion provides a scalar S2M accumulator. */
  template <typename Expansion, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_S2M &
                          !ExpansionTraits<Expansion>::has_vector_S2M>::type
  apply(const Expansion& K,
        SourceIter s_begin, SourceIter s_end, ChargeIter c_begin,
        const typename Expansion::point_type& center,
        typename Expansion::multipole_type& M) {
    for ( ; s_begin != s_end; ++s_begin, ++c_begin)
      apply(K, *s_begin, *c_begin, center, M);
  }

  /** Unpack from Context and apply */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox) {
    apply(c.expansion(),
          c.source_begin(sbox), c.source_end(sbox), c.charge_begin(sbox),
          sbox.center(),
          c.multipole(sbox));
  }
};

/** Public S2M dispatcher */
struct S2M {
  /** Forward to S2M_Helper::apply */
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::source_type& source,
                           const typename Expansion::charge_type& charge,
                           const typename Expansion::point_type& center,
                           typename Expansion::multipole_type& M) {
    typedef S2M_Helper<ExpansionTraits<Expansion>::has_S2M> S2M_H;
    S2M_H::apply(K, source, charge, center, M);
  }

  /** Forward to S2M_Helper::eval */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "S2M:"
              << "\n  " << sbox << std::endl;
#endif
    FMMTL_LOG("S2M");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef S2M_Helper<expansion_traits::has_S2M> S2M_H;
    S2M_H::eval(c, sbox);
  }
};
