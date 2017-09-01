#pragma once
/** @file S2L.hpp
 * @brief Dispatch methods for the S2L stage
 *
 */

#include <type_traits>

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

/** Default behavior gives a warning -- using non-existent method */
template <bool has_s2l>
struct S2L_Helper {
  template <typename... Args>
  inline static void apply(Args&&...) {
    std::cerr << "WARNING: Expansion does not have a correct S2L!\n";
  }
  template <typename... Args>
  inline static void eval(Args&&...) {
    apply();
  }
};

/** Expansion has an S2L method (scalar/vector) to dispatch to */
template <>
struct S2L_Helper<true> {
  /** The Expansion provides a vector S2L accumulator. */
  template <typename Expansion, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_vector_S2L>::type
  apply(const Expansion& K,
        SourceIter s_begin, SourceIter s_end, ChargeIter c_begin,
        const typename Expansion::point_type& center,
        typename Expansion::local_type& L) {
    K.S2L(s_begin, s_end, c_begin, center, L);
  }

  /** The Expansion provides a scalar S2L accumulator. */
  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_S2L &
                          !ExpansionTraits<Expansion>::has_vector_S2L>::type
  apply(const Expansion& K,
        const typename Expansion::source_type& source,
        const typename Expansion::charge_type& charge,
        const typename Expansion::point_type& center,
        typename Expansion::local_type& L) {
    K.S2L(source, charge, center, L);
  }

  /** The Expansion provides a scalar S2L accumulator. */
  template <typename Expansion, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_S2L &
                          !ExpansionTraits<Expansion>::has_vector_S2L>::type
  apply(const Expansion& K,
        SourceIter s_begin, SourceIter s_end, ChargeIter c_begin,
        const typename Expansion::point_type& center,
        typename Expansion::local_type& L) {
    for ( ; s_begin != s_end; ++s_begin, ++c_begin)
      apply(K, *s_begin, *c_begin, center, L);
  }

  /** Unpack from Context and apply */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox) {
    apply(c.expansion(),
          c.source_begin(sbox), c.source_end(sbox), c.charge_begin(sbox),
          tbox.center(),
          c.local(tbox));
  }
};

/** Public S2L dispatcher */
struct S2L {
  /** Forward to S2L_Helper::apply */
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::source_type& source,
                           const typename Expansion::charge_type& charge,
                           const typename Expansion::point_type& center,
                           typename Expansion::local_type& L) {
    typedef S2L_Helper<ExpansionTraits<Expansion>::has_S2L> S2L_H;
    S2L_H::apply(K, source, charge, center, L);
  }

  /** Forward to S2L_Helper::eval */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "S2L:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("S2L");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef S2L_Helper<expansion_traits::has_S2L> S2L_H;
    S2L_H::eval(c, sbox, tbox);
  }
};
