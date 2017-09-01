#pragma once
/** @file M2T.hpp
 * @brief Dispatch methods for the M2T stage
 *
 */

#include <type_traits>

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

/** Default behavior gives a warning -- using non-existent method */
template <bool has_m2t>
struct M2T_Helper {
  template <typename... Args>
  inline static void apply(Args&&...) {
    std::cerr << "WARNING: Expansion does not have a correct M2T!\n";
  }
  template <typename... Args>
  inline static void eval(Args&&...) {
    apply();
  }
};

/** Expansion has an M2T method (scalar/vector) to dispatch to */
template <>
struct M2T_Helper<true> {
  /** The Expansion provides a vector M2T accumulator. */
  template <typename Expansion, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_vector_M2T>::type
  apply(const Expansion& K,
        const typename Expansion::multipole_type& M,
        const typename Expansion::point_type& center,
        TargetIter t_begin, TargetIter t_end, ResultIter r_begin) {
    K.M2T(M, center, t_begin, t_end, r_begin);
  }

  /** The Expansion provides a scalar M2T accumulator. */
  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_M2T &
                          !ExpansionTraits<Expansion>::has_vector_M2T>::type
  apply(const Expansion& K,
        const typename Expansion::multipole_type& M,
        const typename Expansion::point_type& center,
        const typename Expansion::target_type& target,
              typename Expansion::result_type& result) {
    K.M2T(M, center, target, result);
  }

  /** The Expansion provides a scalar M2T accumulator. */
  template <typename Expansion, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_M2T &
                          !ExpansionTraits<Expansion>::has_vector_M2T>::type
  apply(const Expansion& K,
        const typename Expansion::multipole_type& M,
        const typename Expansion::point_type& center,
        TargetIter t_begin, TargetIter t_end, ResultIter r_begin) {
    for ( ; t_begin != t_end; ++t_begin, ++r_begin)
      apply(K, M, center, *t_begin, *r_begin);
  }

  /** Unpack from Context and apply */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox) {
    apply(c.expansion(),
          c.multipole(sbox), sbox.center(),
          c.target_begin(tbox), c.target_end(tbox), c.result_begin(tbox));
  }
};

/** Public M2T dispatcher */
struct M2T {
  /** Forward to M2T_Helper::apply */
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::multipole_type& M,
                           const typename Expansion::point_type& center,
                           const typename Expansion::target_type& target,
                           typename Expansion::result_type& result) {
    typedef M2T_Helper<ExpansionTraits<Expansion>::has_M2T> M2T_H;
    M2T_H::apply(K, M, center, target, result);
  }

  /** Forward to M2T_Helper::eval */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "M2T:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("M2T");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef M2T_Helper<expansion_traits::has_M2T> M2T_H;
    M2T_H::eval(c, sbox, tbox);
  }
};
