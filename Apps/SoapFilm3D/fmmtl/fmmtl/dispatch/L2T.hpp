#pragma once
/** @file L2T.hpp
 * @brief Dispatch methods for the L2T stage
 *
 */

#include <type_traits>

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

/** Default behavior gives a warning -- using non-existent method */
template <bool has_l2t>
struct L2T_Helper {
  template <typename... Args>
  inline static void apply(Args&&...) {
    std::cerr << "WARNING: Expansion does not have a correct L2T!\n";
  }
  template <typename... Args>
  inline static void eval(Args&&...) {
    apply();
  }
};

/** Expansion has an L2T method (scalar/vector) to dispatch to */
template <>
struct L2T_Helper<true> {
  /** The Expansion provides a vector L2T accumulator. */
  template <typename Expansion, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_vector_L2T>::type
  apply(const Expansion& K,
        const typename Expansion::local_type& L,
        const typename Expansion::point_type& center,
        TargetIter t_begin, TargetIter t_end, ResultIter r_begin) {
    K.L2T(L, center, t_begin, t_end, r_begin);
  }

  /** The Expansion provides a scalar L2T accumulator. */
  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_L2T &
                          !ExpansionTraits<Expansion>::has_vector_L2T>::type
  apply(const Expansion& K,
        const typename Expansion::local_type& L,
        const typename Expansion::point_type& center,
        const typename Expansion::target_type& target,
              typename Expansion::result_type& result) {
    K.L2T(L, center, target, result);
  }

  /** The Expansion provides a scalar L2T accumulator. */
  template <typename Expansion, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_L2T &
                          !ExpansionTraits<Expansion>::has_vector_L2T>::type
  apply(const Expansion& K,
        const typename Expansion::local_type& L,
        const typename Expansion::point_type& center,
        TargetIter t_begin, TargetIter t_end, ResultIter r_begin) {
    for ( ; t_begin != t_end; ++t_begin, ++r_begin)
      apply(K, L, center, *t_begin, *r_begin);
  }

  /** Unpack from Context and apply */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& tbox) {
    apply(c.expansion(),
          c.local(tbox),
          tbox.center(),
          c.target_begin(tbox), c.target_end(tbox),
          c.result_begin(tbox));
  }
};

// Public L2T dispatcher
struct L2T {
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::local_type& L,
                           const typename Expansion::point_type& center,
                           const typename Expansion::target_type& target,
                           typename Expansion::result_type& result) {
    typedef L2T_Helper<ExpansionTraits<Expansion>::has_L2T> L2T_H;
    L2T_H::apply(K, L, center, target, result);
  }

  /** Forward to L2T_Helper::eval */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "L2T:"
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("L2T");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef L2T_Helper<expansion_traits::has_L2T> L2T_H;
    L2T_H::eval(c, tbox);
  }
};
