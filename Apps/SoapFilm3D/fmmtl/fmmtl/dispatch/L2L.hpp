#pragma once
/** @file L2L.hpp
 * @brief Dispatch methods for the L2L stage
 *
 */

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

/** Default behavior gives a warning -- using non-existent method */
template <bool has_l2l>
struct L2L_Helper {
  template <typename... Args>
  inline static void apply(Args&&...) {
    std::cerr << "WARNING: Expansion does not have a correct L2L!\n";
  }
  template <typename... Args>
  inline static void eval(Args&&...) {
    apply();
  }
};

/** Expansion has an L2L method to dispatch to */
template <>
struct L2L_Helper<true> {
  /** The Expansion provides an L2L accumulator */
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::local_type& source,
                           typename Expansion::local_type& target,
                           const typename Expansion::point_type& translation) {
    K.L2L(source, target, translation);
  }

  /** Unpack from Context and apply */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& sbox,
                          const typename Context::target_box_type& tbox) {
    apply(c.expansion(),
          c.local(sbox),
          c.local(tbox),
          tbox.center() - sbox.center());
  }
};

// Public L2L dispatcher
struct L2L {
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::local_type& source,
                           typename Expansion::local_type& target,
                           const typename Expansion::point_type& translation) {
    typedef L2L_Helper<ExpansionTraits<Expansion>::has_L2L> L2L_H;
    L2L_H::apply(K, source, target, translation);
  }

  /** Unwrap data from Context and dispatch to the L2L evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& sbox,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "L2L:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("L2L");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef L2L_Helper<expansion_traits::has_L2L> L2L_H;
    L2L_H::eval(c, sbox, tbox);
  }
};
