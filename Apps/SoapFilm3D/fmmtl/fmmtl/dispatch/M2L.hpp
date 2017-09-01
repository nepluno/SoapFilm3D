#pragma once
/** @file M2L.hpp
 * @brief Dispatch methods for the M2L stage
 *
 */

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

/** Default behavior gives a warning -- using non-existent method */
template <bool has_m2l>
struct M2L_Helper {
  template <typename... Args>
  inline static void apply(Args&&...) {
    std::cerr << "WARNING: Expansion does not have a correct M2L!\n";
  }
  template <typename... Args>
  inline static void eval(Args&&...) {
    apply();
  }
};

/** Expansion has an M2L method to dispatch to */
template <>
struct M2L_Helper<true> {
  /** The Expansion provides an M2L accumulator */
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::multipole_type& M,
                           typename Expansion::local_type& L,
                           const typename Expansion::point_type& translation) {
    K.M2L(M, L, translation);
  }

  /** Unpack from Context and apply */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox) {
    apply(c.expansion(),
          c.multipole(sbox),
          c.local(tbox),
          tbox.center() - sbox.center());
  }
};

// Public M2L dispatcher
struct M2L {
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::multipole_type& M,
                           typename Expansion::local_type& L,
                           const typename Expansion::point_type& translation) {
    typedef M2L_Helper<ExpansionTraits<Expansion>::has_M2L> M2L_H;
    M2L_H::apply(K, M, L, translation);
  }

  /** Forward to M2L_Helper::eval */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "M2L:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("M2L");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef M2L_Helper<expansion_traits::has_M2L> M2L_H;
    M2L_H::eval(c, sbox, tbox);
  }
};
