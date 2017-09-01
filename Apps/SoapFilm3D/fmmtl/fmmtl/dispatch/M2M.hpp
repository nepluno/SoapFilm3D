#pragma once
/** @file M2M.hpp
 * @brief Dispatch methods for the M2M stage
 *
 */

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

/** Default behavior gives a warning -- using non-existent method */
template <bool has_m2m>
struct M2M_Helper {
  template <typename... Args>
  inline static void apply(Args&&...) {
    std::cerr << "WARNING: Expansion does not have a correct M2M!\n";
  }
  template <typename... Args>
  inline static void eval(Args&&...) {
    apply();
  }
};

/** Expansion has an M2M method to dispatch to */
template <>
struct M2M_Helper<true> {
  /** The Expansion provides an M2M accumulator */
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::multipole_type& source,
                           typename Expansion::multipole_type& target,
                           const typename Expansion::point_type& translation) {
    K.M2M(source, target, translation);
  }

  /** Unpack from Context and apply */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::source_box_type& tbox) {
    apply(c.expansion(),
          c.multipole(sbox),
          c.multipole(tbox),
          tbox.center() - sbox.center());
  }
};

// Public M2M dispatcher
struct M2M {
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::multipole_type& source,
                           typename Expansion::multipole_type& target,
                           const typename Expansion::point_type& translation) {
    typedef M2M_Helper<ExpansionTraits<Expansion>::has_M2M> M2M_H;
    M2M_H::apply(K, source, target, translation);
  }

  /** Unwrap data from Context and dispatch to the M2M evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::source_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "M2M:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("M2M");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef M2M_Helper<expansion_traits::has_M2M> M2M_H;
    M2M_H::eval(c, sbox, tbox);
  }
};
