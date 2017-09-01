#pragma once
/** @file INITM.hpp
 * @brief Dispatch methods for the initializing a multipole expansion
 */

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

/** Default behavior uses the default constructor */
template <bool has_initm>
struct INITM_Helper {
  template <typename multipole_type>
  inline static void apply(multipole_type& M) {
    M = multipole_type();
  }

  template <typename Expansion>
  inline static void apply(const Expansion&,
                           typename Expansion::multipole_type& M,
                           const typename Expansion::point_type&,
                           unsigned) {
    apply(M);
  }

  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& box) {
    apply(c.multipole(box));
  }
};

/** Expansion has an INITM method to dispatch to */
template <>
struct INITM_Helper<true> {
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           typename Expansion::multipole_type& M,
                           const typename Expansion::point_type& extents,
                           unsigned level) {
    K.init_multipole(M, extents, level);
  }

  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& box) {
    apply(c.expansion(),
          c.multipole(box),
          box.extents(),
          box.level());
  }
};


struct INITM {
  /** Forward to INITM_Helper::apply */
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           typename Expansion::multipole_type& M,
                           const typename Expansion::point_type& extents,
                           unsigned level) {
    typedef INITM_Helper<ExpansionTraits<Expansion>::has_init_multipole> INITM_H;
    INITM_H::apply(K, M, extents, level);
  }

  /** Forward to INITM_Helper::eval */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& box)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "INITM:"
              << "\n  " << box << std::endl;
#endif
    FMMTL_LOG("INITM");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef INITM_Helper<expansion_traits::has_init_multipole> INITM_H;
    INITM_H::eval(c, box);
  }
};
