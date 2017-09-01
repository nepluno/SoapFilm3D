#pragma once
/** @file INITL.hpp
 * @brief Dispatch methods for the initializing a local expansion
 */

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

/** Default behavior uses the default constructor */
template <bool has_initl>
struct INITL_Helper {
  template <typename local_type>
  inline static void apply(local_type& L) {
    L = local_type();
  }

  template <typename Expansion>
  inline static void apply(const Expansion&,
                           typename Expansion::local_type& L,
                           const typename Expansion::point_type&,
                           unsigned) {
    apply(L);
  }

  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& box) {
    apply(c.local(box));
  }
};

/** Expansion has an INITL method to dispatch to */
template <>
struct INITL_Helper<true> {
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           typename Expansion::local_type& L,
                           const typename Expansion::point_type& extents,
                           unsigned level) {
    K.init_local(L, extents, level);
  }

  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& box) {
    apply(c.expansion(),
          c.local(box),
          box.extents(),
          box.level());
  }
};


struct INITL {
  /** Forward to INITL_Helper::apply */
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           typename Expansion::local_type& L,
                           const typename Expansion::point_type& extents,
                           unsigned level) {
    typedef INITL_Helper<ExpansionTraits<Expansion>::has_init_local> INITL_H;
    INITL_H::apply(K, L, extents, level);
  }

  /** Forward to INITL_Helper::eval */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& box)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "INITL:"
              << "\n  " << box << std::endl;
#endif
    FMMTL_LOG("INITL");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef INITL_Helper<expansion_traits::has_init_local> INITL_H;
    INITL_H::eval(c, box);
  }
};
