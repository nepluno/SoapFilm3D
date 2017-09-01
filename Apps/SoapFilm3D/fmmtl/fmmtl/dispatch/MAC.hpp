#pragma once
/** @file MAC.hpp
 * @brief Dispatch methods for the Multipole Acceptance Criteria
 */

#include <type_traits>

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

class MAC {
  /** If no other MAC dispatcher matches */
  template <typename Context>
  inline static
  typename std::enable_if<!ExpansionTraits<typename Context::expansion_type>::has_dynamic_MAC, bool>::type
  eval_mac(const Context& c,
           const typename Context::source_box_type& sbox,
           const typename Context::target_box_type& tbox) {
    return c.mac(sbox,tbox);
  }

  template <typename Context>
  inline static
  typename std::enable_if<ExpansionTraits<typename Context::expansion_type>::has_dynamic_MAC, bool>::type
  eval_mac(const Context& c,
           const typename Context::source_box_type& sbox,
           const typename Context::target_box_type& tbox) {
    return c.mac(sbox,tbox) && c.expansion().MAC(c.multipole(sbox), c.local(tbox));
  }

 public:

  /** Unwrap data from Context and dispatch to the MAC evaluator
   */
  template <typename Context>
  inline static bool eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "MAC:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("MAC");

    return MAC::eval_mac(c, sbox, tbox);
  }
};
