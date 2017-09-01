#pragma once
/** @file Expansion
 *
 * CRTP base class for Expansions.
 * This allows the future extension of Expansions with methods that
 * need not be implemented within the Expansion.
 * Additionally, the traits system could be incorporated here to give the
 * library improved access and insight to these classes.
 *
 * At the moment, these may not be fully necessary and are subject to change.
 */

#include <utility>

// Protect against C++11 when compiling with nvcc
#if !defined(FMMTL_KERNEL)
# include <type_traits>
# include <boost/range/has_range_iterator.hpp>
# include "fmmtl/KernelMatrix.hpp"
#endif

#include "fmmtl/config.hpp"
#include "fmmtl/meta/kernel_traits.hpp"


namespace fmmtl {


template <class Kernel, class DerivedExpansion>
struct Expansion
    : public Kernel {
  // The type of kernel this expansion is representing
  typedef Kernel            kernel_type;
  // The type of the derived expansion implementation
  typedef DerivedExpansion  expansion_type;

  // Default constructor
  Expansion() {}

  // Forward the construction on a Kernel to the copy constructor
  Expansion(const Kernel& K)
      : Kernel(K) {}

  /** Cast down to the derived expansion */
  expansion_type& expansion() {
    return static_cast<expansion_type&>(*this);
  }
  /** Cast down to the derived expansion */
  const expansion_type& expansion() const {
    return static_cast<const expansion_type&>(*this);
  }
  /** Cast up to the base kernel */
  kernel_type& kernel() {
    return static_cast<kernel_type&>(*this);
  }
  /** Cast up to the base kernel */
  const kernel_type& kernel() const {
    return static_cast<const kernel_type&>(*this);
  }

  /** Make sure the op() is not overridden */
  using Kernel::operator();

#if !defined(FMMTL_KERNEL)
  typedef typename KernelTraits<kernel_type>::target_type target_type;
  typedef typename KernelTraits<kernel_type>::source_type source_type;

  // Range constructor
  // Creates a kernel_matrix from a set of sources and targets
  template <class TR, class SR>
  typename std::enable_if<
    !std::is_same<TR, target_type>::value &&
    !std::is_same<SR, source_type>::value &&
    boost::has_range_iterator<TR>::value &&
    boost::has_range_iterator<SR>::value,
  kernel_matrix<expansion_type, TR, SR> >::type
  operator()(const TR& t,
             const SR& s) const {
    // TODO: Optimize the kernel matrix on symmetric source/target sets
    if (&s == &t)
      return kernel_matrix<expansion_type,TR,SR>(expansion(), t);
    else
      return kernel_matrix<expansion_type,TR,SR>(expansion(), t, s);
  }
#endif
};

} // end namespace fmmtl
