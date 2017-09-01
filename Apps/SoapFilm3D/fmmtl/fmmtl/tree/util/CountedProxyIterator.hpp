#pragma once

#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/counting_iterator.hpp>

namespace fmmtl {

/** CountedProxyIterator
 * @brief A random-access iterator for indexed proxy objects.
 * A counted iterator over a type construted with an index-pointer pair.
 * @tparam T       The value object that supports construction: T(I, Friend*)
 * @tparam Friend  The friend class and pointer type
 * @tparam I       The index type
 *
 * Note: Since this dereferences to a value rather than a reference,
 *       it does not fully satisfy the random-access-iterator concept. Thus,
 *       this should not be implemented with boost::transform_iterator.
 */
template <typename T, typename Friend, typename I = std::size_t>
struct CountedProxyIterator
    : boost::iterator_adaptor<CountedProxyIterator<T,Friend,I>,  // Derived
                              boost::counting_iterator<I>,       // BaseType
                              T,                                 // Value
                              std::random_access_iterator_tag,   // Category
                              T>                                 // Reference
{
  typedef I size_type;
  //! Construct an invalid iterator
  CountedProxyIterator() {}
  //! The index of this iterator
  size_type index() const {
    return *(this->base());
  }
 private:
  Friend* p_;
  friend Friend;
  CountedProxyIterator(I idx, Friend* p)
      : CountedProxyIterator::iterator_adaptor(boost::counting_iterator<I>(idx)),
        p_(p) {
  }
  friend class boost::iterator_core_access;
  T dereference() const {
    return T(index(), p_);
  }
};

} // end namespace fmmtl
