#pragma once
/** @file BoundingBox.hpp
 * @brief Define the BoundingBox class for ND bounding boxes. */

#include <iostream>
#include <algorithm>
#include <cmath>

#include "fmmtl/meta/dimension.hpp"

namespace fmmtl {

/** @class BoundingBox
 * @brief Class representing ND bounding boxes that always contain at least
 *        one point.
 *
 * A BoundingBox is a ND volume.
 * Its fundamental operations are
 *   contains(), which tests whether a point is in the volume, and
 *   operator|=(), which extends the volume to ensure it contains a point.
 *
 * BoundingBoxes are implemented as ND rectangular cuboids whose
 * sides are aligned with the principal axes.
 *
 * The point_type of a BoundingBox satisfies the concept
 * concept point_type {
 *   point_type(const point_type&);            // Copy constructible
 *   point_type& operator=(const point_type&); // Assignable
 *   [constexpr] unsigned size();              // Dimension of the point
 *   value_type& operator[](unsigned i);       // mutable access to coordinate i
 * };
 * TODO: require only begin()/end()?
 * and value_type is comparable and assignable.
 */
template <typename POINT>
class BoundingBox {
 public:
  typedef POINT point_type;

  /** Construct the minimal bounding box containing @a p.
   * @post contains(@a p) && min() == @a p && max() == @a p */
  explicit BoundingBox(const point_type& p = point_type())
      : min_(p), max_(p) {
  }
  /** Construct the minimal bounding box containing @a p1 and @a p2.
   * @post contains(@a p1) && contains(@a p2) */
  BoundingBox(const point_type& p1, const point_type& p2)
      : min_(p1), max_(p1) {
    *this |= p2;
  }
  /** Construct a bounding box containing the points in [first, last).
   * @pre first != last
   */
  template <typename IT>
  BoundingBox(IT first, IT last)
      : min_(*first), max_(min_) {
    FMMTL_ASSERT(first != last);
    insert(++first, last);
  }

  /** Test if the bounding box is empty (contains exactly one point). */
  bool empty() const {
    for (unsigned i = 0; i != min_.size(); ++i)
      if (min_[i] != max_[i])
        return false;
    return true;
  }

  /** Return the minimum corner of the bounding box.
   * @post empty() || contains(min()) */
  const point_type& min() const {
    return min_;
  }

  /** Return the maximum corner of the bounding box.
   * @post contains(max()) */
  const point_type& max() const {
    return max_;
  }

  /** Test if point @a p is in the bounding box. */
  bool contains(const point_type& p) const {
    for (unsigned i = 0; i != min_.size(); ++i)
      if (p[i] < min_[i] || p[i] > max_[i])
        return false;
    return true;
  }

  /** Test if @a b is entirely within this bounding box.
   * @returns true if all @a p with @a b.contains(@a p) implies contains(@a p) */
  bool contains(const BoundingBox& b) const {
    for (unsigned i = 0; i != min_.size(); ++i)
      if (b.min_[i] < min_[i] || b.min_[i] > max_[i] ||
          b.max_[i] < min_[i] || b.max_[i] > max_[i])
        return false;
    return true;
  }

  /** Test if @a b intersects this bounding box.
   * @returns true if there exists @a p such that
   *            contains(@a p) && b.contains(@a p) */
  bool intersects(const BoundingBox& b) const {
    for (unsigned i = 0; i != min_.size(); ++i)
      if (b.min_[i] > max_[i] || b.max_[i] < min_[i])
        return false;
    return true;
  }

  /** Extend the bounding box to contain @a p.
   * @post contains(@a p) is true
   * @post For all @a x with old contains(@a x),
             then new contains(@a x) is true. */
  BoundingBox& operator|=(const point_type& p) {
    for (unsigned i = 0; i != min_.size(); ++i) {
      if (p[i] < min_[i])  min_[i] = p[i];
      if (p[i] > max_[i])  max_[i] = p[i];
    }
    return *this;
  }

  /** Extend the bounding box to contain @a b.
   * @post contains(@a b)
   * @post For all @a x with old contains(@a x) or @a b.contains(@a x),
   *         then new contains(@a x) is true. */
  BoundingBox& operator|=(const BoundingBox& b) {
    return (*this |= b.min()) |= b.max();
  }

  /** Extend the bounding box to contain the points in [first, last).
   * @post For all @a p in [@a first, @a last), contains(@a p) is true.
   * @post For all @a x with old contains(@a x),
   *         then new contains(@a x) is true. */
  template <typename IT>
  BoundingBox& insert(IT first, IT last) {
    for ( ; first != last; ++first)
      *this |= *first;
    return *this;
  }

 private:
  // RI: max_[i] >= min_[i] for all 0 <= i < min_.size()
  // RI: max_.size() == min_.size()
  point_type min_;
  point_type max_;
};


/** Return a bounding box that contains @a b and @a p. */
template <typename P>
BoundingBox<P> operator|(BoundingBox<P> b, const P& p) {
  return b |= p;
}
/** Return the union of @a b1 and @a b2. */
template <typename P>
BoundingBox<P> operator|(BoundingBox<P> b1, const BoundingBox<P>& b2) {
  return b1 |= b2;
}

/** Compute minimum box-to-point squared distance
 */
template <typename P>
double norm_2_sq(const BoundingBox<P>& bb, const P& p) {
  double result = 0;    // XXX: double
  for (unsigned i = 0; i < p.size(); ++i) {
    auto low = bb.min()[i] - p[i];
    low = low < 0 ? 0 : low;
    auto hi  = p[i] - bb.max()[i];
    hi = hi < 0 ? 0 : hi;
    result += (low + hi) * (low + hi);
  }
  return result;
}


/** Write a BoundingBox to an output stream.
 *
 * An empty bounding box is written as "[]".
 * A nonempty bounding box is written as "[min:max] (dim)".
 */
template <typename P>
inline std::ostream& operator<<(std::ostream& s, const BoundingBox<P>& b) {
  const unsigned dim = b.min().size();
  s << '[' << b.min()[0];
  for (unsigned i = 1; i != dim; ++i)  s << ", " << b.min()[i];
  s << " : " << b.max()[0];
  for (unsigned i = 1; i != dim; ++i)  s << ", " << b.max()[i];
  s << "] (" << b.max()[0] - b.min()[0];
  for (unsigned i = 1; i != dim; ++i)  s << ", " << b.max()[i] - b.min()[i];
  return s << ')';
}

} // end namespace fmmtl
