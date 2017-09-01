#pragma once
/** @file SingleTraversal
 * @brief Generic single-tree traversals for finding boxes
 * that satisfy some criteria.
 */

#include "fmmtl/tree/TreeRange.hpp"

namespace fmmtl {

struct DefaultVisit {
  template <typename Box>
  ChildRange<Box> operator()(const Box& b) const {
    return {b};
  }
};


/** @brief Traverse a tree
 *
 * // Whether the box has any children and can be traversed
 * concept Box {
 *   bool is_leaf() const;
 * }
 *
 * // Whether the box should be considered for traversal
 * const Prune {
 *   bool operator()(Box);
 * }
 *
 * // No longer able to recurse, evaluate box
 * const Base {
 *   void operator()(Box);
 * }
 *
 * // Iterable range of child boxes
 * // Defaults to [box.child_begin(), box.child_end())
 * const Visit {
 *   BoxRange operator()(Box);  // Iterable range of child boxes
 * }
 */
template <typename Box,
          typename Prune,
          typename Base,
          typename Visit = DefaultVisit>
void traverse(const Box& b,
              const Prune& prune,
              const Base& base_case,
              const Visit& visit_order = Visit()) {
  if (prune(b))
    return;

  if (b.is_leaf()) {
    base_case(b);
  } else {
    for (const Box& child : visit_order(b))
      traverse(child, prune, base_case, visit_order);
  }
}

}
