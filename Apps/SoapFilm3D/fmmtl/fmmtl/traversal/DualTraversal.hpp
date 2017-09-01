#pragma once
/** @file DualTraversal
 * @brief Generic dual-tree traversals for finding pairs of boxes
 * that satisfy some criteria.
 */

// TODO: Deprecate in favor of Traversal.hpp?

#include <utility>
#include <tuple>
#include <queue>
#include <stack>

namespace fmmtl {

struct depth_first {};
struct breadth_first {};

template <class Traversal, class T>
struct traversal_impl;

template <class T>
struct traversal_impl<depth_first,T> {
  std::stack<T> q;
  bool empty() const { return q.empty(); }
  T& next() { return q.top(); }
  const T& next() const { return q.top(); }
  void pop() { return q.pop(); }
  void push(const T& v) { return q.push(v); }
  template <class... Args>
  void emplace(Args&&... args) {
    return q.emplace(std::forward<Args>(args)...);
  }
};

template <class T>
struct traversal_impl<breadth_first,T> {
  std::queue<T> q;
  bool empty() const { return q.empty(); }
  T& next() { return q.front(); }
  const T& next() const { return q.front(); }
  void pop() { return q.pop(); }
  void push(const T& v) { return q.push(v); }
  template <class... Args>
  void emplace(Args&&... args) {
    return q.emplace(std::forward<Args>(args)...);
  }
};

/** Generic dual-tree traversal with near- and far-field operators.
 *
 * Recursively split the boxes (by volume) until either
 *   1) @a far_eval returns true
 *   2) both boxes are leaves
 * If @a far_eval returns false, then the largest box is split and recursed on.
 * If the boxes cannot be split (both are leaves and have no children), then
 * @near_field is called on the pair.
 *
 * concept SourceBox/TargetBox {
 *   box_iterator child_begin() const;
 *   box_iterator child_end() const;
 *   bool is_leaf() const;   // XXX, abstract on base case?
 *   double radius_sq() const;  // XXX, abstract on splitting criteria?
 * }
 *
 * concept FarEvaluator {
 *   // Returns true if @a s and @a t may be considered "far"
 *   // and should not be recursed on.
 *   bool operator()(SourceBox s, TargetBox t);
 * }
 * concept NearEvaluator {
 *   // Operator to apply to two boxes that have not been proven to be "far"
 *   // but cannot be recursed on.
 *   void operator(SourceBox, TargetBox);
 * }
 */
template <typename TraversalOrder = breadth_first,
          class SourceBox, class TargetBox,
          class NearEvaluator, class FarEvaluator>
inline void traverse_nearfar(SourceBox sbox, TargetBox tbox,
                             NearEvaluator& near_eval, FarEvaluator& far_eval) {
  // Queue based traversal
  typedef std::pair<SourceBox, TargetBox> BoxPair;
  traversal_impl<TraversalOrder, BoxPair> pairQ;

  // Initialize
  if (!far_eval(sbox, tbox))
    pairQ.emplace(sbox, tbox);

  // Loop until empty
  while (!pairQ.empty()) {
    std::tie(sbox, tbox) = pairQ.next();
    pairQ.pop();

    const char code = (sbox.is_leaf() << 1) | (tbox.is_leaf() << 0);
    switch (code) {
      case 0: {             // sbox and tbox are not leaves
        // Split the larger of the two into children and interact
        if (sbox.radius_sq() > tbox.radius_sq()) {
          case 1:           // tbox is a leaf, sbox is not a leaf
            // Split the source box into children
            auto c_end = sbox.child_end();
            for (auto cit = sbox.child_begin(); cit != c_end; ++cit) {
              SourceBox cbox = *cit;
              if (!far_eval(cbox, tbox))
                pairQ.emplace(cbox, tbox);
            }
        } else {
          case 2:           // sbox is a leaf, tbox is not a leaf
            // Split the target box into children
            auto c_end = tbox.child_end();
            for (auto cit = tbox.child_begin(); cit != c_end; ++cit) {
              TargetBox cbox = *cit;
              if (!far_eval(sbox, cbox))
                pairQ.emplace(sbox, cbox);
            }
        }
      } continue;

      case 3: {             // sbox and tbox are leaves
        near_eval(sbox, tbox);
      } continue;
    } // end switch
  } // end while
}


/** Alternative implementation dispatching all logic to the evaluator.
 *
 * If the evaluator returns
 *  0: Base case, neither box is split
 *  1: Split the source box and recurse on all pairs
 *  2: Split the target box and recurse on all pairs
 *  3: Split both boxes and recurse on all pairs
 *
 * concept SourceBox/TargetBox {
 *   box_iterator child_begin() const;
 *   box_iterator child_end() const;
 * }
 * concept Evaluator {
 *   // Returns 0, 1, 2, or 3 as defined above.
 *   int operator()(SourceBox s, TargetBox t);
 * }
 */
template <class TraversalOrder = breadth_first,
          class SourceBox, class TargetBox,
          class Evaluator>
inline void traverse_if(SourceBox sbox, TargetBox tbox, Evaluator& eval) {

  // Queue based traversal
  typedef std::pair<SourceBox, TargetBox> BoxPair;
  traversal_impl<TraversalOrder, BoxPair> pairQ;

  // Initialize
  pairQ.emplace(sbox, tbox);

  while (!pairQ.empty()) {
    std::tie(sbox, tbox) = pairQ.next();
    pairQ.pop();

    switch(eval(sbox, tbox)) {
      case 0: {
      } continue;
      case 1: {
        // Split the source box into children
        auto c_end = sbox.child_end();
        for (auto cit = sbox.child_begin(); cit != c_end; ++cit)
          pairQ.emplace(*cit, tbox);       // traverse_if(*cit, tbox, eval)?
      } continue;
      case 2: {
        // Split the target box into children
        auto c_end = tbox.child_end();
        for (auto cit = tbox.child_begin(); cit != c_end; ++cit)
          pairQ.emplace(sbox, *cit);       // traverse_if(sbox, *cit, eval)?
      } continue;
      case 3: {
        // Split both into children
        auto cs_end = sbox.child_end();
        auto ct_end = tbox.child_end();
        auto ct_beg = tbox.child_begin();
        for (auto csit = sbox.child_begin(); csit != cs_end; ++csit) {
          for (auto ctit = ct_beg; ctit != ct_end; ++ctit) {
            pairQ.emplace(*csit, *ctit);   // traverse_if(*csit, *ctit, eval)?
          }
        }
      } continue;
    } // end switch
  } // end while
}


} // end namespace fmmtl
