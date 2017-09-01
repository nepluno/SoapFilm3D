#pragma once

#include "fmmtl/dispatch/Dispatchers.hpp"


/** @brief Process the boxes from bottom to top
 * concept Tree {
 *   unsigned levels();                       // Levels in the tree, root:0
 *   box_iterator box_begin(unsigned level);  // Iterator range to boxes of level
 *   box_iterator box_end(unsigned level);
 * }
 * concept Evaluator {
 *   void operator()(typename Tree::box_type& b);   // Process box b
 * }
 */
struct UpwardPass {
  template <typename Tree, class Evaluator>
	inline static void eval(Tree& tree, Evaluator& eval) {
		// For the lowest level up to the highest level
		for (int l = tree.levels()-1; l >= 0; --l) {
			// For all boxes at this level
			auto b_end = tree.box_end(l);
      for (auto bit = tree.box_begin(l); bit != b_end; ++bit) {
        auto box = *bit;
        eval(box);
			}
		}
	}
};

/** Helper for computing the multipole for a box and all sub-boxes
 */
struct ComputeM {
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox) {
    INITM::eval(c, sbox);

    if (sbox.is_leaf()) {
      // Compute the multipole from the box's sources
      S2M::eval(c, sbox);
    } else {
      auto c_end = sbox.child_end();
      for (auto cit = sbox.child_begin(); cit != c_end; ++cit) {
        typename Context::source_box_type child = *cit;

        // Recursively initialize the multipole
        ComputeM::eval(c, child);
        // Accumulate the child's multipole into sbox
        M2M::eval(c, child, sbox);
      }
    }
  }
};
