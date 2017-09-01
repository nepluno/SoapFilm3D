#pragma once

#include "fmmtl/executor/Evaluator.hpp"

#include "fmmtl/traversal/Upward.hpp"
#include "fmmtl/traversal/DualTraversal.hpp"
#include "fmmtl/traversal/Downward.hpp"

#include "fmmtl/dispatch/Dispatchers.hpp"
#include "fmmtl/tree/TreeRange.hpp"
#include "fmmtl/meta/kernel_traits.hpp"


template <class Context>
class EvalTraverse
    : public EvaluatorBase<Context>
{
  typedef typename Context::source_box_type source_box;
  typedef typename Context::target_box_type target_box;

 public:

  EvalTraverse(Context&) {
  }

  void execute(Context& c) {
    // Initialize all the multipoles and locals (not all may be needed)
    auto s_end = c.source_tree().box_end();
    for (auto bi = c.source_tree().box_begin(); bi != s_end; ++bi)
      INITM::eval(c, *bi);
    auto t_end = c.target_tree().box_end();
    for (auto bi = c.target_tree().box_begin(); bi != t_end; ++bi)
      INITL::eval(c, *bi);

    // Perform the upward pass (not all may be needed)
    auto up_dispatch = [&c](const source_box& box) {
      if (box.is_leaf()) {
        // If leaf, make S2M calls
        S2M::eval(c, box);
      } else {
        // If not leaf, then for all the children M2M
        auto c_end = box.child_end();
        for (auto cit = box.child_begin(); cit != c_end; ++cit)
          M2M::eval(c, *cit, box);
      }
    };
    UpwardPass::eval(c.source_tree(), up_dispatch);

    // Perform the source-target box interactions
    auto far_dispatch = [&c](const source_box& s, const target_box& t) {
      if (MAC::eval(c,s,t)) {
        M2L::eval(c,s,t);
        return true;
      }
      return false;
    };
    auto near_dispatch = [&c](const source_box& s, const target_box& t) {
      S2T::eval(c,s,t,S2T::ONE_SIDED());
    };
    fmmtl::traverse_nearfar(c.source_tree().root(), c.target_tree().root(),
                            near_dispatch, far_dispatch);

    // Perform the downward pass (not all may be needed)
    auto down_dispatch = [&c](const target_box& box) {
      if (box.is_leaf()) {
        // If leaf, make L2T calls
        L2T::eval(c, box);
      } else {
        // If not leaf, then for all children L2L
        auto c_end = box.child_end();
        for (auto cit = box.child_begin(); cit != c_end; ++cit)
          L2L::eval(c, box, *cit);
      }
    };
    DownwardPass::eval(c.target_tree(), down_dispatch);
  }
};


template <class Context, class Options>
EvaluatorBase<Context>* make_eval_traverse(Context& c, Options&) {
  return new EvalTraverse<Context>(c);
}
