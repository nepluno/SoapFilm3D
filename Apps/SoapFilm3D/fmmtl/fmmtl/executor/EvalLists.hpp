#pragma once

#include "fmmtl/executor/Evaluator.hpp"

#include "fmmtl/traversal/Upward.hpp"
#include "fmmtl/traversal/DualTraversal.hpp"
#include "fmmtl/traversal/Downward.hpp"

#include "fmmtl/dispatch/Dispatchers.hpp"
#include "fmmtl/tree/TreeRange.hpp"
#include "fmmtl/meta/kernel_traits.hpp"


template <class Context>
class EvalLists
    : public EvaluatorBase<Context>
{
  typedef typename Context::source_box_type source_box;
  typedef typename Context::target_box_type target_box;

  BatchNear<Context> near_batch_;
  BatchFar<Context> far_batch_;

  struct UpDispatch {
    Context& c_;
    UpDispatch(Context& c) : c_(c) {}

    inline void operator()(const source_box& box) {
      if (ExpansionTraits<typename Context::expansion_type>::has_M2M) {
        if (box.is_leaf()) {
          // If leaf, make S2M calls
          S2M::eval(c_, box);
        } else {
          // If not leaf, then for all the children M2M
          for (auto&& cbox : children(box))
            M2M::eval(c_, cbox, box);
        }
      } else {
        // Rely solely on S2M
        S2M::eval(c_, box);
      }
    }
  };
  struct DownDispatch {
    Context& c_;
    DownDispatch(Context& c) : c_(c) {}

    inline void operator()(const target_box& box) {
      if (ExpansionTraits<typename Context::expansion_type>::has_L2L) {
        if (box.is_leaf()) {
          // If leaf, make L2T calls
          L2T::eval(c_, box);
        } else {
          // If not leaf, then for all children L2L
          for (auto&& cbox : children(box))
            L2L::eval(c_, box, cbox);
        }
      } else {
        // Rely solely on L2T
        L2T::eval(c_, box);
      }
    }
  };

 public:

  EvalLists(Context& c) {
    // Construct functors for dispatched near and far operators
    auto far_batcher = [&c,this](const source_box& s, const target_box& t) {
      if (MAC::eval(c,s,t)) {
        far_batch_.insert(s,t);
        return true;
      }
      return false;
    };
    auto near_batcher = [this](const source_box& s, const target_box& t) {
      near_batch_.insert(s,t);
    };
    // Determine the box interactions
    fmmtl::traverse_nearfar(c.source_tree().root(), c.target_tree().root(),
                            near_batcher, far_batcher);
  }

  void execute(Context& c) {
    // Launch the p2p early (potentially asynchronously?)
    near_batch_.execute(c);

    // Initialize all the multipoles and locals (not all may be needed)
    for (auto&& sbox : boxes(c.source_tree()))
      INITM::eval(c, sbox);
    for (auto&& tbox : boxes(c.target_tree()))
      INITL::eval(c, tbox);

    // Perform the upward pass (not all may be needed)
    // TODO: Use far_field batch to only initialize and compute used multipoles
    if (ExpansionTraits<typename Context::expansion_type>::has_S2M) {
      UpDispatch up(c);
      UpwardPass::eval(c.source_tree(), up);
    }

    // Perform the source-target box interactions
    far_batch_.execute(c);

    // Perform the downward pass (not all may be needed)
    // TODO: Use far_field batch to only initialize and compute used locals
    if (ExpansionTraits<typename Context::expansion_type>::has_L2T) {
      DownDispatch down(c);
      DownwardPass::eval(c.target_tree(), down);
    }
  }
};


template <class Context, class Options>
EvaluatorBase<Context>* make_eval_lists(Context& c, Options&) {
  return new EvalLists<Context>(c);
}
