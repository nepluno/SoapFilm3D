#pragma once

#include "fmmtl/executor/EvalLists.hpp"
#include "fmmtl/executor/EvalTraverse.hpp"

#include "fmmtl/meta/kernel_traits.hpp"

template <typename Context, typename Options>
EvaluatorBase<Context>* make_evaluator(Context& c, Options& opts) {
  // Determine the type of Evaluator
  // Statically from the type of Options
  // Dynamically from the Options input

  // For now
  if (ExpansionTraits<typename Context::expansion_type>::has_dynamic_MAC)
    return make_eval_traverse(c, opts);
  else
    return make_eval_lists(c, opts);
}
