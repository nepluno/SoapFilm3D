#pragma once

#include <vector>

// Don't think this needs to be virtual...
template <typename Context>
struct EvaluatorBase {
  virtual ~EvaluatorBase() {};
  virtual void execute(Context&) = 0;
};


template <typename Context>
class EvaluatorCollection
    : public EvaluatorBase<Context> {
  typedef Context context_type;

  //! Evaluator algorithms to apply
  std::vector<EvaluatorBase<context_type>*> evals_;

 public:

  virtual ~EvaluatorCollection() {
    for (auto eval : evals_)
      delete eval;
  }

  void insert(EvaluatorBase<context_type>* eval) {
    if (eval)
      evals_.push_back(eval);
  }

  void execute(context_type& context) const {
    for (auto eval : evals_)
      eval->execute(context);
  }
};
