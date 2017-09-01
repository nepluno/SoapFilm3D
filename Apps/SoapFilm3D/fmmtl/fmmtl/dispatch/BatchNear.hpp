#pragma once

#include <cmath>
#include <vector>

#include "fmmtl/dispatch/S2T/S2T_Compressed.hpp"

/** A lazy S2T evaluator which saves a list of pairs of boxes
 * That are sent to the S2T dispatcher on demand.
 */
template <typename Context>
class BatchNear {
  //! Kernel type
  typedef typename Context::kernel_type kernel_type;
  //! Kernel value type
  typedef typename Context::kernel_value_type kernel_value_type;

  //! Type of box
  typedef typename Context::source_box_type source_box_type;
  typedef typename Context::target_box_type target_box_type;

  std::vector<target_box_type> target_box_list;
  std::vector<std::vector<source_box_type>> source_boxes;
  unsigned box_pair_count;

  //! Box list for S2T interactions    TODO: could further compress these...
  //typedef std::pair<source_box_type, target_box_type> box_pair;
  //std::vector<box_pair> p2p_list;

  // For now, only use for GPU...
  S2T_Compressed<kernel_type>* p2p_compressed;

 public:
  BatchNear()
      : box_pair_count(0), p2p_compressed(nullptr) {
  }
  ~BatchNear() {
    delete p2p_compressed;
  }

  /** Insert a source-target box interaction to the interaction list */
  void insert(const source_box_type& s, const target_box_type& t) {
    if (source_boxes.size() <= t.index()) {
      source_boxes.resize(t.index() + 1);
      target_box_list.push_back(t);
    } else if (source_boxes[t.index()].empty()) {
      target_box_list.push_back(t);
    }

    source_boxes[t.index()].push_back(s);
    ++box_pair_count;

    //p2p_list.push_back(std::make_pair(s,t));
  }

  /** Compute all interations in the interaction list */
  void execute(Context& c) {
    FMMTL_LOG("S2T Batch");
#if defined(FMMTL_WITH_CUDA)        // XXX: Dispatch this
    if (p2p_compressed == nullptr)
      p2p_compressed =
          S2T_Compressed<kernel_type>::make(c, target_box_list, source_boxes, box_pair_count);
    p2p_compressed->execute(c);
#else
    auto t_end = target_box_list.end();
#pragma omp parallel for
    for (auto ti = target_box_list.begin(); ti < t_end; ++ti) {
      target_box_type& tb = *ti;
      auto s_end = source_boxes[tb.index()].end();
      for (auto si = source_boxes[tb.index()].begin(); si != s_end; ++si) {
        S2T::eval(c, *si, tb, S2T::ONE_SIDED());
      }
    }
#endif
  }

  /*
  class S2T_Matrix
      : public EvaluatorBase<Context> {
    ublas::compressed_matrix<kernel_value_type> A;

   public:
    virtual void execute(Context& c) {
      // printf("EvalLocalSparse::execute(Context&)\n");

      typedef typename Context::charge_type charge_type;
      ublas::vector<charge_type> charges(c.source_tree().bodies());
      std::copy(c.charge_begin(), c.charge_end(), charges.begin());

      // Call the matvec
      typedef typename Context::result_type result_type;
      ublas::vector<result_type> results = ublas::prod(A, charges);

      // Accumulate results
      std::transform(results.begin(), results.end(),
                     c.result_begin(), c.result_begin(),
                     std::plus<result_type>());
    }
  };
  */


  /** Convert the interaction list to an interaction matrix
   * by evaluating all of the elements.
   */
  /*
  ublas::compressed_matrix<kernel_value_type> to_matrix(Context& bc) {
    auto first_source = bc.source_begin();
    auto first_target = bc.target_begin();

    // Interaction list for each target body
    unsigned rows = bc.target_tree().bodies();
    unsigned cols = 0;
    unsigned nnz = 0;
    std::vector<std::vector<unsigned>> csr(rows);

    for (const box_pair& b2b : p2p_list) {
      const source_box_type& box1 = b2b.first;
      const target_box_type& box2 = b2b.second;

      auto source_end = bc.source_end(box1);
      auto target_end = bc.target_end(box2);
      for (auto t = bc.target_begin(box2); t != target_end; ++t) {
        // Row
        unsigned i = t - first_target;
        std::vector<unsigned>& csri = csr[i];

        for (auto s = bc.source_begin(box1); s != source_end; ++s) {
          // Column
          unsigned j = s - first_source;

          //FMMTL_ASSERT(std::find(csri.begin(), csri.end(), j) == csri.end());
          csri.push_back(j);
          ++nnz;
          cols = std::max(cols, j);
        }
      }
    }
    ++cols;

    // The precomputed interaction matrix
    ublas::compressed_matrix<kernel_value_type> m(rows, cols, nnz);

    typedef typename kernel_type::source_type source_type;
    typedef typename kernel_type::target_type target_type;
    for (unsigned i = 0; i < csr.size(); ++i) {
      // Insert elements to compressed_matrix in-order for efficiency
      std::sort(csr[i].begin(), csr[i].end());
      const target_type& target = first_target[i];

      for (unsigned j : csr[i]) {
        const source_type& source = first_source[j];
        m.push_back(i, j, bc.kernel()(target, source));
      }
    }

    return m;
  }
  */
};
