#pragma once

#include <vector>

#include "M2T.hpp"
#include "M2L.hpp"
#include "S2L.hpp"

#include "fmmtl/meta/kernel_traits.hpp"

/** A lazy far-field evaluator which saves a list of pairs of boxes
 * That are sent to the M2T/M2L/S2L dispatcher on demand.
 */
template <typename Context>
class BatchFar {
  //! Type of box
  typedef typename Context::source_box_type source_box_type;
  typedef typename Context::target_box_type target_box_type;

  //! CSR-like storage of box pairs
  // XXX: Need to also protect against parent-box races!!
  std::vector<target_box_type> target_box_list;
  std::vector<std::vector<source_box_type>> source_boxes;

 public:

  /** Insert a source-target box interaction to the interaction list */
  void insert(const source_box_type& s, const target_box_type& t) {
    if (source_boxes.size() <= t.index()) {
      source_boxes.resize(t.index() + 1);
      target_box_list.push_back(t);
    } else if (source_boxes[t.index()].empty()) {
      target_box_list.push_back(t);
    }

    source_boxes[t.index()].push_back(s);
  }

  /** Compute all interations in the interaction list */
  void execute(Context& c) {
    FMMTL_LOG("FarBatch");
    // Choose which operators are available and dispatch to it
    // XXX: Hacky version

    if (ExpansionTraits<typename Context::expansion_type>::has_L2T) {
      // Need iterator style for OMP compatibility.
#pragma omp parallel for
      for (auto ti = target_box_list.begin(); ti < target_box_list.end(); ++ti) {
        target_box_type& tb = *ti;
        auto s_end = source_boxes[tb.index()].end();
        for (auto si = source_boxes[tb.index()].begin(); si != s_end; ++si) {
          // A hacky adaption on the operator graph
          // TODO: Actually measure/autotune these
          if (ExpansionTraits<typename Context::expansion_type>::has_M2L) {
            M2L::eval(c, *si, tb);
          } else if (ExpansionTraits<typename Context::expansion_type>::has_L2T) {
            S2L::eval(c, *si, tb);
          }
        }
      }
    } else if (ExpansionTraits<typename Context::expansion_type>::has_S2M) {
      // XXX: can't run in parallel with the current data structure
      for (auto tb : target_box_list)
        for (auto sb : source_boxes[tb.index()])
          M2T::eval(c, sb, tb);
    }
  }
};
