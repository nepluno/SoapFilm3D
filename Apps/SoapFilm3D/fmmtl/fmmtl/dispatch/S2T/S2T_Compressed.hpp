#pragma once
/** @brief Header file for the GPU_S2T class.
 *
 * Note: This header file may be compiled with nvcc and must use C++03.
 */

#include <iostream>
#include <vector>

#include "fmmtl/config.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

template <typename Kernel>
class S2T_Compressed {
 public:
  FMMTL_IMPORT_KERNEL_TRAITS(Kernel);

  // Supporting data
  void* data_;

  // Device data for S2T computation
  std::pair<unsigned,unsigned>* target_ranges_;
  unsigned* source_range_ptrs_;
  std::pair<unsigned,unsigned>* source_ranges_;

  // Device source and target arrays
  source_type* sources_;
  target_type* targets_;

  S2T_Compressed();

	S2T_Compressed(std::vector<std::pair<unsigned,unsigned> >& target_ranges,
                 std::vector<unsigned>& target_ptrs,
                 std::vector<std::pair<unsigned,unsigned> >& source_ranges,
                 const std::vector<source_type>& sources,
                 const std::vector<target_type>& targets);

  ~S2T_Compressed();

  // Convenience function
  template <class Context>
  void execute(Context& c) {
    return execute(c.kernel(),
                   c.charge_begin(), c.charge_end(),
                   c.result_begin(), c.result_end());
  }

  template <typename ChargeIter, typename ResultIter>
  void execute(const Kernel& K,
               ChargeIter cfirst, ChargeIter clast,
               ResultIter rfirst, ResultIter rlast) {
    // TODO: Ugh, iterator type hiding via copy
    std::vector<charge_type> charges(cfirst, clast);
    std::vector<result_type> results(rfirst, rlast);

    execute(K, charges, results);

    std::copy(results.begin(), results.end(), rfirst);
  }

  void execute(const Kernel& K,
               const std::vector<charge_type>& charges,
               std::vector<result_type>& results);

  static void execute(const Kernel& K,
                      const std::vector<source_type>& s,
                      const std::vector<charge_type>& c,
                      const std::vector<target_type>& t,
                      std::vector<result_type>& r);

  /** Construct a S2T_Compressed object by taking
   * associated source ranges and target ranges and constructing a compressed
   * representation.
   *
   * @param srfirst,srlast  A range of source ranges
   * @param trfirst         A range of target ranges
   *                          The source/target ranges are associated
   * @param sources         The sources that the source ranges map into
   * @param targets         The targets that the target ranges map into
   *
   * @pre Target ranges are disjoint. No two target ranges overlap.
   *
   * @note Creates a CSR-like compressed representation of the blocked matrix
   *
   * TODO: Clean up...
   */
  template <class SourceRangeIter, class TargetRangeIter>
  static
  S2T_Compressed<Kernel>*
  make(SourceRangeIter srfirst, SourceRangeIter srlast,
       TargetRangeIter trfirst,
       const std::vector<source_type>& sources,
       const std::vector<target_type>& targets) {
    unsigned num_targets = targets.size();
    //unsigned num_sources = sources.size();
    unsigned num_box_pairs = srlast - srfirst;

    // Interaction list for each target box
    // (target_first,target_last) -> {(source_first, source_last), ...}
    // TODO: faster?
    typedef std::pair<unsigned, unsigned> upair;
    std::vector<std::vector<upair> > target2sources(num_targets);
    // A list of target ranges we've seen: {(target_first, target_last), ...}
    std::vector<upair> target_ranges;

    for ( ; srfirst != srlast; ++srfirst, ++trfirst) {
      upair s_range = *srfirst;
      upair t_range = *trfirst;

      unsigned i_begin = t_range.first;
      unsigned i_end   = t_range.second;

      unsigned j_begin = s_range.first;
      unsigned j_end   = s_range.second;

      // If this is the first time we've seen this target range, record it
      if (target2sources[i_begin].empty())
        target_ranges.push_back(upair(i_begin, i_end));

      // Record this source range with this target range
      target2sources[i_begin].push_back(upair(j_begin,j_end));
    }

    unsigned num_target_ranges = target_ranges.size();

    // Construct a compressed interaction list
    std::vector<unsigned> target_ptr(num_target_ranges + 1);
    target_ptr[0] = 0;
    std::vector<upair> source_ranges(num_box_pairs);
    std::vector<upair>::iterator source_ranges_curr = source_ranges.begin();

    // For all the target ranges
    for (unsigned k = 0; k < num_target_ranges; ++k) {
      // Copy the source ranges that interact with the kth target range
      unsigned i_begin = target_ranges[k].first;
      source_ranges_curr = std::copy(target2sources[i_begin].begin(),
                                     target2sources[i_begin].end(),
                                     source_ranges_curr);

      // Record the stop index
      target_ptr[k+1] = source_ranges_curr - source_ranges.begin();
    }

    // Sanity checking
    FMMTL_ASSERT(target_ptr.back() == source_ranges.size());
    FMMTL_ASSERT(source_ranges_curr == source_ranges.end());

    return new S2T_Compressed<Kernel>(target_ranges,
                                      target_ptr,
                                      source_ranges,
                                      sources,
                                      targets);
  }


  template <class Context>
  static
  S2T_Compressed<typename Context::kernel_type>*
  make(Context& c,
       const std::vector<typename Context::target_box_type> t_boxes,
       const std::vector<std::vector<typename Context::source_box_type> > s_boxes,
       int num_box_pairs) {
    typedef typename Context::target_box_type target_box_type;
    typedef typename Context::source_box_type source_box_type;

    typedef std::pair<unsigned, unsigned> upair;

    // Interaction list for each target box
    // (target_first,target_last) -> {(source_first, source_last), ...}

    typename Context::source_iterator first_source = c.source_begin();
    typename Context::target_iterator first_target = c.target_begin();

    std::vector<upair> target_ranges;
    target_ranges.reserve(t_boxes.size());
    std::vector<unsigned> target_ptr;
    target_ptr.reserve(num_box_pairs+1);
    std::vector<upair> source_ranges;
    source_ranges.reserve(num_box_pairs);

    // Construct a compressed interaction list
    target_ptr.push_back(0);
    for (unsigned k = 0; k < t_boxes.size(); ++k) {
      // Copy the source ranges that interact with the kth target range
      const target_box_type& t = t_boxes[k];
      target_ranges.push_back(upair(std::distance(first_target, c.target_begin(t)),
                                    std::distance(first_target, c.target_end(t))));

      const std::vector<source_box_type>& s_list = s_boxes[t.index()];
      for (unsigned n = 0; n < s_list.size(); ++n) {
        const source_box_type& s = s_list[n];
        source_ranges.push_back(upair(std::distance(first_source, c.source_begin(s)),
                                      std::distance(first_source, c.source_end(s))));
      }

      // Record the stop index
      target_ptr.push_back(std::distance(source_ranges.begin(), source_ranges.end()));
    }

    // Sanity checking
    FMMTL_ASSERT(target_ptr.back() == source_ranges.size());

    // Copy the source and target ranges into contiguous vectors
    std::vector<source_type> sources(c.source_begin(), c.source_end());
    std::vector<target_type> targets(c.target_begin(), c.target_end());

    return new S2T_Compressed<Kernel>(target_ranges,
                                      target_ptr,
                                      source_ranges,
                                      sources,
                                      targets);
  }

  /*
  template <class Context, class BoxPairIter>
  static
  S2T_Compressed<typename Context::kernel_type>*
  make(Context& c, BoxPairIter first, BoxPairIter last) {
    typename Context::source_iterator first_source = c.source_begin();
    typename Context::source_iterator first_target = c.target_begin();

    unsigned num_targets = c.target_tree().bodies();
    //unsigned num_sources = c.source_tree().bodies();
    unsigned num_box_pairs = last - first;

    // Interaction list for each target box
    // (target_first,target_last) -> {(source_first, source_last), ...}
    // TODO: faster?
    typedef std::pair<unsigned, unsigned> upair;
    std::vector<std::vector<upair> > target2sources(num_targets);
    // A list of target ranges we've seen: {(target_first, target_last), ...}
    std::vector<upair> target_ranges;

    for ( ; first != last; ++first) {
      const typename BoxPairIter::value_type& bpair = *first;
      const typename Context::source_box_type& source_box = bpair.first;
      const typename Context::target_box_type& target_box = bpair.second;

      // Target boxes need to be leaf boxes because the computations are
      // grouped by disjoint target ranges
      // TODO: Generalize?
      FMMTL_ASSERT(target_box.is_leaf());

			// Target range
			unsigned i_begin = c.target_begin(target_box) - first_target;
      unsigned i_end   = c.target_end(target_box) - first_target;

			// Source range
			unsigned j_begin = c.source_begin(source_box) - first_source;
			unsigned j_end   = c.source_end(source_box) - first_source;

      // If this is the first time we've seen this target range
      if (target2sources[i_begin].empty())
        target_ranges.push_back(upair(i_begin, i_end));

      target2sources[i_begin].push_back(upair(j_begin,j_end));
    }

    unsigned num_target_ranges = target_ranges.size();

    // Construct a compressed interaction list
    std::vector<unsigned> target_ptr(num_target_ranges + 1);
    target_ptr[0] = 0;
    std::vector<upair> source_ranges(num_box_pairs);
    std::vector<upair>::iterator source_ranges_curr = source_ranges.begin();

    // For all the target ranges
    for (unsigned k = 0; k < num_target_ranges; ++k) {
      // Copy the source ranges that interact with the kth target range
      unsigned i_begin = target_ranges[k].first;
      source_ranges_curr = std::copy(target2sources[i_begin].begin(),
                                     target2sources[i_begin].end(),
                                     source_ranges_curr);

      // Record the stop index
      target_ptr[k+1] = source_ranges_curr - source_ranges.begin();
    }

    // Sanity checking
    FMMTL_ASSERT(target_ptr.back() == source_ranges.size());
    FMMTL_ASSERT(source_ranges_curr == source_ranges.end());

    // Copy the source and target ranges into contiguous vectors
    std::vector<source_type> sources(c.source_begin(), c.source_end());
    std::vector<target_type> targets(c.target_begin(), c.target_end());

    return new S2T_Compressed<Kernel>(target_ranges,
                                      target_ptr,
                                      source_ranges,
                                      sources,
                                      targets);
  }
  */
};


#define FMMTL_KERNEL
#if defined(FMMTL_KERNEL)   // If compiling the .kern, include implementations
#  if defined(__CUDACC__)    // If compiling the .kern with nvcc
#    include "fmmtl/dispatch/S2T/S2T_Compressed.cu"
#  else                      // If not compiling the .kern with nvcc
#    include "fmmtl/dispatch/S2T/S2T_Compressed.cpp"
#  endif
#endif
#undef FMMTL_KERNEL
