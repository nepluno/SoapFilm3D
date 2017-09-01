#pragma once
/** @file NDTree
 * @brief General class representing a {1D,2D,3D,4D}-Tree.
 */

#include <vector>
#include <algorithm>

#include <iostream>
#include <iomanip>

#include <boost/range.hpp>
#include <boost/iterator/permutation_iterator.hpp>

#include "fmmtl/tree/util/CountedProxyIterator.hpp"

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/tree/BoundingBox.hpp"
#include "fmmtl/tree/MortonCoder.hpp"

namespace fmmtl {
using boost::has_range_iterator;


/** In-place bucket sort using counting sort
 *
 * @param[in] first,last Iterator pair to the sequence to be bucketed
 * @param[in] num_buckets The number of buckets to be used
 * @param[in] map Mapping of elements in [first,last) to [0,num_buckets)
 * @returns vector of iterators:
 *          [result[i],result[i+1]) is the range of the ith bucket.
 *
 * @pre For all i in [first,last), 0 <= map(*i) < num_buckets.
 * @post For all i and j such that first <= i < j < last,
 *       then 0 <= map(*i) <= map(*j) < num_buckets.
 * @tparam Iterator models a random-access mutable iterator
 */
template <typename Iterator, typename BucketMap>
std::vector<Iterator> bucket_sort(Iterator first, Iterator last,
                                  unsigned num_buckets, BucketMap map) {
  typedef typename std::iterator_traits<Iterator>::value_type value_type;

  std::vector<Iterator> start(num_buckets+1, first);

  for ( ; first != last; ++first)
    ++start[1 + map(*first)];

  for (unsigned k = 2; k <= num_buckets; ++k)
    start[k] += (start[k-1] - start[0]);

  std::vector<Iterator> end = start;

  // For each bin
  for (unsigned curr_bin = 0; curr_bin < num_buckets; ++curr_bin) {
    first = start[curr_bin];
    last  = end[1+curr_bin];

    while (first != last) {
      value_type& v = *first++;

      // Swap until an element comes back to this bin
      unsigned bin;
      while (curr_bin != (bin = map(v)))
        std::swap(v, *start[bin]++);
    }
  }

  return end;
}


//! Class for tree structure
template <unsigned DIM>
class NDTree {
  // Predeclarations
  struct Box;
  struct Body;
  struct BoxData;

  // Morton coder to order the points
  typedef MortonCoder<DIM> coder_type;
  typedef typename coder_type::code_type code_type;

  // The type of this tree
  typedef NDTree<DIM> tree_type;

 public:
  //! The type of indices and integers in this tree
  typedef unsigned size_type;

  //! The spacial point type used for centers and extents
  typedef Vec<DIM,double> point_type;

  //! Public type declarations
  typedef Box           box_type;
  typedef Body          body_type;
  using box_iterator  = CountedProxyIterator<Box,  const NDTree, size_type>;
  using body_iterator = CountedProxyIterator<Body, const NDTree, size_type>;

 private:
  // Tree representation

  // Morton coder to transform points wrt a bounding box
  coder_type coder_;
  // Morton code for each point
  std::vector<code_type> mc_;
  // Permutation: permute_[i] is the current idx of originally ith point
  std::vector<size_type> permute_;
  // Vector of data describing a box
  std::vector<BoxData> box_data_;
  // level_offset_[i] and level_offset_[i+1] is the start and end of level i
  std::vector<size_type> level_offset_;

  struct BoxData {
    // Index of the first body in this box
    size_type body_begin_;
    // Index of one-past-last body in this box
    size_type body_end_;

    // Index of first child box
    size_type cbox_begin_;
    // Index of one-past-last child box
    size_type cbox_end_;

    // Index of the parent box
    size_type parent_;

    // Precomputed center
    point_type center_;

    BoxData(size_type bb, size_type be, size_type cbb, size_type cbe,
            size_type p, const point_type& c)
        : body_begin_(bb), body_end_(be),
          cbox_begin_(cbb), cbox_end_(cbe), parent_(p),
          center_(c) {
    }

    bool is_leaf() const {
      return cbox_begin_ == size_type(-1);
    }
  };

  struct Body {
    /** Construct an invalid Body */
    Body() {}
    //! The original order this body was seen
    size_type number() const {
      return tree_->permute_[idx_];
    }
    //! The current order of this body
    size_type index() const {
      return idx_;
    }
    code_type morton_index() const {
      return tree_->mc_[idx_];
    }
   private:
    size_type idx_;
    tree_type* tree_;
    friend body_iterator;
    Body(size_type idx, const tree_type* tree)
        : idx_(idx), tree_(const_cast<tree_type*>(tree)) {
      FMMTL_ASSERT(idx_ < tree_->size());
    }
    friend class NDTree;
  };

  // A tree-aligned box
  struct Box {
    typedef typename tree_type::box_iterator  box_iterator;
    typedef typename tree_type::body_iterator body_iterator;

    //! Construct an invalid Box
    Box() {}
    //! The index of this box
    size_type index() const {
      return idx_;
    }
    //! The level of this box (root level is 0)
    size_type level() const {
      const auto it = std::upper_bound(tree_->level_offset_.begin(),
                                       tree_->level_offset_.end(),
                                       idx_);
      return it - tree_->level_offset_.begin() - 1;
    }

    //! The dimension of each side of this box
    point_type extents() const {
      const BoundingBox<point_type> bb = tree_->coder_.bounding_box();
      return (bb.max() - bb.min()) / (1 << level());
    }
    //! The squared radius of this box
    double radius_sq() const {
      return norm_2_sq(extents()) / 4.0;
    }
    //! The center of this box
    const point_type& center() const {
      return data().center_;
    }

    //! The parent box of this box
    Box parent() const {
      return Box(data().parent_, tree_);
    }

    //! True if this box is a leaf and has no children
    bool is_leaf() const {
      return data().is_leaf();
    }
    //! The begin iterator to the child boxes contained in this box
    box_iterator child_begin() const {
      FMMTL_ASSERT(!is_leaf());
      return box_iterator(data().cbox_begin_, tree_);
    }
    //! The end iterator to the child boxes contained in this box
    box_iterator child_end() const {
      FMMTL_ASSERT(!is_leaf());
      return box_iterator(data().cbox_end_, tree_);
    }
    //! The number of children this box has
    size_type num_children() const {
      return std::distance(child_begin(), child_end());
    }

    //! The begin iterator to the bodies contained in this box
    body_iterator body_begin() const {
      return body_iterator(data().body_begin_, tree_);
    }
    //! The end iterator to the bodies contained in this box
    body_iterator body_end() const {
      return body_iterator(data().body_end_, tree_);
    }
    //! The number of bodies this box contains
    size_type num_bodies() const {
      return std::distance(body_begin(), body_end());
    }

    //! Equality comparison operator
    bool operator==(const Box& b) const {
      FMMTL_ASSERT(tree_ == b.tree_);
      return idx_ == b.idx_;
    }
    //! Comparison operator for std:: containers and algorithms
    bool operator<(const Box& b) const {
      FMMTL_ASSERT(tree_ == b.tree_);
      return idx_ < b.idx_;
    }

    //! Write a Box to an output stream
    friend std::ostream& operator<<(std::ostream& s,
                                    const box_type& b) {
      size_type num_bodies = b.num_bodies();
      size_type first_body = b.body_begin()->index();
      size_type last_body = first_body + num_bodies - 1;

      return s << "Box " << b.index()
               << " (L" << b.level() << ", P" << b.parent().index()
               << ", " << num_bodies << (num_bodies == 1 ? " body" : " bodies")
               << " " << first_body << "-" << last_body
               << "): " << b.center() << " - " << b.extents();
    }

   private:
    size_type idx_;
    tree_type* tree_;
    friend box_iterator;
    Box(size_type idx, const tree_type* tree)
        : idx_(idx), tree_(const_cast<tree_type*>(tree)) {
    }
    BoxData& data() const {
      return tree_->box_data_[idx_];
    }
    friend class NDTree;
  };

 public:

  /** Construct a tree encompassing a bounding box
   * and insert a range of points */
  template <typename Range>
  NDTree(const Range& rng, size_type n_crit = 256,
         typename std::enable_if<has_range_iterator<Range>::value>::type* = 0)
      : NDTree(rng.begin(), rng.end(), n_crit) {
  }

  /** Construct an tree encompassing a bounding box
   * and insert a range of points */
  template <typename PointIter>
  NDTree(PointIter first, PointIter last, size_type n_crit = 256)
      : coder_(get_boundingbox(first, last)) {
    insert(first, last, n_crit);
  }

  /** Return the Bounding Box that this NDTree encompasses */
  BoundingBox<point_type> bounding_box() const {
    return coder_.bounding_box();
  }

  /** Return the center of this NDTree */
  point_type center() const {
    return coder_.center();
  }

  /** The number of bodies contained in this tree */
  size_type size() const {
    return permute_.size();
  }
  /** The number of bodies contained in this tree */
  size_type bodies() const {
    return size();
  }

  /** The number of boxes contained in this tree */
  size_type boxes() const {
    return box_data_.size();
  }

  /** The number of boxes contained in level L of this tree */
  size_type boxes(size_type L) const {
    return level_offset_[L+1] - level_offset_[L];
  }

  /** The maximum level of any box in this tree */
  size_type levels() const {
    return level_offset_.size() - 1;
  }

  /** Returns true if the box is contained in this tree, false otherwise */
  bool contains(const box_type& box) const {
    return this == box.tree_;
  }
  /** Returns true if the body is contained in this tree, false otherwise */
  bool contains(const body_type& body) const {
    return this == body.tree_;
  }

  /** Return the root box of this tree */
  box_type root() const {
    return Box(0, this);
  }
  /** Return a box given its index */
  box_type box(const size_type idx) const {
    FMMTL_ASSERT(idx < box_data_.size());
    return Box(idx, this);
  }
  /** Return a body given its index */
  body_type body(const size_type idx) const {
    FMMTL_ASSERT(idx < size());
    return Body(idx, this);
  }
  /** Return an iterator to the first body in this tree */
  body_iterator body_begin() const {
    return body_iterator(0, this);
  }
  /** Return an iterator one past the last body in this tree */
  body_iterator body_end() const {
    return body_iterator(bodies(), this);
  }
  /** Return an iterator to the first box in this tree */
  box_iterator box_begin() const {
    return box_iterator(0, this);
  }
  /** Return an iterator one past the last box in this tree */
  box_iterator box_end() const {
    return box_iterator(boxes(), this);
  }
  /** Return an iterator to the first box at level L in this tree
   * @pre L < levels()
   */
  box_iterator box_begin(size_type L) const {
    FMMTL_ASSERT(L < levels());
    return box_iterator(level_offset_[L], this);
  }
  /** Return an iterator one past the last box at level L in this tree
   * @pre L < levels()
   */
  box_iterator box_end(size_type L) const {
    FMMTL_ASSERT(L < levels());
    return box_iterator(level_offset_[L+1], this);
  }

  template <typename RandomAccessIter>
  struct body_permuted_iterator {
    typedef typename std::vector<size_type>::const_iterator permute_iter;
    typedef boost::permutation_iterator<RandomAccessIter, permute_iter> type;
  };

  /** Tranform (permute) an iterator so its traversal follows the same order as
   * the bodies contained in this tree
   */
  template <typename RandomAccessIter>
  typename body_permuted_iterator<RandomAccessIter>::type
  body_permute(RandomAccessIter it, const body_iterator& bi) const {
    return boost::make_permutation_iterator(it, permute_.cbegin() + bi.index());
  }

  /** Tranform (permute) an iterator so its traversal follows the same order as
   * the bodies contained in this tree
   *
   * Specialized for bi = body_begin().
   */
  template <typename RandomAccessIter>
  typename body_permuted_iterator<RandomAccessIter>::type
  body_permute(RandomAccessIter it) const {
    return body_permute(it, body_begin());
  }

  /** Write an NDTree to an output stream */
  friend std::ostream& operator<<(std::ostream& s,
                                  const tree_type& t) {
    struct {
      std::ostream& print(std::ostream& s,
                          const box_type& box) {
        s << std::string(2*box.level(), ' ') << box;
        if (!box.is_leaf())
          for (auto ci = box.child_begin(); ci != box.child_end(); ++ci)
            print(s << "\n", *ci);
        return s;
      }
    } recursive_box;

    return recursive_box.print(s, t.root());
  }

 private:
  //! TODO: Make dynamic and public?
  //! Uses incremental bucket sorting
  template <typename PointIter>
  void insert(PointIter p_first, PointIter p_last, size_type NCRIT) {
    FMMTL_LOG("Tree Insert");

    // Create a code-idx pair vector
    typedef std::pair<code_type, size_type> code_pair;
    std::vector<code_pair> codes;
    // If iterators are random access, we can reserve space efficiently
    // Compile-time predicate!
    if (std::is_same<typename std::iterator_traits<PointIter>::iterator_category,
        std::random_access_iterator_tag>::value)
      codes.reserve(std::distance(p_first, p_last));

    // Initialize code-index pairs
    for (size_type idx = 0; p_first != p_last; ++p_first, ++idx)
      codes.emplace_back(coder_.code(*p_first), idx);

    // Allocate representation
    mc_.reserve(codes.size());
    permute_.reserve(codes.size());

    // Push the root box which contains all points
    box_data_.emplace_back(0, codes.size(), -1, -1, 0, center());
    level_offset_.push_back(0);
    level_offset_.push_back(1);

    const size_type max_children = (1 << DIM);

    // For each level
    for (size_type L = 1; L < level_offset_.size(); ++L) {

      // Define the bucketer for this next level
      const size_type shift = DIM * (coder_type::levels() - L);
      auto bucketer = [=] (const code_pair& v) {
        return (v.first >> shift) & (max_children-1);
      };

      // Construct the next level  TODO: 2 or 3 at a time
      for (size_type k = level_offset_[L-1]; k < level_offset_[L]; ++k) {

        // If this box has few enough points, go to next box
        if (box_data_[k].body_end_ - box_data_[k].body_begin_ <= NCRIT)
          continue;

        // Else, split this box

        // Get the box data
        auto code_begin = codes.begin() + box_data_[k].body_begin_;
        auto code_end   = codes.begin() + box_data_[k].body_end_;

        // Sort the points in this box into the "bucket" children
        auto off = bucket_sort(code_begin, code_end, max_children, bucketer);

        // Record the child begin idx
        box_data_[k].cbox_begin_ = box_data_.size();

        // For each bucket
        for (size_type c = 0; c < max_children; ++c) {
          // If this child contains points
          if (off[c+1] != off[c]) {
            // Add the child
            box_data_.emplace_back(off[c]   - codes.begin(),   // Body begin idx
                                   off[c+1] - codes.begin(),   // Body end idx
                                   -1, -1,                     // Children idxs
                                   k,                          // Parent idx
                                   get_center(off[c]->first, L)); // Center
          }
        }

        // Record the child end idx
        box_data_[k].cbox_end_ = box_data_.size();
      }

      // Record the level end
      if (box_data_.size() > level_offset_.back())
        level_offset_.push_back(box_data_.size());
    }

    // Extract the code, permutation vector, and sorted point
    for (auto& c : codes) {
      mc_.push_back(c.first);
      permute_.push_back(c.second);
    }
  }

  /** Get the center of the box that
   * morton code @a c is contained in at level @level
   */
  point_type get_center(code_type c, size_type level) {
    // Mask for boxes of this level
    code_type mask = (code_type(1) << (DIM*(coder_type::levels()-level))) - 1;
    return coder_.center(c & ~mask /*cmin*/, c | mask /*cmax*/);
  }

  template <typename PointIter>
  BoundingBox<point_type> get_boundingbox(PointIter first, PointIter last) {
    // Construct a bounding box
    BoundingBox<point_type> bb(first, last);
    // Determine the size of the maximum side
    point_type extents =  bb.max() - bb.min();
    double max_side = *std::max_element(extents.begin(), extents.end());
    double radius = (1.0+1e-6) * max_side / 2.0;
    point_type center =  (bb.max() + bb.min()) / 2;
    // Make it square and add some wiggle room   TODO: Generalize on square
    return BoundingBox<point_type>(center - radius, center + radius);
  }

  // Just making sure for now
  NDTree(const NDTree&) {};
  void operator=(const NDTree&) {};
};

} // end namespace fmmtl
