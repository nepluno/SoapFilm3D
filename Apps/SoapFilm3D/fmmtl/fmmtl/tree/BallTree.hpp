#pragma once
/** @file KDTree
 * @brief A binary tree constructed by splitting each box by a median.
 */

#include <vector>
#include <algorithm>
#include <type_traits>
#include <iterator>

#include <iostream>
#include <iomanip>

#include <boost/range.hpp>
#include <boost/iterator/permutation_iterator.hpp>

#include "fmmtl/tree/util/CountedProxyIterator.hpp"

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/tree/BoundingSphere.hpp"

#include "fmmtl/numeric/bits.hpp"

namespace fmmtl {
using boost::has_range_iterator;


//! Class for tree structure
template <unsigned DIM>
class BallTree {
  // Predeclarations
  struct Box;
  struct Body;
  struct BoxData;

  // The type of this tree
  typedef BallTree<DIM> tree_type;

 public:
  //! The type of indices and integers in this tree
  typedef unsigned size_type;

  //! The spacial point type used for centers and extents
  typedef Vec<DIM,double> point_type;

  //! Public type declarations
  typedef Box           box_type;
  typedef Body          body_type;
  using box_iterator  = CountedProxyIterator<Box,  const BallTree, size_type>;
  using body_iterator = CountedProxyIterator<Body, const BallTree, size_type>;

 private:
  // Tree representation

  // Permutation: permute_[i] is the current idx of originally ith point
  std::vector<size_type> permute_;
  // Vector of data describing a box
  std::vector<BoxData> box_data_;

  struct BoxData {
    typedef BoundingSphere<point_type> bounding_sphere_type;

    // Index of the first body in this box
    size_type body_begin_;
    // Index of one-past-last body in this box
    size_type body_end_;
    // Bounding sphere
    bounding_sphere_type bounding_sphere_;

    // Constructor
    BoxData(size_type bb, size_type be, const bounding_sphere_type& bs)
        : body_begin_(bb), body_end_(be), bounding_sphere_(bs) {}
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
   private:
    size_type idx_;
    tree_type* tree_;
    friend body_iterator;
    Body(size_type idx, const tree_type* tree)
        : idx_(idx), tree_(const_cast<tree_type*>(tree)) {
      FMMTL_ASSERT(idx_ < tree_->size());
    }
    friend class BallTree;
  };

  // A box of the tree
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
      return std::log2(idx_+1);  // XXX: Slow? do this with bits
    }

    //! The dimension of each side of this box
    point_type extents() const {
      point_type e;
      e.fill(2*std::sqrt(radius_sq()));
      return e;
    }
    //! The squared radius of this box
    double radius_sq() const {
      return data().bounding_sphere_.radius_sq();
    }
    //! The center of this box
    point_type center() const {
      return data().bounding_sphere_.center();
    }

    //! The parent box of this box
    Box parent() const {
      FMMTL_ASSERT(!(*this == tree_->root()));
      return Box((idx_-1)/2, tree_);
    }

    //! True if this box is a leaf and has no children
    bool is_leaf() const {
      return 2*idx_+1 >= tree_->box_data_.size();
    }
    //! The begin iterator to the child boxes contained in this box
    box_iterator child_begin() const {
      FMMTL_ASSERT(!is_leaf());
      return box_iterator(2*idx_+1, tree_);
    }
    //! The end iterator to the child boxes contained in this box
    box_iterator child_end() const {
      FMMTL_ASSERT(!is_leaf());
      return box_iterator(2*idx_+3, tree_);
    }
    //! The number of children this box has
    static constexpr size_type num_children() {
      return 2;
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
    friend std::ostream& operator<<(std::ostream& s, const box_type& b) {
      size_type num_bodies = b.num_bodies();
      size_type first_body = b.body_begin()->index();
      size_type last_body = first_body + num_bodies - 1;
      size_type parent_idx = b.index()==0 ? 0 : b.parent().index();

      return s << "Box " << b.index()
               << " (L" << b.level() << ", P" << parent_idx
               << ", " << num_bodies << (num_bodies == 1 ? " body" : " bodies")
               << " " << first_body << "-" << last_body
               << "): Center: " << b.center()
               << "; Radius: " << std::sqrt(b.radius_sq());
    }

   private:
    size_type idx_;
    tree_type* tree_;
    friend box_iterator;
    Box(size_type idx, const tree_type* tree)
        : idx_(idx), tree_(const_cast<tree_type*>(tree)) {
      FMMTL_ASSERT(idx_ < tree_->boxes());
    }
    inline BoxData& data() const {
      return tree_->box_data_[idx_];
    }
    friend class BallTree;
  };

 public:

  /** Construct a tree encompassing a bounding box
   * and insert a range of points */
  template <typename Range>
  BallTree(const Range& rng, size_type n_crit = 256,
           typename std::enable_if<has_range_iterator<Range>::value>::type* = 0)
      : BallTree(rng.begin(), rng.end(), n_crit) {
  }

  /** Construct an tree encompassing a bounding box
   * and insert a range of points */
  template <typename PointIter>
  BallTree(PointIter first, PointIter last, size_type n_crit = 256) {
    insert(first, last, n_crit);
  }

  /** Return the Bounding Box that this BallTree encompasses */
  BoundingSphere<point_type> bounding_sphere() const {
    return box_data_[0].bounding_sphere_;
  }

  /** Return the center of this BallTree */
  point_type center() const {
    return root().center();
  }

  /** The number of bodies contained in this tree */
  inline size_type size() const {
    return permute_.size();
  }
  /** The number of bodies contained in this tree */
  inline size_type bodies() const {
    return size();
  }

  /** The number of boxes contained in this tree */
  inline size_type boxes() const {
    return box_data_.size();
  }

  /** The number of boxes contained in level L of this tree */
  inline size_type boxes(size_type L) const {
    return (1 << L);
  }

  /** The maximum level of any box in this tree */
  inline size_type levels() const {
    return std::log2(boxes());  // XXX: Slow?
  }

  /** Returns true if the box is contained in this tree, false otherwise */
  inline bool contains(const box_type& box) const {
    return this == box.tree_;
  }
  /** Returns true if the body is contained in this tree, false otherwise */
  inline bool contains(const body_type& body) const {
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
    return box_iterator((1 << L) - 1, this);
  }
  /** Return an iterator one past the last box at level L in this tree
   * @pre L < levels()
   */
  box_iterator box_end(size_type L) const {
    FMMTL_ASSERT(L < levels());
    return box_iterator((1 << (L+1)) - 1, this);
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

  /** Write a BallTree to an output stream */
  friend std::ostream& operator<<(std::ostream& s, const BallTree& t) {
    struct {
      std::ostream& print(std::ostream& s, const box_type& box) {
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
  template <typename PointIter>
  void insert(PointIter p_first, PointIter p_last, size_type NCRIT) {
    FMMTL_LOG("BallTree Insert");

    assert(p_first != p_last);

    // Create a point-idx pair vector
    typedef typename std::iterator_traits<PointIter>::value_type point_i_type;

    // XXX: Generalize?
    static_assert(std::is_same<point_i_type, point_type>::value,
                  "PointIter value_type must be point_type");

    typedef std::pair<point_type, unsigned> point_t;

    std::vector<point_t> point;
    // If iterators are random access, we can reserve space efficiently
    if (std::is_same<typename std::iterator_traits<PointIter>::iterator_category,
        std::random_access_iterator_tag>::value)
      point.reserve(std::distance(p_first, p_last));

    // Copy to point-idx pair
    unsigned idx = 0;
    for (PointIter pi = p_first; pi != p_last; ++pi, ++idx) {
      point.emplace_back(*pi, idx);
    }
    permute_.reserve(point.size());

    // Helper std::pair projection operator
    struct pair2point {
      point_type& operator()(point_t& p) const { return p.first; }
    };
    using proj = boost::transform_iterator<pair2point, decltype(point.begin())>;

    // The number of leaf boxes that will be created
    // (Smallest power of two greater than or equal to ceil(N/NCRIT))
    unsigned leaves = ceil_pow_2((point.size() + NCRIT - 1) / NCRIT);
    unsigned levels = std::log2(leaves);

    // Reserve the number of boxes that will be added
    box_data_.reserve(2*leaves - 1);

    // Push the root sphere which contains all points
    box_data_.emplace_back(0, size_type(point.size()),
                           approx_bounding_sphere(proj(point.begin()),
                                                  proj(point.end())));

    // For every box that is created
    unsigned end_k = (1 << levels) - 1;
    for (unsigned k = 0; k < end_k; ++k) {
      // The range of points in this box
      auto p_begin = point.begin() + box_data_[k].body_begin_;
      auto p_end   = point.begin() + box_data_[k].body_end_;

      // Invariant of helper functions
      assert(p_begin != p_end);

      // Compute the furthest point from the center
      const point_type& p1 = furthest_point_from(proj(p_begin), proj(p_end),
                                                 box(k).center());
      // Compute the furthest point from point1
      const point_type& p2 = furthest_point_from(proj(p_begin), proj(p_end),
                                                 p1);

      // Partition on the median: balanced tree, but worse spheres
      auto p_mid = p_begin + (p_end - p_begin) / 2;
      point_type r = p2 - p1;
      std::nth_element(p_begin, p_mid, p_end,
                       [&](const point_t& a, const point_t& b) {
           return inner_prod(a.first-p1, r) < inner_prod(b.first-p1, r);
        });

      /*
      // Partition on proximity: better spheres, potentially unbalanced
      // Would need to account for the possibility of empty boxes
      auto p_mid = std::partition(p_begin, p_end, [&] (const point_t& a) {
          return norm_2_sq(a.first-p1) < norm_2_sq(a.first-p2);
        });
      */

      // Record the child boxes
      size_type mid = p_mid - point.begin();
      box_data_.emplace_back(box_data_[k].body_begin_, mid,
                             approx_bounding_sphere(proj(p_begin), proj(p_mid)));

      box_data_.emplace_back(mid, box_data_[k].body_end_,
                             approx_bounding_sphere(proj(p_mid), proj(p_end)));
    }

    // Assert no re-allocation
    assert(box_data_.size() <= 2*leaves-1);

    // Extract the permutation idx
    for (const point_t& p : point)
      permute_.push_back(p.second);
  }

  /** Calculate the approximate bounding sphere of a set of points
   * using the Points Closest to Furthest Pair approach.
   */
  template <typename PointIter>
  BoundingSphere<point_type> approx_bounding_sphere(PointIter first,
                                                    PointIter last) {
    auto n = std::distance(first, last);
    point_type center = std::accumulate(first, last, point_type{}) / n;

    double r_sq = std::accumulate(first, last, double{},
                                  [&] (double r_sq, point_type& p) {
                                    return std::max(r_sq, norm_2_sq(center-p));
                                  });

    return {center, r_sq};
  }

  /** Find the point in a set of points that is furthest away from a target point
   */
  template <typename PointIter>
  const point_type& furthest_point_from(PointIter first,
                                        PointIter last,
                                        const point_type& target) {
    using pair = std::pair<double, const point_type*>;
    return *std::accumulate(first, last, pair{},
                            [&] (const pair& r, const point_type& p) {
                              return std::max(r, pair{norm_2_sq(target-p), &p});
                            }).second;
  }

  // Just making sure for now
  BallTree(const BallTree&) {};
  void operator=(const BallTree&) {};
};


} // end namespace fmmtl
