#pragma once

#include <iterator>
#include <vector>


/** Maps boxes and box iterators to data and data iterators
 */
template <typename T, typename Tree>
struct BoxBind {
  typedef std::vector<T> container_type;
  typedef typename container_type::iterator iterator;
  typedef typename container_type::const_iterator const_iterator;

  container_type data;

  BoxBind(const Tree&, unsigned size)
      : data(size) {
  }

  T& operator[](const typename Tree::box_type& box) {
    return data[box.index()];
  }
  const T& operator[](const typename Tree::box_type& box) const {
    return data[box.index()];
  }

  iterator operator[](const typename Tree::box_iterator& bi) {
    return std::begin(data) + (*bi).index();
  }
  const_iterator operator[](const typename Tree::box_iterator& bi) const {
    return std::begin(data) + (*bi).index();
  }
};

template <typename T, typename Tree>
BoxBind<T,Tree> make_box_binding(const Tree& tree) {
  return {tree, tree.boxes()};
}


/** Maps bodies and body iterators to data and data iterators
 */
template <typename T, typename Tree>
struct BodyBind {
  typedef std::vector<T> container_type;
  typedef typename container_type::iterator iterator;
  typedef typename container_type::const_iterator const_iterator;

  //const Tree& tree;
  container_type data;

  struct BodyDataRange {
    iterator b, e;
    iterator       begin()       { return b; }
    const_iterator begin() const { return b; }
    iterator       end()         { return e; }
    const_iterator end()   const { return e; }
    std::size_t    size()  const { return end() - begin(); }
  };

  // Construct without permuting data or initialization data
  BodyBind(std::size_t n)
      : data(n) {
  }
  // Construct by permuting some initialization data based on the tree
  template <typename Iterator>
  BodyBind(const Tree& tree, Iterator data_it)
      : data(tree.body_permute(data_it, tree.body_begin()),
             tree.body_permute(data_it, tree.body_end())) {
  }

  iterator       begin()       { return std::begin(data); }
  const_iterator begin() const { return std::begin(data); }
  iterator       end()         { return std::end(data); }
  const_iterator end()   const { return std::end(data); }

  T& operator[](const typename Tree::body_type& body) {
    return data[body.index()];
  }
  const T& operator[](const typename Tree::body_type& body) const {
    return data[body.index()];
  }

  iterator operator[](const typename Tree::body_iterator& bi) {
    return std::begin(data) + bi.index();
  }
  const_iterator operator[](const typename Tree::body_iterator& bi) const {
    return std::begin(data) + bi.index();
  }

  BodyDataRange operator[](const typename Tree::box_type& box) {
    return {operator[](box.body_begin()), operator[](box.body_end())};
  }
};

// Specialization for iterators
template <typename Tree, typename Iterator>
BodyBind<typename std::iterator_traits<Iterator>::value_type, Tree>
make_body_binding(const Tree& tree, Iterator range) {
  return {tree, range};
}

// Any Range with std::begin
template <typename Tree, typename Range>
auto
make_body_binding(const Tree& tree, const Range& range)
    -> decltype(make_body_binding(tree, std::begin(range))) {
  return make_body_binding(tree, std::begin(range));
}

// Any Type without initial data (avoid permutation, just allocate)
template <typename Type, typename Tree>
BodyBind<Type, Tree>
make_body_binding(const Tree& tree) {
  return {tree.bodies()};
}
