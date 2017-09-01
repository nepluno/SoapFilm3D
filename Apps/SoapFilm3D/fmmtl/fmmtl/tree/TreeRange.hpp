#pragma once
/** Helper utilities for range-based for loops with trees
 */

namespace fmmtl {

template <class T>
struct BoxRange {
  const T& t;
  auto begin() -> decltype(t.box_begin()) {
    return t.box_begin();
  }
  auto end() -> decltype(t.box_end()) {
    return t.box_end();
  }
};

template <class T>
BoxRange<T> boxes(const T& t) {
  return {t};
}


template <class T>
struct BoxLevelRange {
  unsigned L;
  const T& t;
  auto begin() -> decltype(t.box_begin(L)) {
    return t.box_begin(L);
  }
  auto end() -> decltype(t.box_end(L)) {
    return t.box_end(L);
  }
};

template <class T>
BoxLevelRange<T> boxes(unsigned L, const T& t) {
  return {L,t};
}


template <class T>
struct ChildRange {
  const T& t;
  auto begin() -> decltype(t.child_begin()) {
    return t.child_begin();
  }
  auto end() -> decltype(t.child_end()) {
    return t.child_end();
  }
};

template <class T>
ChildRange<T> children(const T& t) {
  return {t};
}


template <class T>
struct BodyRange {
  const T& t;
  auto begin() -> decltype(t.body_begin()) {
    return t.body_begin();
  }
  auto end() -> decltype(t.body_end()) {
    return t.body_end();
  }
};

template <class T>
BodyRange<T> bodies(const T& t) {
  return {t};
}



} // end namespace fmmtl
