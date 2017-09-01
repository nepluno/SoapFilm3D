#pragma once

namespace fmmtl {

template <typename T>
struct dimension;

template <>
struct dimension<double> {
  const static std::size_t value = 1;
};

template <>
struct dimension<float> {
  const static std::size_t value = 1;
};

} // end namepsace fmmtl
