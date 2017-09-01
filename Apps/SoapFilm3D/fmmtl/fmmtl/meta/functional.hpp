#pragma once

namespace fmmtl {

// Identity function object
struct identity {
  template <typename T>
  constexpr auto operator()(T&& v) const noexcept
      -> decltype(std::forward<T>(v)) {
    return std::forward<T>(v);
  }
};

// Null function object
struct nullop {
  template <typename T>
  constexpr void operator()(T&&) const noexcept {}
};

}
