#pragma once

#include <iostream>

#include "fmmtl/dispatch/S2T/S2T_Compressed.hpp"

template <typename Kernel>
S2T_Compressed<Kernel>::S2T_Compressed() {
}

template <typename Kernel>
S2T_Compressed<Kernel>::S2T_Compressed(
    std::vector<std::pair<unsigned,unsigned> >&,
    std::vector<unsigned>&,
    std::vector<std::pair<unsigned,unsigned> >&,
    const std::vector<source_type>&,
    const std::vector<target_type>&) {
}

template <typename Kernel>
S2T_Compressed<Kernel>::~S2T_Compressed() {
  // TODO: delete
}

template <typename Kernel>
void S2T_Compressed<Kernel>::execute(
    const Kernel&,
    const std::vector<charge_type>&,
    std::vector<result_type>&) {
  std::cerr << "ERROR: Calling unimplemented S2T_compressed CPU" << std::endl;
}

template <typename Kernel>
void S2T_Compressed<Kernel>::execute(
    const Kernel&,
    const std::vector<source_type>&,
    const std::vector<charge_type>&,
    const std::vector<target_type>&,
    std::vector<result_type>&) {
  std::cerr << "ERROR: Calling unimplemented S2T_compressed CPU" << std::endl;
}
