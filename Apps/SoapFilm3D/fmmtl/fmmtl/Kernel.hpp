#pragma once
/** @file Kernel
 *
 * CRTP base classes for Kernels.
 * This allows the future extension of Kernels with methods that
 * need not be implemented within the Kernel.
 * Additionally, the traits system could be incorporated here to give the
 * library improved access and insight to these classes.
 *
 * At the moment, these may not be fully necessary and are subject to change.
 */

#include "fmmtl/config.hpp"

#include "fmmtl/dispatch/S2T/S2T_Compressed.hpp"

namespace fmmtl {


template <class DerivedKernel>
struct Kernel {
  // TODO: Do the template instantiations for linking automatically...?
  // TODO: Vectorized and/or GPU evaluations for Kernels?
};

// Template instantiations for external compilation and linking
// TODO: Remove?
#define FMMTL_KERNEL_EXTRAS(kernel) template class S2T_Compressed<kernel>

} // end namespace fmmtl
