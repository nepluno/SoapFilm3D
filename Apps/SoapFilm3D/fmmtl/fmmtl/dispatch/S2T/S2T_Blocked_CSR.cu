#pragma once

#include <thrust/system/cuda/detail/detail/uninitialized.h>
using thrust::system::cuda::detail::detail::uninitialized_array;

#include "fmmtl/meta/kernel_traits.hpp"

/** CSR Blocked S2T in CUDA
 * @brief  Computes the kernel matrix-vector product using a blocked CSR-like
 *     format. Each target range has a range of source ranges to compute.
 *
 * @param[in] K                 The kernel to generate matrix elements.
 * @param[in] target_range      Maps blockIdx.x to pair<uint,uint> representing
                                 the [start,end) of targets for this threadblock.
 * @param[in] source_range_ptr  Maps blockIdx.x to pair<uint,uint> representing
                                 the [start,end) of the source ranges interacting
                                 the target range for this threadblock.
 * @param[in] source_range      Maps each index of the source_range_ptr range to
                                 a [start,end) of a source range to interact
 *                               with each target of this threadblock..
 *
 * @param[in]     source  Array of sources to index into
 * @param[in]     charge  Array of charges associated with sources to index into
 * @param[in]     target  Array of targets to index into
 * @param[in,out] result  Array of results associated with targets to accumulate
 *
 * @pre For all k, target_range[k].second - target_range[k].first <= blockDim.x
 *       -- One target/thread = each target range is smaller than the blocksize
 */
template <unsigned BLOCKDIM,
          typename Kernel,
          typename Indexable1,
          typename Indexable2,
          typename Indexable3>
__global__ void
blocked_p2p(const Kernel K,
            Indexable1 target_range,
            Indexable2 source_range_ptr,
            Indexable3 source_range,
            const typename KernelTraits<Kernel>::source_type* source,
            const typename KernelTraits<Kernel>::charge_type* charge,
            const typename KernelTraits<Kernel>::target_type* target,
                  typename KernelTraits<Kernel>::result_type* result) {
  typedef typename KernelTraits<Kernel>::source_type source_type;
  typedef typename KernelTraits<Kernel>::charge_type charge_type;
  typedef typename KernelTraits<Kernel>::target_type target_type;
  typedef typename KernelTraits<Kernel>::result_type result_type;

  typedef thrust::pair<unsigned, const unsigned> upair;

  // Allocate shared memory -- prevent initialization of non-POD
  __shared__ uninitialized_array<source_type,BLOCKDIM> sh_s;
  __shared__ uninitialized_array<charge_type,BLOCKDIM> sh_c;

  // Get the target range this block is responsible for
  upair t_range = target_range[blockIdx.x];
  // The target index this thread is responsible for
  t_range.first += threadIdx.x;

  // Get the range of source ranges this block is responsible for
  upair s_range_ptr = source_range_ptr[blockIdx.x];

  // Each thread is assigned to one target in the target range
  result_type r = result_type();
  target_type t = ((t_range.first < t_range.second)
                   ? target[t_range.first] : target_type());

  // For each source range
  for ( ; s_range_ptr.first < s_range_ptr.second; ++s_range_ptr.first) {
    // Get the source range
    upair s_range = source_range[s_range_ptr.first];

    // For each chuck of sources
    for ( ; s_range.first < s_range.second; s_range.first += BLOCKDIM) {

      // Read up to blockDim.x sources into shared memory
      unsigned n = min(s_range.second - s_range.first, BLOCKDIM);
      if (threadIdx.x < n) {
        sh_s[threadIdx.x] = source[s_range.first + threadIdx.x];
        sh_c[threadIdx.x] = charge[s_range.first + threadIdx.x];
      }
      __syncthreads();

      // Each target computes its interaction with each source in smem
      if (t_range.first < t_range.second) {
        do {  // Note that n >= 1 by @pre for (s_range.first ...
          --n;
          r += K(t, sh_s[n]) * sh_c[n];
        } while (n != 0);
      }
      __syncthreads();   // TODO: Unroll to prevent an extra __syncthreads()?
    }
  }

  if (t_range.first < t_range.second)
    result[t_range.first] += r;
}
