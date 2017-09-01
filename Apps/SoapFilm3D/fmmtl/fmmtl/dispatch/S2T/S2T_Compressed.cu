#pragma once

#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

#include "fmmtl/config.hpp"

#include "fmmtl/dispatch/S2T/S2T_Compressed.hpp"
#include "fmmtl/dispatch/S2T/S2T_Blocked_CSR.cu"

struct Data {
  unsigned num_sources;
  unsigned num_targets;
  unsigned num_blocks;
  Data(unsigned s, unsigned t, unsigned b)
      : num_sources(s),
        num_targets(t),
        num_blocks(b) {
  }
};


template <typename T>
inline T* gpu_new(unsigned n) {
  return thrust::raw_pointer_cast(thrust::device_malloc<T>(n));
}

template <typename Container>
inline typename Container::value_type* gpu_copy(const Container& c) {
  typedef typename Container::value_type c_value;
  // Allocate
  thrust::device_ptr<c_value> dptr = thrust::device_malloc<c_value>(c.size());
  // Copy
  //thrust::uninitialized_copy(c.begin(), c.end(), dptr);
  thrust::copy(c.begin(), c.end(), dptr);
  // Return
  return thrust::raw_pointer_cast(dptr);
}

template <typename T>
inline void gpu_free(T* p) {
  thrust::device_free(thrust::device_pointer_cast<void>(p));
}

template <typename Kernel>
S2T_Compressed<Kernel>::S2T_Compressed()
    : data_(0) {
}

template <typename Kernel>
S2T_Compressed<Kernel>::S2T_Compressed(
    std::vector<std::pair<unsigned,unsigned> >& target_ranges,
    std::vector<unsigned>& source_range_ptrs,
    std::vector<std::pair<unsigned,unsigned> >& source_ranges,
    const std::vector<source_type>& sources,
    const std::vector<target_type>& targets)
    : data_(new Data(sources.size(), targets.size(), target_ranges.size())),
      target_ranges_(gpu_copy(target_ranges)),
      source_range_ptrs_(gpu_copy(source_range_ptrs)),
      source_ranges_(gpu_copy(source_ranges)),
      sources_(gpu_copy(sources)),
      targets_(gpu_copy(targets)) {
}

template <typename Kernel>
S2T_Compressed<Kernel>::~S2T_Compressed() {
  delete reinterpret_cast<Data*>(data_);
  gpu_free(target_ranges_);
  gpu_free(source_range_ptrs_);
  gpu_free(source_ranges_);
  gpu_free(sources_);
  gpu_free(targets_);
}

/** A functor that indexes an array as one type but returns another type */
template <typename T1, typename T2>
class tricky_cast {
  T1* a_;
 public:
  __host__ __device__
  tricky_cast(T1* a) : a_(a) {}
  __host__ __device__
  T2 operator[](unsigned blockidx) const {
    return *((T2*)(a_ + blockidx));
  }
};

template <typename Kernel>
void S2T_Compressed<Kernel>::execute(
    const Kernel& K,
    const std::vector<charge_type>& charges,
    std::vector<result_type>& results) {
  typedef Kernel kernel_type;
  typedef typename kernel_type::source_type source_type;
  typedef typename kernel_type::target_type target_type;
  typedef typename kernel_type::charge_type charge_type;
  typedef typename kernel_type::result_type result_type;

  // XXX: Using a device_vector here was giving "floating point exceptions"...
  // XXX: device_vector doesn't like the Vec?
  charge_type* d_charges = gpu_copy(charges);
  result_type* d_results = gpu_copy(results);

  Data* data = reinterpret_cast<Data*>(data_);

  // TODO: set tpb to ncrit
  const unsigned num_tpb    = 256;
  const unsigned num_blocks = data->num_blocks;

#if defined(FMMTL_DEBUG)
  std::cout << "Launching GPU Kernel: (blocks, threads/block) = ("
            << num_blocks << ", " << num_tpb << ")" << std::endl;
#endif

  typedef thrust::pair<unsigned,unsigned> upair;

  // Launch kernel <<<grid_size, block_size>>>
  blocked_p2p<num_tpb><<<num_blocks,num_tpb>>>(
      K,
      target_ranges_,
      tricky_cast<unsigned, upair>(source_range_ptrs_),
      source_ranges_,
      sources_,
      //thrust::raw_pointer_cast(d_charges.data()),
      d_charges,
      targets_,
      d_results);
      //thrust::raw_pointer_cast(d_results.data()));
  FMMTL_CUDA_CHECK;

  // Copy results back
  thrust::device_ptr<result_type> d_results_ptr = thrust::device_pointer_cast(d_results);
  thrust::copy(d_results_ptr, d_results_ptr + results.size(), results.begin());

  gpu_free(d_results);
  gpu_free(d_charges);
}


/** A functor that maps blockidx -> (target_begin,target_end) */
template <unsigned BLOCKDIM>
class block_range {
  unsigned N_;
 public:
  __host__ __device__
  block_range(unsigned N) : N_(N) {}
  __host__ __device__
  thrust::pair<unsigned,unsigned> operator[](unsigned blockidx) const {
    return thrust::make_pair(blockidx * BLOCKDIM,
                             min(blockidx * BLOCKDIM + BLOCKDIM, N_));
  }
};

/** A functor that returns a constant */
template <typename T>
class constant {
  T value_;
 public:
  __host__ __device__
  constant(T value) : value_(value) {}
  __host__ __device__
  T operator[](unsigned) const {
    return value_;
  }
};

template <typename Kernel>
void
S2T_Compressed<Kernel>::execute(const Kernel& K,
                                const std::vector<source_type>& s,
                                const std::vector<charge_type>& c,
                                const std::vector<target_type>& t,
                                std::vector<result_type>& r) {
  typedef Kernel kernel_type;
  typedef typename kernel_type::source_type source_type;
  typedef typename kernel_type::target_type target_type;
  typedef typename kernel_type::charge_type charge_type;
  typedef typename kernel_type::result_type result_type;

  source_type* d_sources = gpu_copy(s);
  charge_type* d_charges = gpu_copy(c);
  target_type* d_targets = gpu_copy(t);
  result_type* d_results = gpu_copy(r);

  // XXX: device_vector doesn't like our vector?
  //thrust::device_vector<source_type> d_sources(s);
  //thrust::device_vector<charge_type> d_charges(c);
  //thrust::device_vector<target_type> d_targets(t);
  //thrust::device_vector<result_type> d_results(r);

  const unsigned num_tpb    = 256;
  const unsigned num_blocks = (t.size() + num_tpb - 1) / num_tpb;

#if defined(FMMTL_DEBUG)
  std::cout << "Launching GPU Kernel: (blocks, threads/block) = ("
            << num_blocks << ", " << num_tpb << ")" << std::endl;
#endif

  typedef thrust::pair<unsigned,unsigned> upair;

  // Launch kernel <<<grid_size, block_size>>>
  blocked_p2p<num_tpb><<<num_blocks, num_tpb>>>(
      K,
      block_range<num_tpb>(t.size()),
      constant<upair>(upair(0,1)),
      constant<upair>(upair(0,s.size())),
      d_sources,
      d_charges,
      d_targets,
      d_results);
      //thrust::raw_pointer_cast(d_sources.data()),
      //thrust::raw_pointer_cast(d_charges.data()),
      //thrust::raw_pointer_cast(d_targets.data()),
      //thrust::raw_pointer_cast(d_results.data()));
  FMMTL_CUDA_CHECK;

  // Copy results back and assign
  thrust::device_ptr<result_type> d_results_ptr = thrust::device_pointer_cast(d_results);
  thrust::copy(d_results_ptr, d_results_ptr + r.size(), r.begin());

  gpu_free(d_sources);
  gpu_free(d_charges);
  gpu_free(d_targets);
  gpu_free(d_results);
}
