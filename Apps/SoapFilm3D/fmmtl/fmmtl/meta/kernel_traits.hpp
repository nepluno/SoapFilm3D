#pragma once

#include <boost/iterator/iterator_adaptor.hpp>

#include <iostream>

#include "fmmtl/meta/func_traits.hpp"


template <typename Kernel>
struct KernelTraits {
 public:
  typedef KernelTraits<Kernel>                    self_type;
  typedef Kernel                                  kernel_type;

  typedef typename kernel_type::source_type       source_type;
  typedef typename kernel_type::target_type       target_type;
  typedef typename kernel_type::charge_type       charge_type;
  typedef typename kernel_type::kernel_value_type kernel_value_type;
  typedef typename kernel_type::result_type       result_type;

  // Kernel evaluation operator, K(t,s)
  HAS_MEM_FUNC(HasEvalOp,
               kernel_value_type, operator(),
               const target_type&, const source_type&);
  static const bool has_eval_op = HasEvalOp<Kernel>::value;
  // Kernel transpose, K.transpose(K(t,s))
  HAS_MEM_FUNC(HasTranspose,
               kernel_value_type, transpose,
               const kernel_value_type&);
  static const bool has_transpose = HasTranspose<Kernel>::value;

 protected:
    // A dummy iterator adaptor to check for templated vectorized methods
  template <typename T>
  struct dumb_iterator : public boost::iterator_adaptor<dumb_iterator<T>,T*> {};

  typedef dumb_iterator<source_type> source_iterator;
  typedef dumb_iterator<charge_type> charge_iterator;
  typedef dumb_iterator<target_type> target_iterator;
  typedef dumb_iterator<result_type> result_iterator;

 public:
  HAS_MEM_FUNC(HasS2Tsymm,
               void, S2T,
               source_iterator, source_iterator, charge_iterator,
               target_iterator, target_iterator, charge_iterator,
               result_iterator, result_iterator);
  static const bool has_vector_S2T_symm = HasS2Tsymm<Kernel>::value;
  HAS_MEM_FUNC(HasS2Tasymm,
               void, S2T,
               source_iterator, source_iterator, charge_iterator,
               target_iterator, target_iterator, result_iterator);
  static const bool has_vector_S2T_asymm = HasS2Tasymm<Kernel>::value;

  friend std::ostream& operator<<(std::ostream& s, const self_type& traits) {
    s << "has_eval_op: "           << traits.has_eval_op           << std::endl;
    s << "has_transpose: "         << traits.has_transpose         << std::endl;
    s << "has_vector_S2T_symm: "   << traits.has_vector_S2T_symm   << std::endl;
    s << "has_vector_S2T_asymm: "  << traits.has_vector_S2T_asymm;
    return s;
  }
};


#include "fmmtl/meta/dimension.hpp"

// Expansion traits
template <typename Expansion>
struct ExpansionTraits
    : public KernelTraits<typename Expansion::kernel_type> {
  typedef Expansion                               expansion_type;
  typedef typename Expansion::kernel_type         kernel_type;
  typedef ExpansionTraits<expansion_type>         self_type;
  typedef KernelTraits<kernel_type>               super_type;

  typedef typename super_type::source_type        source_type;
  typedef typename super_type::target_type        target_type;
  typedef typename super_type::charge_type        charge_type;
  typedef typename super_type::kernel_value_type  kernel_value_type;
  typedef typename super_type::result_type        result_type;

  typedef typename super_type::source_iterator    source_iterator;
  typedef typename super_type::charge_iterator    charge_iterator;
  typedef typename super_type::target_iterator    target_iterator;
  typedef typename super_type::result_iterator    result_iterator;

  typedef typename expansion_type::point_type     point_type;
  static const std::size_t dimension = fmmtl::dimension<point_type>::value;

  typedef typename expansion_type::multipole_type multipole_type;
  typedef typename expansion_type::local_type     local_type;

  // Initializers
  HAS_MEM_FUNC(HasInitMultipole,
               void, init_multipole,
               multipole_type&, const point_type&, unsigned);
  static const bool has_init_multipole = HasInitMultipole<Expansion>::value;
  HAS_MEM_FUNC(HasInitLocal,
               void, init_local,
               local_type&, const point_type&, unsigned);
  static const bool has_init_local = HasInitLocal<Expansion>::value;

  // S2P, T2P
  HAS_MEM_FUNC(HasS2P,
               point_type, S2P,
               const source_type&);
  static const bool has_S2P = HasS2P<Expansion>::value;
  HAS_MEM_FUNC(HasT2P,
               point_type, T2P,
               const target_type&);
  static const bool has_T2P = HasT2P<Expansion>::value;

  // S2M
  HAS_MEM_FUNC(HasScalarS2M,
               void, S2M,
               const source_type&, const charge_type&,
               const point_type&, multipole_type&);
  static const bool has_scalar_S2M = HasScalarS2M<Expansion>::value;
  HAS_MEM_FUNC(HasVectorS2M,
               void, S2M,
               source_iterator, source_iterator, charge_iterator,
               const point_type&, multipole_type&);
  static const bool has_vector_S2M = HasVectorS2M<Expansion>::value;
  static const bool has_S2M = (has_scalar_S2M || has_vector_S2M);

  // S2L
  HAS_MEM_FUNC(HasScalarS2L,
               void, S2L,
               const source_type&, const charge_type&,
               const point_type&, local_type&);
  static const bool has_scalar_S2L = HasScalarS2L<Expansion>::value;
  HAS_MEM_FUNC(HasVectorS2L,
               void, S2L,
               source_iterator, source_iterator, charge_iterator,
               const point_type&, local_type&);
  static const bool has_vector_S2L = HasVectorS2L<Expansion>::value;
  static const bool has_S2L = (has_scalar_S2L || has_vector_S2L);

  // M2M
  HAS_MEM_FUNC(HasM2M,
               void, M2M,
               const multipole_type&, multipole_type&, const point_type&);
  static const bool has_M2M = HasM2M<Expansion>::value;

  // M2L
  HAS_MEM_FUNC(HasM2L,
               void, M2L,
               const multipole_type&, local_type&, const point_type&);
  static const bool has_M2L = HasM2L<Expansion>::value;

  // MAC
  HAS_MEM_FUNC(HasDynMAC,
               bool, MAC,
               const multipole_type&, const local_type&);
  static const bool has_dynamic_MAC = HasDynMAC<Expansion>::value;

  // L2L
  HAS_MEM_FUNC(HasL2L,
               void, L2L,
               const local_type&, local_type&, const point_type&);
  static const bool has_L2L = HasL2L<Expansion>::value;

  // M2T
  HAS_MEM_FUNC(HasScalarM2T,
               void, M2T,
               const multipole_type&, const point_type&,
               const target_type&, result_type&);
  static const bool has_scalar_M2T = HasScalarM2T<Expansion>::value;
  HAS_MEM_FUNC(HasVectorM2T,
               void, M2T,
               const multipole_type&, const point_type&,
               target_iterator, target_iterator, result_iterator);
  static const bool has_vector_M2T = HasVectorM2T<Expansion>::value;
  static const bool has_M2T = (has_scalar_M2T || has_vector_M2T);

  // L2T
  HAS_MEM_FUNC(HasScalarL2T,
               void, L2T,
               const local_type&, const point_type&,
               const target_type&, result_type&);
  static const bool has_scalar_L2T = HasScalarL2T<Expansion>::value;
  HAS_MEM_FUNC(HasVectorL2T,
               void, L2T,
               const local_type&, const point_type&,
               target_iterator, target_iterator, result_iterator);
  static const bool has_vector_L2T = HasVectorL2T<Expansion>::value;
  static const bool has_L2T = (has_scalar_L2T || has_vector_L2T);

  friend std::ostream& operator<<(std::ostream& s, const self_type& traits) {
    s << static_cast<super_type>(traits)                     << std::endl;
    s << "has_init_multipole: " << traits.has_init_multipole << std::endl;
    s << "has_init_local: "     << traits.has_init_local     << std::endl;
    s << "has_S2M: "            << traits.has_S2M            << std::endl;
    s << "  has_scalar_S2M: "   << traits.has_scalar_S2M     << std::endl;
    s << "  has_vector_S2M: "   << traits.has_vector_S2M     << std::endl;
    s << "has_S2L: "            << traits.has_S2L            << std::endl;
    s << "  has_scalar_S2L: "   << traits.has_scalar_S2L     << std::endl;
    s << "  has_vector_S2L: "   << traits.has_vector_S2L     << std::endl;
    s << "has_M2M: "            << traits.has_M2M            << std::endl;
    s << "has_M2L: "            << traits.has_M2L            << std::endl;
    s << "has_L2L: "            << traits.has_L2L            << std::endl;
    s << "has_M2T: "            << traits.has_M2T            << std::endl;
    s << "  has_scalar_M2T: "   << traits.has_scalar_M2T     << std::endl;
    s << "  has_vector_M2T: "   << traits.has_vector_M2T     << std::endl;
    s << "has_L2T: "            << traits.has_M2T            << std::endl;
    s << "  has_scalar_L2T: "   << traits.has_scalar_L2T     << std::endl;
    s << "  has_vector_L2T: "   << traits.has_vector_L2T     << std::endl;
    s << "has_dynamic_MAC: "    << traits.has_dynamic_MAC;
    return s;
  }
};

#define FMMTL_IMPORT_KERNEL_TRAITS(K)                                     \
  typedef typename KernelTraits<K>::kernel_type        kernel_type;       \
  typedef typename KernelTraits<K>::kernel_value_type  kernel_value_type; \
  typedef typename KernelTraits<K>::source_type        source_type;       \
  typedef typename KernelTraits<K>::target_type        target_type;       \
  typedef typename KernelTraits<K>::charge_type        charge_type;       \
  typedef typename KernelTraits<K>::result_type        result_type

#define FMMTL_IMPORT_EXPANSION_TRAITS(E)                                     \
  typedef typename ExpansionTraits<E>::kernel_type        kernel_type;       \
  typedef typename ExpansionTraits<E>::kernel_value_type  kernel_value_type; \
  typedef typename ExpansionTraits<E>::source_type        source_type;       \
  typedef typename ExpansionTraits<E>::target_type        target_type;       \
  typedef typename ExpansionTraits<E>::charge_type        charge_type;       \
  typedef typename ExpansionTraits<E>::result_type        result_type;       \
  typedef typename ExpansionTraits<E>::expansion_type     expansion_type;    \
  typedef typename ExpansionTraits<E>::multipole_type     multipole_type;    \
  typedef typename ExpansionTraits<E>::local_type         local_type;        \
  typedef typename ExpansionTraits<E>::point_type         point_type
