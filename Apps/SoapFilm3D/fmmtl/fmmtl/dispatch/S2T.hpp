#pragma once
/** @file S2T.hpp
 * @brief Dispatch methods for the S2T stage
 *
 */

#include <iterator>
#include <type_traits>

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

#if !defined(P2P_BLOCK_SIZE)
#  define P2P_BLOCK_SIZE 128
#endif
#if !defined(P2P_NUM_THREADS)
#  define P2P_NUM_THREADS std::thread::hardware_concurrency()
#endif

// Definition of direct block-evaluation schemes
// TODO: Parallel dispatching option, trivial iterator detection
namespace fmmtl {
namespace detail {

/** Dual-Evaluation dispatch when K.transpose does not exist */
template <typename Kernel,
          typename Source, typename Charge,
          typename Target, typename Result>
inline
typename std::enable_if<!KernelTraits<Kernel>::has_transpose>::type
symm_eval(const Kernel& K,
          const Source& p1, const Charge& c1, Result& r1,
          const Target& p2, const Charge& c2, Result& r2)
{
  r1 += K(p1,p2) * c2;
  r2 += K(p2,p1) * c1;
}

/** Dual-Evaluation dispatch when K.transpose does exist */
template <typename Kernel,
          typename Source, typename Charge,
          typename Target, typename Result>
inline
typename std::enable_if<KernelTraits<Kernel>::has_transpose>::type
symm_eval(const Kernel& K,
          const Source& p1, const Charge& c1, Result& r1,
          const Target& p2, const Charge& c2, Result& r2)
{
  typedef typename KernelTraits<Kernel>::kernel_value_type kernel_value_type;

  const kernel_value_type k12 = K(p1,p2);
  r1 += k12 * c2;
  r2 += K.transpose(k12) * c1;
}

/** Asymmetric block P2P evaluation */
template <typename Kernel,
          typename SourceIter, typename ChargeIter,
          typename TargetIter, typename ResultIter>
inline void
block_eval(const Kernel& K,
           SourceIter s_first, SourceIter s_last, ChargeIter c_first,
           TargetIter t_first, TargetIter t_last, ResultIter r_first)
{
  typedef typename std::iterator_traits<SourceIter>::value_type source_type;
  typedef typename std::iterator_traits<TargetIter>::value_type target_type;
  typedef typename std::iterator_traits<ChargeIter>::value_type charge_type;
  typedef typename std::iterator_traits<ResultIter>::value_type result_type;

  for ( ; t_first != t_last; ++t_first, ++r_first) {
    const target_type& t = *t_first;
    result_type& r       = *r_first;

    SourceIter si = s_first;
    ChargeIter ci = c_first;
    for ( ; si != s_last; ++si, ++ci)
      r += K(t,*si) * (*ci);
  }
}

/** Symmetric off-diagonal block P2P evaluation */
template <typename Kernel,
          typename SourceIter, typename ChargeIter,
          typename TargetIter, typename ResultIter>
inline void
block_eval(const Kernel& K,
           SourceIter p1_first, SourceIter p1_last,
           ChargeIter c1_first, ResultIter r1_first,
           TargetIter p2_first, TargetIter p2_last,
           ChargeIter c2_first, ResultIter r2_first)
{
  typedef typename std::iterator_traits<SourceIter>::value_type source_type;
  typedef typename std::iterator_traits<TargetIter>::value_type target_type;
  typedef typename std::iterator_traits<ChargeIter>::value_type charge_type;
  typedef typename std::iterator_traits<ResultIter>::value_type result_type;

  for ( ; p1_first != p1_last; ++p1_first, ++c1_first, ++r1_first) {
    const source_type& pi = *p1_first;
    const charge_type& ci = *c1_first;
    result_type& ri       = *r1_first;

    TargetIter p2i = p2_first;
    ChargeIter c2i = c2_first;
    ResultIter r2i = r2_first;
    for ( ; p2i != p2_last; ++p2i, ++c2i, ++r2i) {
      const target_type& pj = *p2i;
      const charge_type& cj = *c2i;
      result_type& rj       = *r2i;

      symm_eval(K, pi, ci, ri, pj, cj, rj);
    }
  }
}

/** Symmetric diagonal block P2P evaluation */
template <typename Kernel,
          typename SourceIter, typename ChargeIter, typename ResultIter>
inline static void
block_eval(const Kernel& K,
           SourceIter p_first, SourceIter p_last,
           ChargeIter c_first, ResultIter r_first)
{
  typedef typename std::iterator_traits<SourceIter>::value_type source_type;
  typedef typename std::iterator_traits<ChargeIter>::value_type charge_type;
  typedef typename std::iterator_traits<ResultIter>::value_type result_type;

  SourceIter ipi = p_first;
  ChargeIter ici = c_first;
  ResultIter iri = r_first;

  for ( ; ipi != p_last; ++ipi, ++ici, ++iri) {
    const source_type& pi = *ipi;
    const charge_type& ci = *ici;
    result_type& ri       = *iri;

    // The diagonal element
    ri += K(pi,pi) * ci;

    // The off-diagonal elements
    SourceIter ipj = p_first;
    ChargeIter icj = c_first;
    ResultIter irj = r_first;
    for ( ; ipj != ipi; ++ipj, ++icj, ++irj) {
      const source_type& pj = *ipj;
      const charge_type& cj = *icj;
      result_type& rj       = *irj;

      symm_eval(K, pi, ci, ri, pj, cj, rj);
    }
  }
}


} // end namespace detail
} // end namespace fmmtl



// TODO: namespace and fix

struct S2T {

  //////////////////////////////////////
  /////// Context Dispatchers //////////
  //////////////////////////////////////

  struct ONE_SIDED {};
  struct TWO_SIDED {};

  /** Asymmetric S2T
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& source,
                          const typename Context::target_box_type& target,
                          const ONE_SIDED&)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "S2T:"
              << "\n  " << source
              << "\n  " << target << std::endl;
#endif
    FMMTL_LOG("S2T 2box asymm");

    fmmtl::detail::block_eval(c.kernel(),
                              c.source_begin(source), c.source_end(source),
                              c.charge_begin(source),
                              c.target_begin(target), c.target_end(target),
                              c.result_begin(target));
  }

  /** Symmetric S2T
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& box1,
                          const typename Context::target_box_type& box2,
                          const TWO_SIDED&)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "S2T:"
              << "\n  " << box1
              << "\n  " << box2 << std::endl;
    std::cout << "S2T:"
              << "\n  " << box2
              << "\n  " << box1 << std::endl;
#endif
    FMMTL_LOG("S2T 2box symm");

    fmmtl::detail::block_eval(c.kernel(),
                              c.source_begin(box1), c.source_end(box1),
                              c.charge_begin(box1), c.result_begin(box1),
                              c.target_begin(box2), c.target_end(box2),
                              c.charge_begin(box2), c.result_begin(box2));
  }

  /** Symmetric S2T
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& box)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "S2T:"
              << "\n  " << box << std::endl;
#endif
    FMMTL_LOG("S2T 1box symm");

    fmmtl::detail::block_eval(c.kernel(),
                              c.source_begin(box), c.source_end(box),
                              c.charge_begin(box), c.result_begin(box),
                              c.target_begin(box), c.target_end(box),
                              c.charge_begin(box), c.result_begin(box));
  }
};
