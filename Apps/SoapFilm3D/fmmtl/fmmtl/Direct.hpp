#pragma once
/** @file Direct.hpp
 * @brief Dispatch methods for direct matrix evaluation
 */

#include <type_traits>

#include <boost/range/has_range_iterator.hpp>

#include "fmmtl/dispatch/S2T.hpp"

namespace fmmtl {

/** Asymmetric matvec
 */
template <typename Kernel,
          typename SourceIter, typename ChargeIter,
          typename TargetIter, typename ResultIter>
inline void
direct(const Kernel& K,
       SourceIter s_first, SourceIter s_last,
       ChargeIter c_first,
       TargetIter t_first, TargetIter t_last,
       ResultIter r_first)
{
  detail::block_eval(K,
                     s_first, s_last, c_first,
                     t_first, t_last, r_first);
}

/** Symmetric matvec, off-diagonal block
 */
template <typename Kernel,
          typename SourceIter, typename ChargeIter, typename ResultIter>
inline void
direct(const Kernel& K,
       SourceIter p1_first, SourceIter p1_last,
       ChargeIter c1_first, ResultIter r1_first,
       SourceIter p2_first, SourceIter p2_last,
       ChargeIter c2_first, ResultIter r2_first)
{
  detail::block_eval(K,
                     p1_first, p1_last, c1_first, r1_first,
                     p2_first, p2_last, c2_first, r2_first);
}

/** Symmetric matvec, diagonal block
 */
template <typename Kernel,
          typename SourceIter, typename ChargeIter, typename ResultIter>
inline
typename std::enable_if<!boost::has_range_iterator<SourceIter>::value>::type
direct(const Kernel& K,
       SourceIter p_first, SourceIter p_last,
       ChargeIter c_first, ResultIter r_first)
{
  detail::block_eval(K,
                     p_first, p_last, c_first, r_first);
}

/** Convenience function for ranges
 */
template <typename Kernel,
          typename SourceRange, typename ChargeRange,
          typename TargetRange, typename ResultRange>
inline
typename std::enable_if<boost::has_range_iterator<SourceRange>::value>::type
direct(const Kernel& K,
       const SourceRange& s, const ChargeRange& c,
       const TargetRange& t, ResultRange& r)
{
  detail::block_eval(K,
                     s.begin(), s.end(), c.begin(),
                     t.begin(), t.end(), r.begin());
}

/** Convenience function for ranges
 */
template <typename Kernel,
          typename STRange, typename ChargeRange, typename ResultRange>
inline void
direct(const Kernel& K,
       const STRange& p, const ChargeRange& c, ResultRange& r)
{
  detail::block_eval(K,
                     p.begin(), p.end(), c.begin(), r.begin());
}


} // end namespace fmmtl
