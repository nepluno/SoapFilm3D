#pragma once
/** @file MortonCoder.hpp
 * @brief Define the MortonCoder class for Z-order-curve values, aka Morton
 *   codes.
 */

#include <cstdint>
#include <climits>

#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/tree/BoundingBox.hpp"

#include "fmmtl/numeric/bits.hpp"

namespace fmmtl {

/** @class MortonCoder
 * @brief Class representing Z-order-curve values, aka Morton codes.
 *
 * The Z-order curve is a space-filling curve: a one-dimensional curve that
 * fills a multi-dimensional space. Space-filling curves offer advantages for
 * representing points in ND space. Points near each other in 3D space are
 * likely to have close Morton codes! So we can store points in a map with
 * Z-order value as key, then iterate over nearby Z-order values to fetch
 * points nearby in space.
 *
 * Unfortunately, it's impossible to reduce ND space to 1D space perfectly:
 * there are some discontinuities in any space-filling curve mapping. But the
 * mapping is still an awesome tool, and some optimizations (see BIGMIN on the
 * Wikipedia page linked below) make them very effective.
 *
 * The MortonCoder class encapsulates a BoundingBox and can be used to translate
 * between spatial points and Morton codes relative to that BoundingBox.
 *
 * A single Morton code corresponds to a rectangular volume within that
 * BoundingBox, called its <em>cell</em>. Each side of the BoundingBox is
 * divided into 2^L equal-sized cells, for a total of 8^L cells.
 *
 * Read more about the Z-order curve here:
 * http://en.wikipedia.org/wiki/Z-order_curve
 *
 * This class computes maps box numbers to point and visa-versa
 * with respect to a bounding box and the number of equal-volume boxes (2^DL).
 * These mappings are performed in O(1) time.
 */
template <unsigned DIM>
struct MortonCoder {
  typedef Vec<DIM,double> point_type;

  // Using a 32-bit unsigned int for the code_type,
  // means we will resolve 32 1D levels, 16 2D levels, 10 3D levels, 8 4D levels.
  typedef uint32_t code_type;

  /** The number of bits per dimension = the maximum number of levels */
  static constexpr unsigned levels() {
    return std::numeric_limits<code_type>::digits / DIM;
  }
  /** The number of cells per side of the bounding box (2^L). */
  static constexpr uint64_t cells_per_side() {
    return uint64_t(1) << levels();
  }
  /** One more than the largest code (2^DL) */
  static constexpr uint64_t end_code() {
    return uint64_t(1) << (DIM * levels());
  }

  /** Construct a MortonCoder with a bounding box. */
  MortonCoder(const BoundingBox<point_type>& bb)
      : pmin_(bb.min()),
        cell_size_((bb.max() - bb.min()) / cells_per_side()) {
    cell_size_ *= (1.0 + 1e-14);  // Inclusive bounding box by small extension
    FMMTL_ASSERT(!bb.empty());
  }

  /** Return the MortonCoder's bounding box. */
  BoundingBox<point_type> bounding_box() const {
    return BoundingBox<point_type>(pmin_, pmin_ + cell_size_ * cells_per_side());
  }

  /** Return the bounding box of the cell with Morton code @a c.
   * @pre c < end_code */
  BoundingBox<point_type> cell(code_type c) const {
    point_type pmin = pmin_ + cell_size_ * deinterleave(c);
    return BoundingBox<point_type>(pmin, pmin + cell_size_);
  }

  /** Return the bounding box of the box with Morton codes @a cmin and @a cmax.
   * @pre cmin,cmax < end_code */
  BoundingBox<point_type> cell(code_type cmin, code_type cmax) const {
    return BoundingBox<point_type>(pmin_ + cell_size_ * deinterleave(cmin),
                                   pmin_ + cell_size_ *(deinterleave(cmax)+1));
  }

  /** Return the center of the Morton coders' bounding box */
  point_type center() const {
    return pmin_ + cell_size_ * (cells_per_side()/2);
  }

  /** Return the center of the box with Morton codes @a cmin and @a cmax.
   * @pre cmin,cmax < end_code */
  point_type center(code_type cmin, code_type cmax) const {
    return pmin_ + (cell_size_/2)*(deinterleave(cmax)+deinterleave(cmin)+1);
  }

  /** Return the Morton code of Point @a p.
   * @pre bounding_box().contains(@a p)
   * @post cell(result).contains(@a p) */
  code_type code(point_type s) const {
    s -= pmin_;
    s /= cell_size_;
    //for (unsigned i = 0; i < DIM; ++i) FMMTL_ASSERT(s[i] < cells_per_side);
    // Interleave the bits of the s[0], s[1], ...
    return interleave(s);
  }

  //private:
  /** The minimum of the MortonCoder bounding box. */
  point_type pmin_;
  /** The extent of a single cell. */
  point_type cell_size_;

  /** Spreads the bits of a number for interleaving */
  static inline unsigned spread_bits(unsigned x);

  /** Interleave the bits of s[0], s[1], ... */
  static inline code_type interleave(const point_type& s) {
    code_type code = code_type(0);
    for (unsigned i = 0; i < DIM; ++i) {
      FMMTL_ASSERT(unsigned(s[i]) < cells_per_side());
      code |= spread_bits(unsigned(s[i])) << i;
    }
    return code;
  }

  /** Compact the bits of a number for deinterleaving */
  static inline unsigned compact_bits(unsigned x);

  /** Deinterleave the bits from @a c into a point. */
  static inline point_type deinterleave(code_type c) {
    typedef typename point_type::value_type value_type;
    point_type p;
    for (unsigned i = 0; i < DIM; ++i)
      p[i] = value_type(compact_bits(c >> i));
    return p;
  }
};


/**************************
 ***** Generic Coder ******
 **************************/

#if 1
/** Spread the bits of an (unsigned) code_type by
 * inserting DIM-1 0s between each bit. (see bits/spread_bits_)
 */
template <unsigned DIM>
inline
typename MortonCoder<DIM>::code_type
MortonCoder<DIM>::spread_bits(typename MortonCoder<DIM>::code_type x) {
  return spread_bits_<DIM>(x);
}

/** Compact the bits of an (unsigned) code_type by
 * keeping every DIM-1 bit. (see bits/compact_bits_)
 */
template <unsigned DIM>
inline
typename MortonCoder<DIM>::code_type
MortonCoder<DIM>::compact_bits(typename MortonCoder<DIM>::code_type x) {
  return compact_bits_<DIM>(x);
}
#endif

/***************************
 ***** Specializations *****
 ***************************/
// Mostly for explication -- the above generates equivalent code

#if 1
/** Spread the bits of a 32-bit number (null op)
 * @param x 32-bit integer
 * @return 32-bit integer
 */
template <>
inline unsigned MortonCoder<1>::spread_bits(unsigned x) {
  return x;
}

/** Compact the bits of a 32-bit number (null op)
 * @param x 32-bit integer
 * @return 32-bit integer
 */
template <>
inline unsigned MortonCoder<1>::compact_bits(unsigned x) {
  return x;
}

/** Spread the bits of a 16-bit number so that there are 0s
 *  in between each bit.
 * @param x 16-bit integer
 * @return integer of form 0b0X0X0X0X0X0X0X0X0X0X0X0X0X0X0X0X,
 * where the X's are the original bits of @a x
 */
template <>
inline unsigned MortonCoder<2>::spread_bits(unsigned x) {
  x = (x | (x << 8)) & 0b00000000111111110000000011111111;
  x = (x | (x << 4)) & 0b00001111000011110000111100001111;
  x = (x | (x << 2)) & 0b00110011001100110011001100110011;
  x = (x | (x << 1)) & 0b01010101010101010101010101010101;
  return x;
}

/** Compact the bits of a 32-bit number by keeping only every-other bit
 * @param x 32-bit integer
 * @return integer of the form 0b0000000000000000XXXXXXXXXXXXXXXX,
 * where the X's are every-other bit of @a x
 */
template <>
inline unsigned MortonCoder<2>::compact_bits(unsigned x) {
  x &= 0b01010101010101010101010101010101;
  x = (x | (x >> 1)) & 0b00110011001100110011001100110011;
  x = (x | (x >> 2)) & 0b00001111000011110000111100001111;
  x = (x | (x >> 4)) & 0b00000000111111110000000011111111;
  x = (x | (x >> 8)) & 0b00000000000000001111111111111111;
  return x;
}

/** Spread the bits of a 10-bit number so that there are two 0s
 *  in between each bit.
 * @param x 10-bit integer
 * @return integer of form 0b0000X00X00X00X00X00X00X00X00X00X,
 * where the X's are the original bits of @a x
 */
template <>
inline unsigned MortonCoder<3>::spread_bits(unsigned x) {
  x = (x | (x << 16)) & 0b00000011000000000000000011111111;
  x = (x | (x <<  8)) & 0b00000011000000001111000000001111;
  x = (x | (x <<  4)) & 0b00000011000011000011000011000011;
  x = (x | (x <<  2)) & 0b00001001001001001001001001001001;
  return x;
}

/** Compact the bits of a 32-bit number by keeping only every third bit
 * @param x 32-bit integer
 * @return integer of the form 0b0000000000000000000000XXXXXXXXXX,
 * where the X's are every third bit of @a x
 */
template <>
inline unsigned MortonCoder<3>::compact_bits(unsigned x) {
  x &= 0b00001001001001001001001001001001;
  x = (x | (x >>  2)) & 0b00000011000011000011000011000011;
  x = (x | (x >>  4)) & 0b00000011000000001111000000001111;
  x = (x | (x >>  8)) & 0b00000011000000000000000011111111;
  x = (x | (x >> 16)) & 0b00000000000000000000001111111111;
  return x;
}

/** Spread the bits of a 8-bit number so that there are three 0s
 *  in between each bit.
 * @param x 8-bit integer
 * @return integer of form 0b000X000X000X000X000X000X000X000X,
 * where the X's are the original bits of @a x
 */
template <>
inline unsigned MortonCoder<4>::spread_bits(unsigned x) {
  x = (x | (x << 12)) & 0b00000000000011110000000000001111;
  x = (x | (x <<  6)) & 0b00000011000000110000001100000011;
  x = (x | (x <<  3)) & 0b00010001000100010001000100010001;
  return x;
}

/** Compact the bits of a 32-bit number by keeping only every fourth bit
 * @param x 32-bit integer
 * @return integer of the form 0b0000000000000000000000XXXXXXXXXX,
 * where the X's are every fourth bit of @a x
 */
template <>
inline unsigned MortonCoder<4>::compact_bits(unsigned x) {
  x &= 0b00010001000100010001000100010001;
  x = (x | (x >>  3)) & 0b00000011000000110000001100000011;
  x = (x | (x >>  6)) & 0b00000000000011110000000000001111;
  x = (x | (x >> 12)) & 0b00000000000000000000000011111111;
  return x;
}
#endif

} // end namespace fmmtl
