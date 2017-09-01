// -*- Mode: c++ -*-
// Adrian Secord, version 1/2006
// Carefully adapted from Denis' CVec3T class.


#ifndef	__CVEC2T_H__
#define	__CVEC2T_H__

#include <iostream>
#include <assert.h>
#include <math.h>

// lines in this file are up to 160 char long

// list of non-member operations:
// IO: operator>>, istream& operator<<
// F = scalar class, V = vector class 
// when used as argument const ref is implied
//
// univariate: void normalize(V&)
//             F    lenSq (V)
//	       F    len   (V)
//             F    l2    (V)
//	       F    l1    (V)
//             F    linfty(V)
//	       V    op-   (V)
//	       V    dir   (V) 
//             int  largestAbsComp(V)
//             int  smallestAbsComp(V)
//
// bivariate:  V&   op=    (V&,V)
//             V&   op+=   (V&,V) 
//             V&   op-=   (V&,V)
//             V&   op*=   (V&,F)
//             V&   op/=   (V&,F)
//             V    op+    (V, V)
//             V    op-    (V, V)
//             V    compMult(V, V)
//             V    compDiv(V, V)
//             V    op*    (V, F) 
//             V    op*    (F, V) 
//             V    op/    (V, F) 
//             V    max    (V, V)
//             V    min    (V, V)
//             F    dot    (V, V)
//             F    dist   (V, V) 
//             F    angle  (V, V) 
//             V    project(V, V)
//          bool    isCollinear(V,V)   
//
// trivariate: V    lerp   (V, V, F) 

// there are no epsilons used anywhere: we do test for precise 0 before
// dividing, and for the range for  acos; in general one may need 
// to test if vectors are close to being equal etc
// accuracy of such tests is problem depedent and does not belong 
// here 
// TODO define == and !=



template <class F>
class CVec2T {
public:
  enum{ X = 0, Y = 1 };

  CVec2T( void )                  {}
  CVec2T( const CVec2T& c )       { v[X] = c.v[X]; v[Y] = c.v[Y]; } // copy constructor
  CVec2T( const F& a, const F& b) { v[X] = a;      v[Y] = b;      } // vector from 2 numbers
  explicit CVec2T( const F* a )            { v[X] = a[X];   v[Y] = a[Y];   } // vector from an array 	

  // This breaks sorting vectors, because they get cast to pointers and then compared.
  operator const F*( void ) const { return &v[X]; }  // cast a vector to an array

  ~CVec2T( void ) {}

  template <class G>
  CVec2T& operator=(const CVec2T<G>& c2) { v[0] =  c2(0); v[1]  = c2(1); return (*this); }



  // access components  
        F& x( void )       { return v[X]; }  
  const F& x( void ) const { return v[X]; }
        F& y( void )       { return v[Y]; }  
  const F& y( void ) const { return v[Y]; }
  
        F& operator() (const unsigned int i)       { assert(i < 2); return v[i]; }
  const F& operator() (const unsigned int i) const { assert(i < 2); return v[i]; }
  
 protected:
  F v[2];
};


// ----------------------------------------------------------------------------------------------------

// input/output
template <class F> std::ostream& operator<<(std::ostream& os, const CVec2T<F>& c) {  return os << c(0) << " " << c(1); }
template <class F> std::istream& operator>>(std::istream& is,       CVec2T<F>& c) {  is >> c(0);  is >> c(1);  return is;}


// ----------------------------------------------------------------------------------------------------
// univariate operations; most could be implemented as member functions, but do
// everything as non-member functions for consistency to simplify formula translation 

template <class F>  void      normalize(CVec2T<F>& c      ) { F l = len(c); assert(l != F(0)); c = c/len(c); }     // make the vector unit length

template <class F>  F         lenSq    (const CVec2T<F>& c) { return c(0)*c(0) + c(1)*c(1); }                      // squared length
template <class F>  F         len      (const CVec2T<F>& c) { return sqrt(lenSq(c)); }                             // length
template <class F>  F         l2       (const CVec2T<F>& c) { return sqrt(lenSq(c)); }                             // synonym
template <class F>  F         l1       (const CVec2T<F>& c) { return fabs(c(0)) + fabs(c(1)); }                    // l1 norm |x| + |y|
template <class F>  F         linfty   (const CVec2T<F>& c) { return std::max(fabs(c(0)),fabs(c(1))); }            // l infinity norm max(|x|,|y|)

template <class F>  CVec2T<F> operator-(const CVec2T<F>& c) { return CVec2T<F>( -c(0), -c(1)); }                   // unary minus
template <class F>  CVec2T<F> dir      (const CVec2T<F>& c) { F l = len(c); assert(l != F(0)); return c/l;}        // unit vector of the same direction


// find the index of the largest and smallest components of the vector
template <class F>  int largestAbsComp (const CVec2T<F>& c1)  {
  F a= fabs(c1(0)), b = fabs(c1(1));
  if      (a >= b) return 0;
  else return 1;
}

template <class F>  int smallestAbsComp (const CVec2T<F>& c1)  {
  F a= fabs(c1(0)), b = fabs(c1(1));
  if      (a <= b) return 0;
  else return 1;
}

// ----------------------------------------------------------------------------------------------------

// using multiple template params below is a mechanism to facilitate generation of
// additional operations combining different types of vectors
// e.g. as needed to support automatic differentition in an efficient way
// (so that constant vectors do not have to be cast to variable vectors)
// however, these templates are never instantiated automatically, because
// overload on return type is not allowed 
// so if we have a template with three parameters F,G,H, for
// operation F op(G g, H h) 
// it is insufficient to resolve F1 vf,vf1, vf2; vf = op(vf1, vf2).
// for this to be resolved, we need an explicit template specialization;
// the specializations for all parameters of the same type is return is 
// added at the end
// If we want to allow operations F0 vf0; F1 vf1; F2 vf2; 
// vf0 = op(vf1,vf2) (assuming op(f1,f2) is defined and can be cast to F0, 
// we need either a partial specialization 
// template <class G, class H> F0 op(G g, H h) { return op<F0,G,H>(g,h); } 
// or (more safely) a complete specialization 
// F0 op(F1 g, F2 h){ return  op<F0,F1,F2>(g,h) }

// ----------------------------------------------------------------------------------------------------
// bivariate operations

//assignment shortcuts

template <class F, class G, class H> 
inline CVec2T<F>& operator+=(CVec2T<G>& c1, const CVec2T<H>& c2) { c1(0) += c2(0); c1(1) += c2(1); return c1; }
template <class F, class G, class H> 
inline CVec2T<F>& operator-=(CVec2T<G>& c1, const CVec2T<H>& c2) { c1(0) -= c2(0); c1(1) -= c2(1); return c1; }
template <class F, class G, class H> 
inline CVec2T<F>& operator*=(CVec2T<G>& c1, const H& s         ) { c1(0) *= s    ; c1(1) *= s    ; return c1; }
template <class F, class G, class H> 
inline CVec2T<F>& operator/=(CVec2T<G>& c1, const H& s         ) { assert(s!=F(0)); 
                                                                   c1(0) /= s    ; c1(1) /= s    ; return c1; }


// bivariate vector arithmetic operations +,-, component division and multiplication,
// multiplication by constant,division by constant, 

// a gcc problem: there is aliasing with a template operator defined in basic_string
// can be either avoided by requiring not to have using namespace std
// before this file is included, or using a different name

//template <class F, class G, class H> 
//inline CVec2T<F> operator+(const CVec2T<G>& c1, const CVec2T<H>& c2) { return CVec2T<F>( c1(0) + c2(0), c1(1) + c2(1), c1(2) + c2(2)); }
template <class F, class G, class H> 
inline CVec2T<F>    opplus(const CVec2T<G>& c1, const CVec2T<H>& c2) { return CVec2T<F>( c1(0) + c2(0), c1(1) + c2(1)); }
template <class F, class G, class H> 
inline CVec2T<F> operator-(const CVec2T<G>& c1, const CVec2T<H>& c2) { return CVec2T<F>( c1(0) - c2(0), c1(1) - c2(1)); }
template <class F, class G, class H> 
inline CVec2T<F> compMult (const CVec2T<G>& c1, const CVec2T<H>& c2) { return CVec2T<F>( c1(0) * c2(0), c1(1) * c2(1)); }
template <class F, class G, class H> 
inline CVec2T<F> compDiv  (const CVec2T<G>& c1, const CVec2T<H>& c2) { assert(c2(0) != F(0) && c2(1) != F(0));
                                                                       return CVec2T<F>( c1(0) / c2(0), c1(1) / c2(1)); }
template <class F, class G, class H> 
inline CVec2T<F> operator*(const CVec2T<G>& c1, const        H&  s ) { return CVec2T<F>( c1(0) * s,     c1(1) * s    ); }
template <class F, class G, class H> 
inline CVec2T<F> operator*(const        G&  s , const CVec2T<H>& c1) { return CVec2T<F>( s     * c1(0), s     * c1(1)); }
template <class F, class G, class H> 
inline CVec2T<F> operator/(const CVec2T<G>& c1, const        H&  s ) { assert(s != F(0) );
                                                                       return CVec2T<F>( c1(0) / s,     c1(1) / s    ); }

// componentwise max/min
template <class F, class G, class H> 
inline CVec2T<F> max      (const CVec2T<G>& c1, const CVec2T<H>& c2) { return CVec2T<F>( max(c1(0),c2(0)), max(c1(1),c2(1))); }
template <class F, class G, class H> 
inline CVec2T<F> min      (const CVec2T<G>& c1, const CVec2T<H>& c2) { return CVec2T<F>( min(c1(0),c2(0)), min(c1(1),c2(1))); }

// geometric operations: dot product, distance len(u-v), angle between vectors,
// collinearity query
// cross product, projection (u dot v) v/|v|^2, 

template <class F, class G, class H> 
inline F dot ( const CVec2T<G>& c1, const CVec2T<H>& c2 ) { return ( c1.x()* c2.x() + c1.y() * c2.y()); } 
template <class F, class G, class H> 
inline F dist( const CVec2T<G>& c1, const CVec2T<H>& c2 ) { return len(c1-c2); } // eliminate?


// angle between two vectors 0.. Pi; conservative version: asserts on
// argument to acos > 1, which may happen because of numerical innacuracy
// however calculations are rearranged to minimize probability of this
// not using truncation to [-1,1] because of issues with automatic differentiation
// probably the right approach to make this robust is to define a special
// version of acos for doubles which does the truncation

template <class F, class G, class H> 
inline F angle( const CVec2T<G>& c1, const CVec2T<H>& c2 ) { 
  F s   = lenSq(c1)*lenSq(c2); assert(s != F(0));
  F dtp = dot(c1,c2); 
  F dps = dtp*dtp/s; assert( dps <= F(1) && dps >= F(-1) );
  return acos(dtp > F(0)? sqrt( dps ):-sqrt(dps) );
}


template <class F, class G, class H> 
inline CVec2T<F> project( const CVec2T<G>& c1, const CVec2T<H>& c2 ) { return c2*dot(c1,c2)/lenSq(c2); }


// ----------------------------------------------------------------------------------------------------

// trivariate geometric operations: linear interpolation (1-t)u + tv, rotation of u
// around v by angle a,  triple product u dot (v cross w)

template <class F, class G, class H, class K>
inline CVec2T<F> lerp      ( const CVec2T<G>& c1, const CVec2T<H>& c2, const K& s )        { return c1 + s*(c2-c1); } //eliminate?

 

// ----------------------------------------------------------------------------------------------------
// specializations for all parameter types equal to the return type


template <class F>  inline CVec2T<F>& operator+= (      CVec2T<F>& c1, const CVec2T<F>& c2) { return  operator+=<F,F,F>( c1,c2); }
template <class F>  inline CVec2T<F>& operator-= (      CVec2T<F>& c1, const CVec2T<F>& c2) { return  operator-=<F,F,F>( c1,c2); }
template <class F>  inline CVec2T<F>& operator*= (      CVec2T<F>& c1, const F&         s ) { return  operator*=<F,F,F>( c1,s ); }
template <class F>  inline CVec2T<F>& operator/= (      CVec2T<F>& c1, const F&         s ) { return  operator/=<F,F,F>( c1,s ); }

template <class F>  inline CVec2T<F>  operator+  (const CVec2T<F>& c1, const CVec2T<F>& c2) { return      opplus<F,F,F>( c1,c2); }
template <class F>  inline CVec2T<F>  operator-  (const CVec2T<F>& c1, const CVec2T<F>& c2) { return   operator-<F,F,F>( c1,c2); }
template <class F>  inline CVec2T<F>  compMult   (const CVec2T<F>& c1, const CVec2T<F>& c2) { return   compMult <F,F,F>( c1,c2); }
template <class F>  inline CVec2T<F>  compDiv    (const CVec2T<F>& c1, const CVec2T<F>& c2) { return   compDiv  <F,F,F>( c1,c2); }
template <class F>  inline CVec2T<F>  operator*  (const CVec2T<F>& c1, const         F&  s) { return   operator*<F,F,F>( c1,s ); }
template <class F>  inline CVec2T<F>  operator*  (const F& s,          const CVec2T<F>& c1) { return   operator*<F,F,F>(  s,c1); }
template <class F>  inline CVec2T<F>  operator/  (const CVec2T<F>& c1, const         F&  s) { return   operator/<F,F,F>( c1,s ); }

template <class F>  inline CVec2T<F>  max        (const CVec2T<F>& c1, const CVec2T<F>& c2) { return         max<F,F,F>( c1,c2); }
template <class F>  inline CVec2T<F>  min        (const CVec2T<F>& c1, const CVec2T<F>& c2) { return         min<F,F,F>( c1,c2); }

template <class F>  inline F          dot        (const CVec2T<F>& c1, const CVec2T<F>& c2) { return         dot<F,F,F>( c1,c2); }
template <class F>  inline F          dist       (const CVec2T<F>& c1, const CVec2T<F>& c2) { return        dist<F,F,F>( c1,c2); }
template <class F>  inline F          angle      (const CVec2T<F>& c1, const CVec2T<F>& c2) { return       angle<F,F,F>( c1,c2); }
template <class F>  inline CVec2T<F>  project    (const CVec2T<F>& c1, const CVec2T<F>& c2) { return     project<F,F,F>( c1,c2); }
template <class F>  inline bool  isCollinear(const CVec2T<F>& c1, const CVec2T<F>& c2) { return isCollinear<F,F,F>( c1,c2); }

template <class F>  inline CVec2T<F>  lerp       (const CVec2T<F>& c1, const CVec2T<F>& c2, const F& s         ) { return lerp      <F,F,F,F>( c1,c2,s ); }


#endif	/* __CVEC2T_H__ */
