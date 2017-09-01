#ifndef _MATRIX3_H_
#define _MATRIX3_H_

#include <cstring>

#include "cvec3t.h"

template <class F> class Mat3T 
{
public:
  // make this an identity matrix
  Mat3T&  setIdentity() {     
    memset(_m,0, 9*sizeof(F));
    (*this)(0,0) = (*this)(1,1) = (*this)(2,2)  = 1.0f;
    return *this;
  }
  // returns an identity matrix  
  static Mat3T Identity() {
    return Mat3T().setIdentity();
  }

  static F doublecontraction(const Mat3T& A, const Mat3T& B) { 
    F fn = 0; 
    for(int i = 0; i < 9; i++) fn += A._m[i]*B._m[i]; 
    return fn;
  }

  Mat3T() { 
    setIdentity();
  }

  // create from elements
  Mat3T(F m00, F m01, F m02,
        F m10, F m11, F m12,
        F m20, F m21, F m22) { 
        Mat3T& M= *this;
        M(0,0) = m00; M(0,1) = m01; M(0,2) = m02; 
        M(1,0) = m10; M(1,1) = m11; M(1,2) = m12; 
        M(2,0) = m20; M(2,1) = m21; M(2,2) = m22; 
      }
  Mat3T(const Mat3T& mat) { 
    memcpy(_m,mat._m,9*sizeof(F));
  }
  
  // cast to array
  // entries are in column-major order, i.e. first 3 elements are the first column
  operator F*() 
  { return _m; }
  // cast to const array
  operator const F*() const { 
    return _m;
  }
  // conversion from array to matrix: require to be explicit to avoid
  // unexpected implicit casts 
  explicit Mat3T(F* m ) { 
    memcpy(_m, m,9*sizeof(F));
  }

#if 0
  // set to matrix for translation transform
  Mat3T& setTranslation( const CVec3T<F>& trans ) { 
    setIdentity(); 
    (*this)(0,3) = trans.x();
    (*this)(1,3) = trans.y();
    (*this)(2,3) = trans.z();
    return *this;
  }
  
  static Mat3T Translation( const CVec3T<F>& trans ) {
    return Mat3T().setTranslation(trans);
  }
#endif
  
  // set to matrix for nonuniform scale transform
  Mat3T& setScale( const CVec3T<F>& scale ) { 
    setIdentity(); 
    (*this)(0,0) = scale.x();
    (*this)(1,1) = scale.y();
    (*this)(2,2) = scale.z();
    return *this;
  }
  
  static Mat3T Scale( const CVec3T<F>& scale ) {
    return Mat3T().setScale(scale);
  }
  // set to a rotation matrix for axis v and angle a; unlike OpenGL the angle is in radians, not degrees!
  Mat3T& setRotation(F a,const CVec3T<F>& v) {
    //Vec3T<F> u = v.dir();
    CVec3T<F> u = v/v.l2();
    F u1 = u.x();
    F u2 = u.y();
    F u3 = u.z();
   
    Mat3T U = Mat3T(u1*u1, u1*u2, u1*u3,
                    u2*u1, u2*u2, u2*u3,
                    u3*u1, u3*u2, u3*u3 );
    Mat3T S = Mat3T( 0, -u3,  u2,  
                    u3,   0, -u1, 
                   -u2,  u1,   0);
  
    (*this) = U + (Identity() - U) * cos(a)  + S * sin(a);
#if 0
    (*this)(3,3) = 1.0;
#endif
    return *this;
  }

  static Mat3T Rotation( F a,const CVec3T<F>& v) {
    return Mat3T().setRotation(a,v);
  }

  
  // (i,j) element access
  // this is the only function for which the internal storage order matters
  F operator()(int i, int j) const { 
    return _m[3*j + i ];
  }
  
  F& operator()(int i, int j)  { 
    return _m[3*j + i ];
  }
  
  Mat3T& operator=(const Mat3T& mat) { 
    memcpy(_m, mat._m, 9*sizeof(F));
    return *this;
  }

  // extract column
  CVec3T<F> col(int i) { 
    return CVec3T<F>( (*this)(0,i), (*this)(1,i), (*this)(2,i));
  }
  // extract row
  CVec3T<F> row(int i) { 
    return CVec3T<F>( (*this)(i,0), (*this)(i,1), (*this)(i,2));
  }
  
  
// set the matrix to the inverse of M; return true if inverse exists  
  bool setInverse(const Mat3T& M) {
    Mat3T A;
    int i, j, k;
    F V;
    
    A = M;
    setIdentity();
   
    
    for (i = 0; i < 3; i++) {
      V = A(i,i);              /* Find the new pivot. */
      k = i;
      for (j = i + 1; j < 3; j++) 
      if (fabs(A(j,i)) > fabs(V)) {
       /* Find maximum on col i, row i+1..n */
       V = A(j,i);
       k = j;
      }
      j = k;
    

      F tmp;
      if (i != j)
      for (k = 0; k < 3; k++) {
       tmp = A(i,k); A(i,k) = A(j,k); A(j,k) = tmp;
       tmp = (*this)(i,k); (*this)(i,k) = (*this)(j,k); (*this)(j,k) = tmp;
      }
    
    
      for (j = i + 1; j < 3; j++) {   /* Eliminate col i from row i+1..n. */
      if(A(j,i) != F(0)) {
       V = A(j,i) / A(i,i);
       
       for (k = 0; k < 3; k++) {
        A(j,k)    -= V * A(i,k);
        (*this)(j,k) -= V * (*this)(i,k);
       }
      }
      
      }
    }
   
    for (i = 2; i >= 0; i--) {             /* Back Substitution. */
      if (A(i,i) == F(0))
      return false;             /* Error. */
    
      for (j = 0; j < i; j++) {   /* Eliminate col i from row 1..i-1. */
      V = A(j,i) / A(i,i);
      
      for (k = 0; k < 3; k++) {
       /* A[j][k] -= V * A[i][k]; */
       (*this)(j,k) -= V * (*this)(i,k);
      }
      }
    }
    
    for (i = 0; i < 3; i++)        /* Normalize the inverse Matrix. */
      for (j = 0; j < 3; j++)
      (*this)(i,j) /= A(i,i);
    
    return true;
  }  
  
  // return inverse; WARNING: returns garbage if the matrix is not invertible
  Mat3T inverse() const { 
  Mat3T M; M.setInverse(*this);
  return M;
  }
// set the matrix to transpose of M
  Mat3T& setTranspose(const Mat3T& M )  { 
  for(int i = 0; i < 3; i++) 
    for(int j = 0; j < 3; j++) 
    (*this)(i,j)= M(j,i);
  return (*this);
  }

  Mat3T transpose() const { 
  Mat3T M; M.setTranspose(*this);
  return M;  
  }
  // matrix vector product Mv, vector regarded as a column  
  CVec3T<F> operator* (const CVec3T<F>& v) const { 
   const Mat3T& M = *this;
   

   return CVec3T<F>( 
        v.x() * M(0,0) + v.y() * M(0,1) + v.z() * M(0,2) ,
        v.x() * M(1,0) + v.y() * M(1,1) + v.z() * M(1,2) ,
        v.x() * M(2,0) + v.y() * M(2,1) + v.z() * M(2,2) 
        );
  }
  // scale a matrix
  Mat3T operator* (F s) const {
   const Mat3T& M = *this;
   return Mat3T(     s*M(0,0), s*M(0,1), s*M(0,2), 
         s*M(1,0), s*M(1,1), s*M(1,2),
         s*M(2,0), s*M(2,1), s*M(2,2)
           );
  }
  // multiply two matrices  
  Mat3T operator* (const Mat3T& M2) const {     
   const Mat3T& M1 = *this;
   return Mat3T(
      M1(0,0)*M2(0,0) + M1(0,1)*M2(1,0) + M1(0,2)*M2(2,0),
      M1(0,0)*M2(0,1) + M1(0,1)*M2(1,1) + M1(0,2)*M2(2,1),
      M1(0,0)*M2(0,2) + M1(0,1)*M2(1,2) + M1(0,2)*M2(2,2),

             
      M1(1,0)*M2(0,0) + M1(1,1)*M2(1,0) + M1(1,2)*M2(2,0),
      M1(1,0)*M2(0,1) + M1(1,1)*M2(1,1) + M1(1,2)*M2(2,1),
      M1(1,0)*M2(0,2) + M1(1,1)*M2(1,2) + M1(1,2)*M2(2,2),
      
      M1(2,0)*M2(0,0) + M1(2,1)*M2(1,0) + M1(2,2)*M2(2,0),
      M1(2,0)*M2(0,1) + M1(2,1)*M2(1,1) + M1(2,2)*M2(2,1),
      M1(2,0)*M2(0,2) + M1(2,1)*M2(1,2) + M1(2,2)*M2(2,2)
             );
  }
  
  Mat3T operator+ (const Mat3T& M2) const {
  const Mat3T& M1 = *this;
  return Mat3T(
         M1(0,0)+M2(0,0), M1(0,1)+M2(0,1), M1(0,2)+M2(0,2),
         M1(1,0)+M2(1,0), M1(1,1)+M2(1,1), M1(1,2)+M2(1,2),
         M1(2,0)+M2(2,0), M1(2,1)+M2(2,1), M1(2,2)+M2(2,2) );
  }
  
  Mat3T operator- (const Mat3T& M2) const {
   const Mat3T& M1 = *this;
   return Mat3T(
          M1(0,0)-M2(0,0), M1(0,1)-M2(0,1), M1(0,2)-M2(0,2),
          M1(1,0)-M2(1,0), M1(1,1)-M2(1,1), M1(1,2)-M2(1,2),
          M1(2,0)-M2(2,0), M1(2,1)-M2(2,1), M1(2,2)-M2(2,2)
          );
   
  }
  
  F trace() const {
    return _m[0] + _m[4] + _m[8];
  }

  // Frobenius norm, i.e. sqrt of the sum of the squares of the entries also defined as trace(M^2)
  F frobnorm() const { 
    F fn = 0; 
    for(int i = 0; i < 9; i++) fn += _m[i]*_m[i]; 
    return sqrt(fn);
  }

  // Determinant
  F det() const {
      const Mat3T& M = *this;
      return   M(0,0) * (M(1,1) * M(2,2) - M(1,2) * M(2,1))
             - M(0,1) * (M(1,0) * M(2,2) - M(1,2) * M(2,0))
             + M(0,2) * (M(1,0) * M(2,1) - M(1,1) * M(2,0));
  }

 private:
  // data stored in column major order for OpenGL compatibility, i.e. first for elements are the first column of the matrix
  F _m[9];
};
template <class F> inline Mat3T<F> operator*(F s, const Mat3T<F>& M) { 
  return M*s;
}
template <class F> inline std::ostream& operator<< ( std::ostream& os, const Mat3T<F>& M ) { 
  os << "[ " << M(0,0)  << " " << M(0,1)  << " " << M(0,2)  << " " <<  "; ";
  os         << M(1,0)  << " " << M(1,1)  << " " << M(1,2)  << " " <<  "; ";
  os         << M(2,0)  << " " << M(2,1)  << " " << M(2,2)  << " " <<  "; ";
  os << "] ";   
return os;
}

template <class F>
    inline Mat3T<F> outerProd( const CVec3T<F>& c1, const CVec3T<F>& c2 )
{ return Mat3T<F>( c1.x() * c2.x(), c1.x() * c2.y(), c1.x() * c2.z(),
                   c1.y() * c2.x(), c1.y() * c2.y(), c1.y() * c2.z(),
                   c1.z() * c2.x(), c1.z() * c2.y(), c1.z() * c2.z() ); }

#endif
