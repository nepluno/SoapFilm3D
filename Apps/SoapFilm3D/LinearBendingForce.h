#ifndef __LINEAR_BENDING_FORCE_H__
#define __LINEAR_BENDING_FORCE_H__

#include <Eigen/Core>
#include "Force.h"
#include <iostream>

class LinearBendingForce : public Force
{
public:
  
  LinearBendingForce( int idx1, int idx2, int idx3, const scalar& alpha, const scalar& beta, const Vector2s& theta0, const scalar& eb1n, const scalar& eb2n );
  
  virtual ~LinearBendingForce();
  
  virtual void addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E );
  
  virtual void addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE );
  
  virtual void addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE );
  
  virtual void addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE );
  
  virtual Force* createNewCopy();
  
  virtual void preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt );
private:
  
  int m_idx1;
  int m_idx2;
  int m_idx3;
  
  Matrix3s m_R;       // rotation matrix
  
  scalar m_alpha;     // stiffness coefficient
  scalar m_beta;      // damping coefficient
  Vector2s m_theta0;    // rest angle
  scalar m_eb1n;      // norm of e1 bar
  scalar m_eb2n;      // norm of e2 bar
  
  Vector3s m_x1;
  Vector3s m_x2;
  Vector3s m_x3;
  
  Vector3s m_L0;
  
  Vector3s m_RL0;
  
  scalar m_c1;
  scalar m_c2;
};

#endif
