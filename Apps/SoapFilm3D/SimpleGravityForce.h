#ifndef __SIMPLE_GRAVITY_FORCE_H__
#define __SIMPLE_GRAVITY_FORCE_H__

#include <Eigen/Core>
#include "Force.h"
#include <iostream>

class SimpleGravityForce : public Force
{
public:

  SimpleGravityForce( const Vector3s& gravity );

  virtual ~SimpleGravityForce();
  
  virtual void addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E );
  
  virtual void addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE );
  
  virtual void addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE );
  
  virtual void addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE );
  
  virtual Force* createNewCopy();
private:
  Vector3s m_gravity;
};

#endif
