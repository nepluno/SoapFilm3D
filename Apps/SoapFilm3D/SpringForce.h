#ifndef __SPRING_FORCE_H__
#define __SPRING_FORCE_H__

#include <Eigen/Core>
#include "Force.h"
#include <iostream>

class SpringForce : public Force
{
public:

  SpringForce( const std::pair<int,int>& endpoints, const scalar& k, const scalar& l0, const scalar& b = 0.0 );

  virtual ~SpringForce();
  
  virtual void addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E );
  
  virtual void addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE );
  
  virtual void addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE );
  
  virtual void addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE );
  
  virtual Force* createNewCopy();

private:
  std::pair<int,int> m_endpoints;
  scalar m_k;
  scalar m_l0;
  scalar m_b;
};

#endif
