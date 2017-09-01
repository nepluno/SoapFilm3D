#ifndef __FORCE_H__
#define __FORCE_H__

#include <Eigen/Core>

#include "MathDefs.h"

class Force
{
public:

  virtual ~Force();
  
  virtual void addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E ) = 0;
  
  virtual void addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE ) = 0;
  
  virtual void addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE ) = 0;

  virtual void addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE ) = 0;
  
  virtual Force* createNewCopy() = 0;
  
  virtual void preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt );
};

#endif
