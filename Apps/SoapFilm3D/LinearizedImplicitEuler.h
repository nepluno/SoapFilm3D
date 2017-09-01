#ifndef __LINEARIZED_IMPLICIT_EULER__
#define __LINEARIZED_IMPLICIT_EULER__

#include <Eigen/Dense>
#include <iostream>

#include "SceneStepper.h"

class LinearizedImplicitEuler : public SceneStepper
{
public:
  LinearizedImplicitEuler();
  
  virtual ~LinearizedImplicitEuler();
  
  virtual bool stepScene( VS3D& scene, scalar dt );
  
  virtual std::string getName() const;

private:
  void zeroFixedDoFs( const VS3D& scene, VectorXs& vec );
};

#endif
