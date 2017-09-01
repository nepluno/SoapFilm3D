#ifndef __SCENE_STEPPER__
#define __SCENE_STEPPER__

#include "MathDefs.h"

class VS3D;

class SceneStepper
{
public:
  virtual ~SceneStepper();
  
  virtual bool stepScene( VS3D& scene, scalar dt ) = 0;
  
  virtual std::string getName() const = 0;
};

#endif
