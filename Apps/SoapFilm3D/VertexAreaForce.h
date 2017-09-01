#ifndef __VERTEX_AREA_FORCE_H__
#define __VERTEX_AREA_FORCE_H__

#include <Eigen/Core>
#include "Force.h"
#include <iostream>

class VS3D;

class VertexAreaForce : public Force
{
public:

  VertexAreaForce(VS3D* parent, const scalar& k);
    
  virtual ~VertexAreaForce();
  
  virtual void addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E );
  
  virtual void addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE );
  
  virtual void addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE );
  
  virtual void addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE );
  
  virtual Force* createNewCopy();
private:
  VS3D* m_parent;
    scalar m_stiffness;
};

#endif
