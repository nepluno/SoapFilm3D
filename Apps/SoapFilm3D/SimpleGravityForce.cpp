#include "SimpleGravityForce.h"

SimpleGravityForce::SimpleGravityForce( const Vector3s& gravity )
: Force()
, m_gravity(gravity)
{
  assert( (m_gravity.array()==m_gravity.array()).all() );
  assert( (m_gravity.array()!=std::numeric_limits<scalar>::infinity()).all() );
}

SimpleGravityForce::~SimpleGravityForce()
{}
  
void SimpleGravityForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size()%3 == 0 );

  // Assume 0 potential is at origin
  for( int i = 0; i < x.size()/3; ++i ) E -= m(3*i)*m_gravity.dot(x.segment<3>(3*i));
}

void SimpleGravityForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%3 == 0 );
  
  for( int i = 0; i < x.size()/3; ++i ) gradE.segment<3>(3*i) -= m(3*i)*m_gravity;
}

void SimpleGravityForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size()%3 == 0 );
  // Nothing to do.
}

void SimpleGravityForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size()%3 == 0 );
  // Nothing to do.
}

Force* SimpleGravityForce::createNewCopy()
{
  return new SimpleGravityForce(*this);
}
