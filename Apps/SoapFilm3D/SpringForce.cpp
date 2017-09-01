#include "SpringForce.h"

SpringForce::SpringForce( const std::pair<int,int>& endpoints, const scalar& k, const scalar& l0, const scalar& b )
: Force()
, m_endpoints(endpoints)
, m_k(k)
, m_l0(l0)
, m_b(b)
{
  assert( m_endpoints.first >= 0 );
  assert( m_endpoints.second >= 0 );
  assert( m_endpoints.first != m_endpoints.second );
  assert( m_k >= 0.0 );
  assert( m_l0 >= 0.0 );
  assert( m_b >= 0.0 );
}

SpringForce::~SpringForce()
{}
  
void SpringForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size()%3 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/3 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/3 );

  scalar l = (x.segment<3>(3*m_endpoints.second)-x.segment<3>(3*m_endpoints.first)).norm();
  E += 0.5*m_k*(l-m_l0)*(l-m_l0);
}

void SpringForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%3 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/3 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/3 );

  // Compute the elastic component
  Vector3s nhat = x.segment<3>(3*m_endpoints.second)-x.segment<3>(3*m_endpoints.first); 
  scalar l = nhat.norm(); 
  assert( l != 0.0 ); 
  nhat /= l;
  Vector3s fdamp = nhat;
  nhat *= m_k*(l-m_l0);
  gradE.segment<3>(3*m_endpoints.first)  -= nhat;
  gradE.segment<3>(3*m_endpoints.second) += nhat;

  // Compute the internal damping
  // Remember we are computing minus the force here
  fdamp *= m_b*fdamp.dot(v.segment<3>(3*m_endpoints.second)-v.segment<3>(3*m_endpoints.first));
  gradE.segment<3>(3*m_endpoints.first)  -= fdamp;
  gradE.segment<3>(3*m_endpoints.second) += fdamp;
}

void SpringForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size()%3 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/3 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/3 );

  // Contribution from elastic component
  Vector3s nhat = x.segment<3>(3*m_endpoints.second)-x.segment<3>(3*m_endpoints.first); 
  scalar l = nhat.norm();
  assert( l != 0 );
  nhat /= l;

  Matrix3s hess;
  hess = nhat*nhat.transpose();
  hess += (l-m_l0)*(Matrix3s::Identity()-hess)/l;
  hess *= m_k;
  
  // Contribution from damping
  Vector3s dv = v.segment<3>(3*m_endpoints.second)-v.segment<3>(3*m_endpoints.first);
  Matrix3s hessB = nhat*dv.transpose();
  hessB.diagonal().array() += nhat.dot(dv);
  hessB -= hessB*(nhat*nhat.transpose());
  hessB *= -m_b/l;
  
  hess -= hessB;
  
  for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_endpoints.first + r, 3 * m_endpoints.first + s, hess(r, s)));
  for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_endpoints.second + r, 3 * m_endpoints.second + s, hess(r, s)));
  for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_endpoints.first + r, 3 * m_endpoints.second + s, -hess(r, s)));
  for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_endpoints.second + r, 3 * m_endpoints.first + s, -hess(r, s)));
}

void SpringForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size()%3 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/3 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/3 );

  // Contribution from damping
  Vector3s nhat = x.segment<3>(3*m_endpoints.second)-x.segment<3>(3*m_endpoints.first); 
  scalar l = nhat.norm();
  assert( l != 0 );
  nhat /= l;
  
  Matrix3s hess;
  hess = m_b*nhat*nhat.transpose();
  
  for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_endpoints.first + r, 3 * m_endpoints.first + s, hess(r, s)));
  for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_endpoints.second + r, 3 * m_endpoints.second + s, hess(r, s)));
  for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_endpoints.first + r, 3 * m_endpoints.second + s, -hess(r, s)));
  for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_endpoints.second + r, 3 * m_endpoints.first + s, -hess(r, s)));
}

Force* SpringForce::createNewCopy()
{
  return new SpringForce(*this);
}
