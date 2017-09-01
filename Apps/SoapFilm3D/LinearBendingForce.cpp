#include "LinearBendingForce.h"

LinearBendingForce::LinearBendingForce( int idx1, int idx2, int idx3, const scalar& alpha, const scalar& beta, const Vector2s& theta0, const scalar& eb1n, const scalar& eb2n )
: Force()
, m_idx1(idx1)
, m_idx2(idx2)
, m_idx3(idx3)
, m_alpha(alpha)
, m_theta0(theta0)
, m_eb1n(eb1n)
, m_eb2n(eb2n)
, m_beta(beta)
{
    assert( idx1 >= 0 );
    assert( idx2 >= 0 );
    assert( idx3 >= 0 );
    assert( idx1 != idx2 );
    assert( idx1 != idx3 );
    assert( idx2 != idx3 );
    assert( m_alpha >= 0.0 );
    assert( m_eb1n >= 0.0 );
    assert( m_eb2n >= 0.0 );
    
    m_c1 = (m_eb1n + m_eb2n) / m_eb1n * 0.5;
    m_c2 = (m_eb1n + m_eb2n) / m_eb2n * 0.5;
    
    m_x2.setZero();
    
    m_x1 = Vector3s(0, -m_eb1n, 0);
    m_x3 = Vector3s(sin(m_theta0(0))*cos(m_theta0(1))*m_eb2n, cos(m_theta0(0))*cos(m_theta0(1))*m_eb2n, sin(m_theta0(1))*m_eb2n);
    
    Vector3s center = (m_x1 + m_x2 + m_x3) / 3.0;
    m_x1 -= center;
    m_x2 -= center;
    m_x3 -= center;
    
    m_L0 = m_x2 * (m_c1 + m_c2) - m_x1 * m_c1 - m_x3 * m_c2;
}

void LinearBendingForce::preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt )
{
    const Vector3s& q1 = x.segment<3>(m_idx1 * 3);
    const Vector3s& q2 = x.segment<3>(m_idx2 * 3);
    const Vector3s& q3 = x.segment<3>(m_idx3 * 3);
    
    // compute COM
    Vector3s com = (q1 + q2 + q3) / 3.0;
    
    // compute COV Matrix
    Matrix3s A = (
                  (q1 - com) * m_x1.transpose() +
                  (q2 - com) * m_x2.transpose() +
                  (q3 - com) * m_x3.transpose()) / 3.0;
    
    const Eigen::JacobiSVD< Matrix3s >& svd = A.jacobiSvd(Eigen::ComputeFullU|Eigen::ComputeFullV);
    m_R = svd.matrixU() * svd.matrixV().transpose();
    
    m_RL0 = m_R * m_L0;
}

LinearBendingForce::~LinearBendingForce()
{}

void LinearBendingForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
    assert( x.size() == v.size() );
    assert( x.size()%3 == 0 );
    assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/3 );
    assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/3 );
    assert( m_idx3 >= 0 );  assert( m_idx3 < x.size()/3 );
    
    Vector3s x1 = x.segment<3>(m_idx1 * 3);
    Vector3s x2 = x.segment<3>(m_idx2 * 3);
    Vector3s x3 = x.segment<3>(m_idx3 * 3);
    
    Vector3s v1 = v.segment<3>(m_idx1 * 3);
    Vector3s v2 = v.segment<3>(m_idx2 * 3);
    Vector3s v3 = v.segment<3>(m_idx3 * 3);
    
    Vector3s L = x2 * (m_c1 + m_c2) - x1 * m_c1 - x3 * m_c2;
    Vector3s dL = L - m_RL0;
    Vector3s Ldot = v2 * (m_c1 + m_c2) - v1 * m_c1 - v3 * m_c2;
    
    E += 0.5 * m_alpha * dL.squaredNorm() + 0.5 * m_beta * Ldot.squaredNorm();
}

void LinearBendingForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
    assert( x.size() == v.size() );
    assert( x.size() == gradE.size() );
    assert( x.size()%3 == 0 );
    assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/3 );
    assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/3 );
    assert( m_idx3 >= 0 );  assert( m_idx3 < x.size()/3 );
    
    
    Vector3s x1 = x.segment<3>(m_idx1 * 3);
    Vector3s x2 = x.segment<3>(m_idx2 * 3);
    Vector3s x3 = x.segment<3>(m_idx3 * 3);
    
    Vector3s L = x2 * (m_c1 + m_c2) - x1 * m_c1 - x3 * m_c2;
    Vector3s dL = L - m_RL0;
    
    gradE.segment<3>(3 * m_idx1) += dL * (-m_c1) * m_alpha;
    gradE.segment<3>(3 * m_idx2) += dL * (m_c1 + m_c2) * m_alpha;
    gradE.segment<3>(3 * m_idx3) += dL * (-m_c2) * m_alpha;
    
    gradE.segment<3>(3 * m_idx1) += m_c1 * m_c1 * v.segment<3>(m_idx1 * 3) * m_beta;
    gradE.segment<3>(3 * m_idx2) += (m_c1 + m_c2) * (m_c1 + m_c2) * v.segment<3>(m_idx2 * 3) * m_beta;
    gradE.segment<3>(3 * m_idx3) += m_c2 * m_c2 * v.segment<3>(m_idx3 * 3) * m_beta;
}

void LinearBendingForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size()%3 == 0 );
    assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/3 );
    assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/3 );
    assert( m_idx3 >= 0 );  assert( m_idx3 < x.size()/3 );
    
    Matrix3s J00 = m_c1 * m_c1 * Matrix3s::Identity() * m_alpha;
    Matrix3s J01 = -m_c1 * (m_c1 + m_c2) * Matrix3s::Identity() * m_alpha;
    Matrix3s J02 = m_c1 * m_c2 * Matrix3s::Identity() * m_alpha;
    Matrix3s J10 = (m_c1 + m_c2) * (-m_c1) * Matrix3s::Identity() * m_alpha;
    Matrix3s J11 = (m_c1 + m_c2) * (m_c1 + m_c2) * Matrix3s::Identity() * m_alpha;
    Matrix3s J12 = (m_c1 + m_c2) * (-m_c2) * Matrix3s::Identity() * m_alpha;
    Matrix3s J20 = m_c2 * m_c1 * Matrix3s::Identity() * m_alpha;
    Matrix3s J21 = -m_c2 * (m_c1 + m_c2) * Matrix3s::Identity() * m_alpha;
    Matrix3s J22 = m_c2 * m_c2 * Matrix3s::Identity() * m_alpha;
    
    for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_idx1 + r, 3 * m_idx1 + s, J00(r, s)));
    for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_idx1 + r, 3 * m_idx2 + s, J01(r, s)));
    for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_idx1 + r, 3 * m_idx3 + s, J02(r, s)));
    for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_idx2 + r, 3 * m_idx1 + s, J10(r, s)));
    for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_idx2 + r, 3 * m_idx2 + s, J11(r, s)));
    for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_idx2 + r, 3 * m_idx3 + s, J12(r, s)));
    for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_idx3 + r, 3 * m_idx1 + s, J20(r, s)));
    for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_idx3 + r, 3 * m_idx2 + s, J21(r, s)));
    for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_idx3 + r, 3 * m_idx3 + s, J22(r, s)));
}

void LinearBendingForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size()%3 == 0 );
    assert( m_idx1 >= 0 );  assert( m_idx1 < x.size()/3 );
    assert( m_idx2 >= 0 );  assert( m_idx2 < x.size()/3 );
    assert( m_idx3 >= 0 );  assert( m_idx3 < x.size()/3 );
    
    Matrix3s J00 = m_c1 * m_c1 * Matrix3s::Identity() * m_beta;
    Matrix3s J11 = (m_c1 + m_c2) * (m_c1 + m_c2) * Matrix3s::Identity() * m_beta;
    Matrix3s J22 = m_c2 * m_c2 * Matrix3s::Identity() * m_beta;
    
    for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_idx1 + r, 3 * m_idx1 + s, J00(r, s)));
    for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_idx2 + r, 3 * m_idx2 + s, J11(r, s)));
    for(int r = 0; r < 3; ++r) for(int s = 0; s < 3; ++s) hessE.push_back(Triplets(3 * m_idx3 + r, 3 * m_idx3 + s, J22(r, s)));
}

Force* LinearBendingForce::createNewCopy()
{
    return new LinearBendingForce(*this);
}
